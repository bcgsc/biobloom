/*
 * BloomMapClassifier.h
 *
 *  Created on: Mar 23, 2016
 *      Author: cjustin
 */

#ifndef BLOOMMAPCLASSIFIER_H_
#define BLOOMMAPCLASSIFIER_H_

#include <string>
#include "boost/shared_ptr.hpp"
#include "Common/Options.h"
#include "Options.h"
#include <google/dense_hash_map>
#include <google/dense_hash_set>
#include "Common/Dynamicofstream.h"
#include <vector>
#include "BioBloomClassifier.h"
#include "bloomfilter/MIBloomFilter.hpp"
#include "bloomfilter/RollingHashIterator.h"
#include "btl_bloomfilter/ntHashIterator.hpp"
#include <iostream>
#include <map>
#include <set>
#include <algorithm>
#include <functional>

using namespace std;

static const std::string base64_chars =
             "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

class BloomMapClassifier {
public:
	explicit BloomMapClassifier(const string &filterFile);
	void filter(const vector<string> &inputFiles);

	void filterPair(const string &file1, const string &file2);

	virtual ~BloomMapClassifier();
private:
	MIBloomFilter<ID> m_filter;
	vector<string> m_fullIDs;
	vector<boost::shared_ptr<vector<ID> > > m_colliIDs;
	size_t m_numRead;

	vector<size_t> m_countTable;
	vector<double> m_freqTable;
	google::dense_hash_map<string, ID> m_idToIndex;

	bool fexists(const string &filename) const {
		ifstream ifile(filename.c_str());
		return ifile.good();
	}

	//from https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c
	inline unsigned nChoosek(unsigned n, unsigned k) {
		if (k > n)
			return 0;
		if (k * 2 > n)
			k = n - k;
		if (k == 0)
			return 1;

		int result = n;
		for (unsigned i = 2; i <= k; ++i) {
			result *= (n - i + 1);
			result /= i;
		}
		return result;
	}

	/*
	 * Computes criteria used for judging a read consisting of:
	 * Position of matches
	 * Number of actually evaluated k-mers
	 * Return count of matching k-mers to set
	 */
	//TODO proof of concept not yet optimized
	//TODO turn into streaming algorithm to terminate early
	inline vector<unsigned> getMatchSignature(const string &seq,
			unsigned &evaluatedSeeds, vector<vector<ID> > &hitsPattern) {
		vector<unsigned> matchPos;
		matchPos.reserve(seq.size() - m_filter.getKmerSize());
		ntHashIterator itr(seq, m_filter.getHashNum(), m_filter.getKmerSize());
		while (itr != itr.end()) {
			vector<ID> results = m_filter.at(*itr);
			if (results.size() > 0) {
				matchPos.push_back(itr.pos());
				hitsPattern.push_back(results);
			}
			++itr;
			++evaluatedSeeds;
		}
		return matchPos;
	}

	/*
	 * Takes match signature and number of hits and computes:
	 * Chance of pattern due to random chance
	 * Expected number of false positives
	 * Removed false positives assuming sequence linearity
	 */
	//TODO proof of concept not yet optimized
	//TODO turn into streaming algorithm to terminate early
	inline vector<unsigned> calcProb(const vector<unsigned> &hitsVector,
			unsigned &rmCount, unsigned evaluatedSeeds, double &probFP, ID id,
			const vector<vector<ID> > &hitsPattern) {

		if (!m_filter.getSeedValues().empty()) {
			cerr << "SS not yet supported" << endl;
			exit(1);
		}

		if (hitsVector.size() == 0) {
			return vector<unsigned>();
		}

		//algorithm overview
		//1) flag cluster of overlapping segments
		//2) within cluster flag run that is the worst section (least part of id)
		//3) using the last start and end point seen, prune reads of that run that contradict
		//4) repeat from step 1 until all conflicting segments are removed

		//index of start of cluster (to hitsVector) and runlength
		vector<unsigned> startIndex;
		startIndex.reserve(hitsVector.size());
		vector<unsigned> starts;
		starts.reserve(hitsVector.size());
		vector<unsigned> counts;
		counts.reserve(hitsVector.size());

		unsigned prevPos = hitsVector[0];
		bool ovlpPresent = false;

		//construct vectors for evaluation
		for (unsigned i = 0; i < hitsVector.size(); ++i) {
			//is part of the same chain
			if ( hitsVector[i] == prevPos + 1) {
				++counts.back();
			}
			//part of new chain
			else {
				//does not overlap
				if ( hitsVector[i]  > prevPos + m_filter.getKmerSize()) {
					if (!ovlpPresent) {
						startIndex.pop_back();
						starts.pop_back();
						counts.pop_back();
					}
					ovlpPresent = false;
				}
				//new chain but overlaps with previous chain (aside from first case)
				else if (prevPos !=  hitsVector[i] ) {
					ovlpPresent = true;
				}

				startIndex.push_back(i);
				starts.push_back(hitsVector[i]);
				counts.push_back(1);
			}
			prevPos =  hitsVector[i];
		}
		if (!ovlpPresent) {
			startIndex.pop_back();
			counts.pop_back();
			starts.pop_back();
		}
		//for tracking changes
		vector<unsigned> countsOld = vector<unsigned>(counts);
		vector<unsigned> startOld = vector<unsigned>(starts);

		for(unsigned i = 0 ; i < startIndex.size(); ++i){
			startOld.push_back(hitsVector[i]);
		}

		if (startIndex.size()) {
			unsigned removedAmount = 1;
			while (removedAmount) {
				removedAmount = rmWorstConflict(startIndex, id, hitsPattern,
						starts, counts);
				rmCount += removedAmount;
			}
		}

		//reconsititute new vector (for debugging purposes)
		vector<unsigned> removedVect;

		for (unsigned i = 0; i < starts.size(); ++i) {
			if (counts[i] != countsOld[i]) {
				for (unsigned j = 0; j < countsOld[i]; ++j) {
					if (starts[i] > (startOld[i] + j)
							|| (starts[i] + counts[i]) <= startOld[i] + j) {
						removedVect.push_back(startOld[i] + j);
					}
				}
			}
		}

		//		//computer final probablity
		double bfFPR = m_filter.getFPR();
		//binomial
		probFP = nChoosek(evaluatedSeeds, hitsVector.size() - rmCount)
				* pow(bfFPR, hitsVector.size() - rmCount)
				* pow((1.0 - bfFPR),
						evaluatedSeeds - (hitsVector.size() - rmCount));
		return removedVect;

	}

	/*
	 * Takes match signature and number of hits and computes:
	 * Chance of pattern due to random chance
	 * Expected number of false positives
	 * Removed false positives assuming sequence linearity
	 */
	//TODO proof of concept not yet optimized
	//TODO turn into streaming algorithm to terminate early
	inline vector<unsigned> calcProb(const vector<unsigned> &hitsVector,
			unsigned &rmCount, unsigned evaluatedSeeds, double &probFP) {

		if (!m_filter.getSeedValues().empty()) {
			cerr << "SS not yet supported" << endl;
			exit(1);
		}

		if(hitsVector.size() == 0){
			return vector<unsigned>();
		}

		//algorithm overview
		//1) flag cluster of overlapping segments
		//2) within cluster flag run that is the worst
		//3) using the last start and end point seen, prune reads of that run that contradict
		//4) repeat from step 1 until all conflicting segments are removed

		//start, runlength
		vector<unsigned> start;
		start.reserve(hitsVector.size());
		vector<unsigned> counts;
		counts.reserve(hitsVector.size());

		vector<unsigned>::const_iterator itr = hitsVector.begin();
		unsigned prevPos = *itr;
		bool ovlpPresent = false;

		//construct vectors for evaluation
		for (; itr != hitsVector.end(); itr++) {
			//is part of the same chain
			if (*itr == prevPos + 1) {
				++counts.back();
			}
			//part of new chain
			else {
				//does not overlap
				if (*itr > prevPos + m_filter.getKmerSize()) {
					if (!ovlpPresent) {
						start.pop_back();
						counts.pop_back();
					}
					ovlpPresent = false;
				}
				//new chain but overlaps with previous chain (aside from first case)
				else if (prevPos != *itr) {
					ovlpPresent = true;
				}

				start.push_back(*itr);
				counts.push_back(1);
			}
			prevPos = *itr;
		}
		if (!ovlpPresent) {
			start.pop_back();
			counts.pop_back();
		}

		//for tracking changes
		vector<unsigned> countsOld = vector<unsigned>(counts);
		vector<unsigned> startOld = vector<unsigned>(start);

		if (start.size()) {
			unsigned removedAmount = 1;
			while (removedAmount) {
//				for (unsigned i = 0; i < start.size(); ++i) {
//					if (counts[i] != 0) {
//						cerr << i << " " << start[i] << " " << counts[i]
//								<< endl;
//					}
//				}
				removedAmount = rmWorstConflict(start, counts);
				rmCount += removedAmount;
//				cerr << removedAmount;
//				cerr << "~\n";
			}
		}

		//reconsititute new vector (for debugging purposes)
		vector<unsigned> removedVect;

		for (unsigned i = 0; i < start.size(); ++i) {
			if (counts[i] != countsOld[i]) {
				for (unsigned j = 0; j < countsOld[i]; ++j) {
					if (start[i] > (startOld[i] + j)
							|| (start[i] + counts[i]) <= startOld[i] + j) {
						removedVect.push_back(startOld[i] + j);
					}
				}
			}
		}

//		//computer final probablity
		double bfFPR = m_filter.getFPR();
		//binomial
		probFP = nChoosek(evaluatedSeeds, hitsVector.size() - rmCount)
				* pow(bfFPR, hitsVector.size() - rmCount)
				* pow((1.0 - bfFPR),
						evaluatedSeeds - (hitsVector.size() - rmCount));
		return removedVect;

	}

	/*
	 * verbose output for debugging purposes
	 * Read Seq
	 * Seed match pattern
	 * Expected false positives
	 * Removed match pattern
	 * Number of removed matches
	 * Probability of false positive due to random chance
	 * Expected false positives
	 * Match pattern index
	 */
	inline void printVerbose(const string &header, const string &seq,
			const vector<unsigned> &hitsVector, unsigned evaluatedKmers,
			unsigned rmCount, const vector<unsigned> &rmMatch, double probFP,
			const vector< vector<ID> > &hitsPattern) {
		cerr << header << ' ' << evaluatedKmers << ' ' << rmCount
				<< ' ' << rmMatch.size() << ' '
				<< (evaluatedKmers * m_filter.getFPR()) << ' ' << probFP
				<< endl;
		cerr << seq << endl;
		cerr << vectToStr(hitsVector, seq) << endl;
		cerr << vectToStr(rmMatch, seq) << endl;
		cerr << vectToStr(hitsVector, hitsPattern, seq);
	}

	inline void printVerbose(const string &header, const string &seq) {
		unsigned evaluatedSeeds = 0;
		unsigned rmCount = 0;
		double probFP = 0.0;
		vector<vector<ID> > hitsPattern;
		vector<unsigned> sig = getMatchSignature(seq, evaluatedSeeds,
								hitsPattern);

		//compute counts
		google::dense_hash_map<ID, unsigned> counts;
		counts.set_empty_key(opt::EMPTY);
		for (vector<vector<ID> >::iterator i = hitsPattern.begin();
				i != hitsPattern.end(); i++) {
			for (vector<ID>::iterator j = i->begin(); j != i->end();
					j++) {
				google::dense_hash_map<ID, unsigned>::iterator tempItr =
						counts.find(*j);
				if (tempItr == counts.end()) {
					counts[*j] = 1;
				} else {
					++(tempItr->second);
				}
			}
		}

		//using the hitsPattern find highest ranking hits
		ID bestID = opt::EMPTY;
		unsigned bestCount = 0;

		for (google::dense_hash_map<ID, unsigned>::const_iterator itr = counts.begin();
				itr != counts.end(); itr++) {
			if(itr->second > bestCount){
				bestID = itr->first;
				bestCount = itr->second;
			}
		}
//
//		//compute match vector that prioritizes current best value
//		vector<unsigned> rmMatch = calcProb(sig, rmCount, evaluatedSeeds,
//				probFP, bestID, hitsPattern);
//		cerr << header << ' ' << evaluatedSeeds << ' ' << rmCount
//				<< ' ' << rmMatch.size() << ' '
//				<< (evaluatedSeeds * m_filter.getFPR()) << ' ' << probFP
//				<< endl;

		cout << header << ' ' << evaluatedSeeds << ' ' << rmCount << ' '
				<< (evaluatedSeeds * m_filter.getFPR()) << ' ' << probFP
				<< ' ' << base64_chars[bestID % 64] << ' ' << bestID << endl;
		cout << seq << endl;
		cout << vectToStr(sig, seq) << endl;
//		cerr << vectToStr(rmMatch, seq) << endl;
		cout << vectToStr(sig, hitsPattern, seq);
	}

	inline string vectToStr(const vector<unsigned> &hitsVector,
			const vector<vector<ID> > &hitsPattern, const string &seq) {
		stringstream ss;
		for (unsigned hVal = 0; hVal < m_filter.getHashNum(); ++hVal) {
			unsigned prevIndex = 1;
			unsigned index = 0;
			//print first stretch of zeros
			if (hitsVector.size()) {
				for (unsigned i = prevIndex; i <= hitsVector[index]; ++i) {
					ss << " ";
				}
				prevIndex = hitsVector[index] + 1;
				for (; index < hitsVector.size(); index++) {
					//print 0s
					for (unsigned i = prevIndex; i < hitsVector[index]; ++i) {
						ss <<  " ";
					}
					//print 1s
					ss << base64_chars[hitsPattern[index][hVal] % 64];
					prevIndex = hitsVector[index] + 1;
				}
				++prevIndex;
			}
			for (unsigned i = prevIndex;
					i <= (seq.size() - m_filter.getKmerSize() + 1); ++i) {
				ss << " ";
			}
			ss << "\n";
		}
		return ss.str();
	}


	inline string vectToStr(const vector<unsigned> &hitsVector,
			const string &seq) {
		stringstream ss;
		unsigned prevIndex = 1;
		vector<unsigned>::const_iterator itr = hitsVector.begin();
		//print first stretch of zeros
		if (itr != hitsVector.end()) {
			for (unsigned i = prevIndex; i <= *itr; ++i) {
				ss << 0;
			}
			prevIndex = *itr + 1;
			for (; itr != hitsVector.end(); itr++) {
				//print 0s
				for (unsigned i = prevIndex; i < *itr; ++i) {
					ss << 0;
				}
				//print 1s
				ss << 1;
				prevIndex = *itr + 1;
			}
			++prevIndex;
		}
		for (unsigned i = prevIndex; i <= (seq.size() - m_filter.getKmerSize() + 1);
				++i) {
			ss << 0;
		}
		return ss.str();
	}

	/*
	 * return count of elements in region that contain ID
	 */
	inline unsigned countID(ID id, const vector<vector<ID> > &hitsPattern,
			unsigned indexStart, unsigned indexEnd) {
		unsigned count = 0;
		for (unsigned i = indexStart; i < indexEnd; ++i) {
			bool noID = false;
			for (unsigned j = 0; j < m_filter.getKmerSize(); ++j) {
				if(hitsPattern[i][j] == id){
					noID = true;
					break;
				}
			}
			count += noID;
		}
		return count;
	}

	/*
	 * mutates worst conflict, shift the start and reducing the tail
	 * ie. Highest number of conflicts and that is very short
	 * Returns the amount removed
	 */
	inline unsigned rmWorstConflict(const vector<unsigned> &startIndex, ID id,
			const vector<vector<ID> > &hitsPattern,
			vector<unsigned> &startValues, vector<unsigned> &counts) {
		//index, conflict count
		vector<unsigned> conflicts(startIndex.size(), 0);
		vector<unsigned> rmLeft(startIndex.size(), 0);
		vector<unsigned> rmRight(startIndex.size(), 0);
		unsigned maxCount = 0;
		unsigned worstConfict = 0;
		//Find elements with most collisions external to its run
		for (unsigned i = 0; i < startIndex.size(); ++i) {
			if (counts[i] > 0) {
				for (unsigned j = i + 1; j < startIndex.size(); ++j) {
					if (counts[j]) {
						unsigned lastPos = startValues[i] + counts[i] - 1
								+ m_filter.getKmerSize();
						if (lastPos <= startValues[j]) {
							break;
						}
						//compute overlap
						unsigned overlap = lastPos - startValues[j];

						//if overlaps are contributed by bestID
						//k-mers overlapping
						conflicts[i] += countID(id, hitsPattern, startIndex[j],
								startIndex[j] + min(overlap, counts[j]));
						conflicts[j] += countID(id, hitsPattern, startIndex[i],
								startIndex[i] + min(overlap, counts[i]));
						rmRight[i] = min(overlap, counts[i]);
						rmLeft[j] = min(overlap, counts[j]);
						if ((conflicts[i] == maxCount
								&& counts[worstConfict] > counts[i])
								|| conflicts[i] > maxCount) {
							maxCount = conflicts[i];
							worstConfict = i;
						}
						if ((conflicts[j] == maxCount
								&& counts[worstConfict] > counts[j])
								|| conflicts[j] > maxCount) {
							maxCount = conflicts[j];
							worstConfict = j;
						}
					}
				}
			}
		}
		unsigned oldCount = counts[worstConfict];
		//shift based on conflict
		startValues[worstConfict] += rmLeft[worstConfict];
		//cut off conflict
		counts[worstConfict] -= rmRight[worstConfict];
		if (counts[worstConfict] < rmLeft[worstConfict]) {
			//completely obliterated
			counts[worstConfict] = 0;
			return oldCount;
		} else {
			counts[worstConfict] -= rmLeft[worstConfict];
			//			cerr << (rmLeft[worstConfict] + rmRight[worstConfict]) << endl;
			return rmLeft[worstConfict] + rmRight[worstConfict];
		}
	}

	/*
	 * mutates worst conflict, shift the start and reducing the tail
	 * ie. Highest number of conflicts and that is very short
	 * Returns the amount removed
	 */
	inline unsigned rmWorstConflict(vector<unsigned> &start,
			vector<unsigned> &counts) {
		//index, conflict count
		vector<unsigned> conflicts(start.size(), 0);
		vector<unsigned> rmLeft(start.size(), 0);
		vector<unsigned> rmRight(start.size(), 0);
		unsigned maxCount = 0;
		unsigned worstConfict = 0;
		//Find elements with most collisions external to its run
		for (unsigned i = 0; i < start.size(); ++i) {
			if (counts[i] > 0) {
				for (unsigned j = i + 1; j < start.size(); ++j) {
					if (counts[j]) {
						unsigned lastPos = start[i] + counts[i] - 1
								+ m_filter.getKmerSize();
						if (lastPos <= start[j]) {
							break;
						}
						//compute overlap
						unsigned overlap = lastPos - start[j];
						//k-mers overlapping
						conflicts[i] += min(overlap, counts[j]);
						conflicts[j] += min(overlap, counts[i]);
						rmRight[i] = min(overlap, counts[i]);
						rmLeft[j] = min(overlap, counts[j]);
						if ((conflicts[i] == maxCount
								&& counts[worstConfict] > counts[i])
								|| conflicts[i] > maxCount) {
							maxCount = conflicts[i];
							worstConfict = i;
						}
						if ((conflicts[j] == maxCount
								&& counts[worstConfict] > counts[j])
								|| conflicts[j] > maxCount) {
							maxCount = conflicts[j];
							worstConfict = j;
						}
					}
				}
			}
		}
		unsigned oldCount = counts[worstConfict];
		//shift based on conflict
		start[worstConfict] += rmLeft[worstConfict];
		//cut off conflict
		counts[worstConfict] -= rmRight[worstConfict];
		if (counts[worstConfict] < rmLeft[worstConfict]) {
			//completely obliterated
			counts[worstConfict] = 0;
			return oldCount;
		} else {
			counts[worstConfict] -= rmLeft[worstConfict];
//			cerr << (rmLeft[worstConfict] + rmRight[worstConfict]) << endl;
			return rmLeft[worstConfict] + rmRight[worstConfict];
		}
	}

	/*
	 * Returns the number of hits
	 */
	inline unsigned evaluateRead(const string &seq,
			google::dense_hash_map<ID, unsigned> &hitCounts) {
		unsigned nonZeroCount = 0;
		if (m_filter.getSeedValues().empty()) {
			if (m_filter.getType() == MIBloomFilter<ID>::MIBFMVAL) {
				for (ntHashIterator itr(seq, m_filter.getHashNum(),
						m_filter.getKmerSize()); itr != itr.end(); ++itr) {
					unsigned misses = 0;
					ID id = m_filter.at(*itr, opt::allowMisses, misses);
					if (id != 0) {
						if (id != opt::COLLI) {
							if (hitCounts.find(id) != hitCounts.end()) {
								++hitCounts[id];
							} else {
								hitCounts[id] = 1;
							}
						}
						++nonZeroCount;
					}
				}
			} else {
				for (ntHashIterator itr(seq, m_filter.getHashNum(),
						m_filter.getKmerSize()); itr != itr.end(); ++itr) {
					unsigned misses = 0;
					vector<ID> ids = m_filter.at(*itr, m_colliIDs,
							opt::allowMisses, misses);
					nonZeroCount += ids.size() > 0;
					for (unsigned i = 0; i < ids.size(); ++i) {
						ID id = ids[i];
						if (hitCounts.find(id) == hitCounts.end()) {
							hitCounts[id] = 0;
						}
						hitCounts[id] += ids.size() > 0;
					}
				}
			}
		} else {
			RollingHashIterator itr(seq, m_filter.getKmerSize(),
					m_filter.getSeedValues());
			if (m_filter.getType() == MIBloomFilter<ID>::MIBFMVAL) {
				while (itr != itr.end()) {
					ID id = m_filter.atBest(*itr, opt::allowMisses);
					if (id != 0) {
						if (id != opt::COLLI) {
							if (hitCounts.find(id) != hitCounts.end()) {
								++hitCounts[id];
							} else {
								hitCounts[id] = 1;
							}
						}
						++nonZeroCount;
					}
					++itr;
				}
			} else {
				while (itr != itr.end()) {
					unsigned misses = 0;
					vector<ID> ids = m_filter.at(*itr, m_colliIDs,
							opt::allowMisses, misses);
					nonZeroCount += misses == 0;
					nonZeroCount += ids.size() > 0;
					for (unsigned i = 0; i < ids.size(); ++i) {
						ID id = ids[i];
						if (hitCounts.find(id) == hitCounts.end()) {
							hitCounts[id] = 0;
						}
						++hitCounts[id];
//						hitCounts[id] += misses == 0;
					}
					++itr;
				}
//				nonZeroCount /= 2;
			}
		}
		return nonZeroCount;
	}

	/*
	 * Returns a vector of hits to a specific ID
	 */
	void convertToHits(const google::dense_hash_map<ID, unsigned> &hitCounts,
			vector<ID> &hits) {
		unsigned bestHit = 0;
		for (google::dense_hash_map<ID, unsigned>::const_iterator i =
				hitCounts.begin(); i != hitCounts.end(); ++i) {
			if (bestHit < i->second) {
				bestHit = i->second;
			}
		}
		for (google::dense_hash_map<ID, unsigned>::const_iterator i =
				hitCounts.begin(); i != hitCounts.end(); ++i) {
			if (bestHit <= i->second + opt::delta) {
				hits.push_back(i->first);
			}
		}
	}

	/*
	 * Returns a vector of best hits to a specific read
	 * Both reads must match
	 */
	unsigned convertToHitsBoth(
			const google::dense_hash_map<ID, unsigned> &hitCounts1,
			const google::dense_hash_map<ID, unsigned> &hitCounts2,
			vector<ID> &hits) {
		unsigned bestHit = 0;
		if (m_filter.getType() == MIBloomFilter<ID>::MIBFMVAL) {
			ID best = opt::EMPTY;
			//values shared between both reads
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				const google::dense_hash_map<ID, unsigned>::const_iterator &j =
						hitCounts2.find(i->first);
				if (j != hitCounts2.end()
						&& bestHit <= (i->second + j->second)) {
					if (bestHit == (i->second + j->second) && best < i->first) {
						continue;
					}
					bestHit = i->second + j->second;
					best = i->first;
				}
			}
			//if ensure pair makes any sense
			if (best != 0) {
				hits.push_back(best);
			}
		} else {
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				unsigned hitCount = i->second;
				google::dense_hash_map<ID, unsigned>::const_iterator itr =
						hitCounts2.find(i->first);
				if (itr != hitCounts2.end()) {
					hitCount += itr->second;
				} else {
					continue;
				}
				if (bestHit < hitCount) {
					bestHit = hitCount;
				}
			}
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				unsigned hitCount = i->second;
				google::dense_hash_map<ID, unsigned>::const_iterator itr =
						hitCounts2.find(i->first);
				if (itr != hitCounts2.end()) {
					hitCount += itr->second;
				} else {
					continue;
				}
				if (bestHit == hitCount) {
					hits.push_back(i->first);
				}
			}
		}
		return bestHit;
	}

	/*
	 * Returns a vector of best hits to a specific read
	 * Only one read needs to match
	 * outputs additional score based on results
	 */
	unsigned convertToHitsOnlyOne(
			const google::dense_hash_map<ID, unsigned> &hitCounts1,
			const google::dense_hash_map<ID, unsigned> &hitCounts2,
			vector<ID> &hits) {
		unsigned bestHit = 0;
		if (m_filter.getType() == MIBloomFilter<ID>::MIBFMVAL) {
			ID best = opt::EMPTY;
			//values shared between both reads
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts2.begin(); i != hitCounts2.end(); ++i) {
				if (bestHit <= i->second) {
					if (bestHit == i->second && best < i->first) {
						continue;
					}
					bestHit = i->second;
					best = i->first;
				}
			}
			//values found in only 1 read try to pick stronger one
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				const google::dense_hash_map<ID, unsigned>::const_iterator &j =
						hitCounts2.find(i->first);
				if (j != hitCounts2.end()) {
					if (bestHit <= (i->second + j->second)) {
						if (bestHit == (i->second + j->second)
								&& best < i->first) {
							continue;
						}
						bestHit = i->second + j->second;
						best = i->first;
					}
				} else {
					if (bestHit <= i->second) {
						if (bestHit == i->second && best < i->first) {
							continue;
						}
						bestHit = i->second;
						best = i->first;
					}
				}
			}
			hits.push_back(best);
		} else {

			//for tracking already found element in first set
			google::dense_hash_set<ID> tempSet;
			tempSet.set_empty_key(opt::EMPTY);

			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				unsigned hitCount = i->second;
				google::dense_hash_map<ID, unsigned>::const_iterator itr =
						hitCounts2.find(i->first);
				if (itr != hitCounts2.end()) {
					tempSet.insert(i->first);
					hitCount += itr->second;
				}
				if (bestHit < hitCount) {
					bestHit = hitCount;
				}
			}
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts2.begin(); i != hitCounts2.end(); ++i) {
				if (bestHit < i->second) {
					bestHit = i->second;
				}
			}

			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				unsigned hitCount = i->second;
				google::dense_hash_map<ID, unsigned>::const_iterator itr =
						hitCounts2.find(i->first);
				if (itr != hitCounts2.end()) {
					hitCount += itr->second;
				}
				if (bestHit <= hitCount + opt::delta) {
					hits.push_back(i->first);
				}
			}
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts2.begin(); i != hitCounts2.end(); ++i) {
				if (tempSet.find(i->first) == tempSet.end()) {
					if (bestHit <= i->second + opt::delta) {
						hits.push_back(i->first);
					}
				}
			}
		}
		return bestHit;
	}
};

#endif /* BLOOMMAPCLASSIFIER_H_ */
