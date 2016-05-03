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
#include "Common/Dynamicofstream.h"
#include <vector>
//#include "DataLayer/FastaReader.h"
#include "BioBloomClassifier.h"
#include "bloomfilter/BloomMapSSBitVec.hpp"
#include "bloomfilter/RollingHashIterator.h"

using namespace std;

class BloomMapClassifier {
public:
	explicit BloomMapClassifier(const string &filterFile);
	void filter(const vector<string> &inputFiles);

	void filterPair(const string &file1, const string &file2);

	virtual ~BloomMapClassifier();
private:
	BloomMapSSBitVec<ID> m_filter;
	vector<string> m_fullIDs;

//	//TODO: REFACTOR WITH BioBloomClassifier
//	inline void printSingleToFile(const string &outputFileName,
//			const kseq_t *rec,
//			google::dense_hash_map<string, boost::shared_ptr<Dynamicofstream> > &outputFiles) {
//		if (opt::outputType == "fa") {
//#pragma omp critical(outputFiles)
//			{
//				(*outputFiles[outputFileName]) << ">" << rec->name.s << "\n"
//						<< rec->seq.s << "\n";
//			}
//		} else {
//#pragma omp critical(outputFiles)
//			{
//				(*outputFiles[outputFileName]) << "@" << rec->name.s << "\n"
//						<< rec->seq.s << "\n+\n" << rec->qual.s << "\n";
//			}
//		}
//	}

	inline void printToTSV(const ID outputID, const string &name,
			Dynamicofstream &outputFile) {
		if (outputID != opt::EMPTY) {
#pragma omp critical(outputFiles)
			{
				outputFile << m_fullIDs[outputID] << "\t" << name << "\n";
			}
		}
	}

//	//TODO: REFACTOR WITH BioBloomClassifier
//	inline void printPairToFile(const string &outputFileName,
//			const kseq_t *rec1, const kseq_t *rec2,
//			google::dense_hash_map<string, boost::shared_ptr<Dynamicofstream> > &outputFiles) {
//		if (opt::outputType == "fa") {
//#pragma omp critical(outputFiles)
//			{
//				(*outputFiles[outputFileName + "_1"]) << ">" << rec1->name.s
//						<< "\n" << rec1->seq.s << "\n";
//				(*outputFiles[outputFileName + "_2"]) << ">" << rec2->name.s
//						<< "\n" << rec2->seq.s << "\n";
//			}
//		} else {
//#pragma omp critical(outputFiles)
//			{
//				(*outputFiles[outputFileName + "_1"]) << "@" << rec1->name.s << "\n"
//						<< rec1->seq.s << "\n+\n" << rec1->qual.s << "\n";
//				(*outputFiles[outputFileName + "_2"]) << "@" << rec2->name.s << "\n"
//						<< rec2->seq.s << "\n+\n" << rec2->qual.s << "\n";
//			}
//		}
//	}

	/*
	 * Returns the number of hits
	 */
	inline unsigned evaluateRead(const string &seq,
			google::dense_hash_map<ID, unsigned> &hitCounts) {
		RollingHashIterator itr(seq, m_filter.getKmerSize(),
				m_filter.getSeedValues());
		unsigned nonZeroCount = 0;

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
		return nonZeroCount;
	}

	/*
	 * Returns a vector of hits to a specific ID, favoring smaller IDs
	 */
	//TODO Currently only returns a single ID
	void convertToHits(const google::dense_hash_map<ID, unsigned> &hitCounts,
			vector<ID> &hits) {
		unsigned bestHit = 0;
		ID best = opt::EMPTY;
		for (google::dense_hash_map<ID, unsigned>::const_iterator i =
				hitCounts.begin(); i != hitCounts.end(); ++i) {
			if (bestHit <= i->second) {
				if(bestHit == i->second && best < i->first){
					continue;
				}
				bestHit = i->second;
				best = i->first;
			}
		}
		hits.push_back(best);
	}

	/*
	 * Returns a vector of hits to a specific read
	 * Both reads much match
	 */
	void convertToHitsBoth(
			const google::dense_hash_map<ID, unsigned> &hitCounts1,
			const google::dense_hash_map<ID, unsigned> &hitCounts2,
			vector<ID> &hits) {
		unsigned bestHit = 0;
		ID best = opt::EMPTY;
		for (google::dense_hash_map<ID, unsigned>::const_iterator i =
				hitCounts1.begin(); i != hitCounts1.end(); ++i) {
			const google::dense_hash_map<ID, unsigned>::const_iterator &j =
					hitCounts2.find(i->first);
			if (j != hitCounts2.end() && bestHit <= (i->second + j->second)) {
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
	}

	/*
	 * Returns a vector of hits to a specific read
	 * Only one read needs to match
	 */
	void convertToHitsOnlyOne(
			const google::dense_hash_map<ID, unsigned> &hitCounts1,
			const google::dense_hash_map<ID, unsigned> &hitCounts2,
			vector<ID> &hits) {
		unsigned bestHit = 0;
		ID best = opt::EMPTY;
		for (google::dense_hash_map<ID, unsigned>::const_iterator i =
				hitCounts2.begin(); i != hitCounts2.end(); ++i) {
			if (bestHit <= i->second) {
				if(bestHit == i->second && best < i->first){
					continue;
				}
				bestHit = i->second;
				best = i->first;
			}
		}
		for (google::dense_hash_map<ID, unsigned>::const_iterator i =
				hitCounts1.begin(); i != hitCounts1.end(); ++i) {
			const google::dense_hash_map<ID, unsigned>::const_iterator &j =
					hitCounts2.find(i->first);
			if (j != hitCounts2.end()) {
				if (bestHit <= (i->second + j->second)) {
					if (bestHit == (i->second + j->second) && best < i->first) {
						continue;
					}
					bestHit = i->second + j->second;
					best = i->first;
				}
			}
			else{
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
	}
};

#endif /* BLOOMMAPCLASSIFIER_H_ */
