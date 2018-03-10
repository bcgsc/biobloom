/*
 * MIBFClassifier.hpp
 *
 *  Created on: Dec 20, 2017
 *      Author: cjustin
 */

#ifndef CHROMIUMMAP_MIBFCLASSIFIER_HPP_
#define CHROMIUMMAP_MIBFCLASSIFIER_HPP_

#include <string>
#include "Options.h"
#include <boost/shared_ptr.hpp>
#include <google/dense_hash_map>
#include <google/dense_hash_set>
#include <vector>
#include "bloomfilter/MIBloomFilter.hpp"
#include "bloomfilter/RollingHashIterator.h"
#include "btl_bloomfilter/ntHashIterator.hpp"
#include <iostream>
#include <boost/math/distributions/binomial.hpp>
#include <tuple>

#include <zlib.h>
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "Common/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;
using boost::math::binomial;

static const std::string base64_chars =
		"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

class MIBFClassifier {
public:
	MIBFClassifier(const string &filterFile) :
			m_filter(MIBloomFilter<ID>(filterFile)), m_numRead(0) {
		//load in ID file
		string idFile = (filterFile).substr(0, (filterFile).length() - 3)
				+ "_ids.txt";
		google::dense_hash_map<ID, string> ids;
		ids.set_empty_key(opt::EMPTY);
		if (!fexists(idFile)) {
			cerr << "Error: " << idFile.c_str()
					<< " File cannot be opened. A corresponding id file is needed."
					<< endl;
			exit(1);
		} else {
			cerr << "Loading ID file: " << idFile.c_str() << endl;
		}
		m_fullIDs.push_back("");
		ifstream idFH(idFile.c_str(), ifstream::in);
		m_idToIndex.set_empty_key(m_fullIDs.back());
		string line;
		getline(idFH, line);
		while (idFH.good()) {
			ID id;
			string name;
			stringstream converter(line);
			converter >> id;
			converter >> name;
			m_fullIDs.push_back(name);
			m_idToIndex[m_fullIDs.back()] = id;
			getline(idFH, line);
		}
		idFH.close();
		if (opt::allowMisses > 0 && m_filter.getSeedValues().size() == 0) {
			cerr << "Allowed miss (-a) should not be used with k-mers" << endl;
			exit(1);
		}

		m_perFrameProb = vector<double>(m_fullIDs.size());
		m_perFrameProbMulti = vector<double>(m_fullIDs.size());
		//needed for later statistical calculation
		MIBloomFilter<ID>::calcFrameProbs(m_filter, m_perFrameProb,
				m_perFrameProbMulti);
	}

	void filterOld(const vector<string> &inputFiles) {
		string outputName = opt::outputPrefix + "_reads.tsv";

		//TODO output fasta?
		if (opt::outputType == "fq") {
			outputName = opt::outputPrefix + "_reads.fq";
		}

		Dynamicofstream readsOutput(outputName);

		//print out header info and initialize variables
		ResultsManager<ID> resSummary(m_fullIDs, false);

		cerr << "Filtering Start" << endl;

		//debugging
		size_t multiCount = 0;
		size_t incorrectCount = 0;

		for (vector<string>::const_iterator it = inputFiles.begin();
				it != inputFiles.end(); ++it) {
			gzFile fp;
			fp = gzopen(it->c_str(), "r");
			if (fp == NULL) {
				cerr << "Cannot open file" << it->c_str() << endl;
				exit(1);
			}
			kseq_t *seq = kseq_init(fp);
			int l;
#pragma omp parallel private(l)
			for (string sequence;;) {
				string name;
				string qual;
				string comment;
#pragma omp critical(sequence)
				{
					l = kseq_read(seq);
					sequence = string(seq->seq.s, seq->seq.l);
					name = string(seq->name.s, seq->name.l);
					qual = string(seq->qual.s, seq->qual.l);
					comment = string(seq->comment.s, seq->comment.l);
				}
				if (l >= 0) {
#pragma omp critical(totalReads)
					if (++m_numRead % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << m_numRead
								<< endl;
					}
					vector<vector<ID> > hitsPattern;
					vector<unsigned> saturation;
					unsigned evaluatedSeeds = 0;
					vector<unsigned> sig = getMatchSignature(sequence,
							evaluatedSeeds, hitsPattern, saturation);
					double pVal = 1.0;
					vector<double> pVals;
					unsigned maxCount = 0;
					vector<ID> signifResults;
					vector<unsigned> signifValues;
					vector<unsigned> fullSignifCounts;
					ID idIndex = evalRead(hitsPattern, evaluatedSeeds, pVal,
							maxCount, signifResults, signifValues,
							fullSignifCounts, pVals);
					resSummary.updateSummaryData(idIndex);

					const string &fullID =
							idIndex == opt::EMPTY ? NO_MATCH :
							idIndex == resSummary.getMultiMatchIndex() ?
									MULTI_MATCH : m_fullIDs.at(idIndex);
					//debugging
#pragma omp critical(cout)
				if (name == "noMatch") {
					cout << pVal << "\t" << "F" << "\t" << name << "\t"
							<< fullID << endl;
				} else {
					cout << pVal << "\t" << "T" << "\t" << name << "\t"
							<< fullID << endl;
				}
					bool match = false;
					unsigned tempCount = 0;
					for (unsigned i = 0; i < signifResults.size(); ++i) {
						//TODO resolve precision error better
						if (fullSignifCounts[i] + 3 >= maxCount) {
							if (name == m_fullIDs.at(signifResults[i])) {
								match = true;
							}
							++tempCount;
						}
					}
					if (tempCount > 1) {
#pragma omp atomic
						++multiCount;
					}

//#pragma omp critical(cout)
					if (!match) {
#pragma omp atomic
						++incorrectCount;
//						cout << m_numRead << "\tCorrectID:" << m_idToIndex[name]
//								<< "\tCorrectName:" << name << "\tPredictedID:"
//								<< m_idToIndex[fullID] << "\tPredictedName:"
//								<< fullID << "\tCorrectID:"
//								<< base64_chars[m_idToIndex[name] % 64]
//								<< "\tpVal:" << log10(pVal) * (-10.0) << "\t"
//								<< endl;
//
//						for (unsigned i = 0; i < signifResults.size(); ++i) {
//							cout << signifResults[i] << ","
//									<< m_fullIDs[signifResults[i]] << ","
//									<< base64_chars[signifResults[i] % 64]
//									<< "," << signifValues[i] << ","
//									<< fullSignifCounts[i] << "\t";
//						}
//						cout << endl;
//						printVerbose(name, comment, sequence, hitsPattern, sig,
//								saturation, idIndex);
					}

					if (idIndex != opt::EMPTY) {
						if (opt::outputType == "fq") {
#pragma omp critical(outputFiles)
							{
								readsOutput << "@" << name << " " << fullID
										<< "\n" << sequence << "\n+\n" << qual
										<< "\n";
							}
						} else {
#pragma omp critical(outputFiles)
							{
								readsOutput << fullID << "\t" << name << "\n";
							}
						}
					}
				} else {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}

		cerr << "Multiple Map Count:" << multiCount << endl;
		cerr << incorrectCount << endl;

		cerr << "Total Reads:" << m_numRead << endl;
		cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv"
				<< endl;

		Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
		summaryOutput << resSummary.getResultsSummary(m_numRead);
		summaryOutput.close();
		cout.flush();
	}

	void filterPairOld(const string &file1, const string &file2) {
		string outputName = opt::outputPrefix + "_reads.tsv";

		//TODO output fasta?
		if (opt::outputType == "fq") {
			outputName = opt::outputPrefix + "_reads.fq";
		}

		Dynamicofstream readsOutput(outputName);

		//print out header info and initialize variables
		ResultsManager<ID> resSummary(m_fullIDs, false);

		if (!m_filter.getSeedValues().empty()) {
			cerr << "Spaced Seed Mode" << endl;
		}

		cerr << "Filtering Start" << endl;

		gzFile fp1;
		gzFile fp2;
		fp1 = gzopen(file1.c_str(), "r");
		if (fp1 == NULL) {
			cerr << "Cannot open file " << file1.c_str() << endl;
			exit(1);
		}
		fp2 = gzopen(file2.c_str(), "r");
		if (fp2 == NULL) {
			cerr << "Cannot open file " << file2.c_str() << endl;
			exit(1);
		}
		kseq_t *seq1 = kseq_init(fp1);
		kseq_t *seq2 = kseq_init(fp2);
		int l1;
		int l2;
#pragma omp parallel private(l1, l2)
		for (string sequence1;;) {
			string sequence2;
			//TODO sanity check if names to make sure both names are the same
			string name1;
			string name2;
			string qual1;
			string qual2;
#pragma omp critical(sequence)
			{
				l1 = kseq_read(seq1);
				sequence1 = string(seq1->seq.s, seq1->seq.l);
				name1 = string(seq1->name.s, seq1->name.l);
				qual1 = string(seq1->qual.s, seq1->qual.l);

				l2 = kseq_read(seq2);
				sequence2 = string(seq2->seq.s, seq2->seq.l);
				name2 = string(seq2->name.s, seq2->name.l);
				qual2 = string(seq2->qual.s, seq2->qual.l);
			}
			if (l1 >= 0 && l2 >= 0) {
#pragma omp critical(totalReads)
				if (++m_numRead % opt::fileInterval == 0) {
					cerr << "Currently Reading Read Number: " << m_numRead
							<< endl;
				}

				google::dense_hash_map<ID, unsigned> hitCounts1;
				google::dense_hash_map<ID, unsigned> hitCounts2;
				hitCounts1.set_empty_key(opt::EMPTY);
				hitCounts2.set_empty_key(opt::EMPTY);

				if (m_numRead > opt::fileInterval) {
					exit(1);
				}

				unsigned bestHit = 0;

				vector<ID> hits;
				ID idIndex = opt::EMPTY;

				idIndex = resSummary.updateSummaryData(hits);
				if (idIndex != opt::EMPTY) {
					const string &fullID =
							idIndex == opt::EMPTY ? NO_MATCH :
							idIndex == resSummary.getMultiMatchIndex() ?
									MULTI_MATCH : m_fullIDs.at(idIndex);
					if (opt::outputType == "fq") {
#pragma omp critical(outputFiles)
						{
							readsOutput << "@" << name1 << " " << fullID << "\n"
									<< sequence1 << "\n+\n" << qual1 << "\n"
									<< "@" << name2 << " " << fullID << "\n"
									<< sequence2 << "\n+\n" << qual2 << "\n";
						}
					} else {
#pragma omp critical(outputFiles)
						{
							readsOutput << fullID << "\t" << name1 << "\t"
									<< bestHit << "\n";
						}
					}
				}
			} else {
				if (l1 != -1 || l2 != -1) {
					cerr << "Terminated without getting to eof at read "
							<< m_numRead << endl;
					exit(1);
				}
				break;
			}
		}
		kseq_destroy(seq1);
		kseq_destroy(seq2);
		gzclose(fp1);
		gzclose(fp2);

		cerr << "Total Reads:" << m_numRead << endl;
		cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv"
				<< endl;

		Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
		summaryOutput << resSummary.getResultsSummary(m_numRead);
		summaryOutput.close();
		cout.flush();
	}

	void filter(const vector<string> &inputFiles) {
		string outputName = opt::outputPrefix + "_reads.tsv";

		//TODO output fasta?
		if (opt::outputType == "fq") {
			outputName = opt::outputPrefix + "_reads.fq";
		}

		Dynamicofstream readsOutput(outputName);

		//print out header info and initialize variables
		ResultsManager<ID> resSummary(m_fullIDs, false);

		cerr << "Filtering Start" << endl;

		for (vector<string>::const_iterator it = inputFiles.begin();
				it != inputFiles.end(); ++it) {
			gzFile fp;
			fp = gzopen(it->c_str(), "r");
			if (fp == NULL) {
				cerr << "Cannot open file" << it->c_str() << endl;
				exit(1);
			}
			kseq_t *seq = kseq_init(fp);
			int l;
#pragma omp parallel private(l)
			for (string sequence;;) {
				string name;
				string qual;
				string comment;
#pragma omp critical(sequence)
				{
					l = kseq_read(seq);
					sequence = string(seq->seq.s, seq->seq.l);
					name = string(seq->name.s, seq->name.l);
					qual = string(seq->qual.s, seq->qual.l);
					comment = string(seq->comment.s, seq->comment.l);
				}
				if (l >= 0) {
#pragma omp critical(totalReads)
					if (++m_numRead % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << m_numRead
								<< endl;
					}

					vector<ID> signifResults = classify(sequence);

					if(signifResults.size() == 0){
						resSummary.updateSummaryData(opt::EMPTY);
					}
					else if(signifResults.size() == 1){
						resSummary.updateSummaryData(signifResults[0]);
					}
					else{
						resSummary.updateSummaryData(resSummary.getMultiMatchIndex());
					}

#pragma omp critical(outputFiles)
					if (opt::outputType == "fq") {
						readsOutput << "@" << name << " ";
						for (unsigned i = 0; i < signifResults.size(); ++i) {
							readsOutput << " " << m_fullIDs[signifResults[i]];
						}
						readsOutput << "\n" << sequence << "\n+\n" << qual
								<< "\n";
					} else {
						for (unsigned i = 0; i < signifResults.size(); ++i) {
							readsOutput << m_fullIDs[signifResults[i]] << "\t"
									<< name << "\n";
						}
					}
				} else {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}

		cerr << "Total Reads:" << m_numRead << endl;
		cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv"
				<< endl;

		Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
		summaryOutput << resSummary.getResultsSummary(m_numRead);
		summaryOutput.close();
		cout.flush();
	}

private:
	MIBloomFilter<ID> m_filter;
	size_t m_numRead;
	vector<string> m_fullIDs;
	google::dense_hash_map<string, ID> m_idToIndex;
	vector<double> m_perFrameProb;
	vector<double> m_perFrameProbMulti;

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
	 * Streamlined classification code
	 */
	inline vector<ID> classify(const string &seq) {
		if (m_filter.getSeedValues().size() > 0) {
			RollingHashIterator itr(seq, m_filter.getKmerSize(),
					m_filter.getSeedValues());
			return MIBloomFilter<ID>::query(m_filter, itr, m_perFrameProb,
					m_perFrameProbMulti, opt::score, opt::score);

		} else {
			ntHashIterator itr(seq, m_filter.getHashNum(),
					m_filter.getKmerSize());
			return MIBloomFilter<ID>::query(m_filter, itr, m_perFrameProb,
					m_perFrameProbMulti, opt::score, opt::score);
		}
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
			unsigned &evaluatedSeeds, vector<vector<ID> > &hitsPattern,
			vector<unsigned> &saturation) {
		vector<unsigned> matchPos;
		matchPos.reserve(seq.size() - m_filter.getKmerSize());

		if (m_filter.getSeedValues().size() > 0) {
			RollingHashIterator itr(seq, m_filter.getKmerSize(),
					m_filter.getSeedValues());
			while (itr != itr.end()) {
				bool saturated = true;
				vector<ID> results = m_filter.at(*itr, saturated,
						opt::allowMisses);
				if (results.size() > 0) {
					matchPos.push_back(itr.pos());
					hitsPattern.push_back(results);
				}
				if (saturated) {
					saturation.push_back(itr.pos());
				}
				++itr;
				++evaluatedSeeds;
			}
		} else {
			ntHashIterator itr(seq, m_filter.getHashNum(),
					m_filter.getKmerSize());
			while (itr != itr.end()) {
				bool saturated = true;
				vector<ID> results = m_filter.at(*itr, saturated,
						opt::allowMisses);
				if (results.size() > 0) {
					matchPos.push_back(itr.pos());
					hitsPattern.push_back(results);
				}
				if (saturated) {
					saturation.push_back(itr.pos());
				}
				++itr;
				++evaluatedSeeds;
			}
		}
		return matchPos;
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
	inline void printVerbose(const string &header, const string &comment,
			const string &seq, const vector<unsigned> &hitsVector,
			unsigned evaluatedKmers, unsigned rmCount,
			const vector<unsigned> &rmMatch, double probFP,
			const vector<vector<ID> > &hitsPattern) {
		cerr << header << ' ' << comment << ' ' << evaluatedKmers << ' '
				<< rmCount << ' ' << rmMatch.size() << ' ' << probFP << endl;
		cerr << seq << endl;
		cerr << vectToStr(hitsVector, seq) << endl;
		cerr << vectToStr(rmMatch, seq) << endl;
		cerr << vectToStr(hitsVector, hitsPattern, seq);
	}

	inline ID evalRead(const vector<vector<ID> > &hitsPattern,
			unsigned evaluatedSeeds, double &pVal, unsigned &maxCount,
			vector<ID> &signifResults, vector<unsigned> &signifCounts,
			vector<unsigned> &fullSignifCounts, vector<double> &pVals) {
		google::dense_hash_map<ID, unsigned> counts;
		google::dense_hash_map<ID, unsigned> fullCounts;
		counts.set_empty_key(opt::EMPTY);
		fullCounts.set_empty_key(opt::EMPTY);
		for (vector<vector<ID> >::const_iterator i = hitsPattern.begin();
				i != hitsPattern.end(); i++) {
			//to determine if already added for this frame
			google::dense_hash_set<ID> tempIDs;
			tempIDs.set_empty_key(opt::EMPTY);
			unsigned count = 0;
			for (vector<ID>::const_iterator j = i->begin(); j != i->end();
					j++) {
				if (*j != opt::EMPTY) {
					if (tempIDs.find(*j) == tempIDs.end()) {
						google::dense_hash_map<ID, unsigned>::iterator tempItr =
								counts.find(*j);
						if (tempItr == counts.end()) {
							counts[*j] = 1;
						} else {
							++(tempItr->second);
						}
						tempIDs.insert(*j);
					}
					++count;
				}
			}
			if (count == m_filter.getHashNum()) {
				google::dense_hash_set<ID> tempIDsFull;
				tempIDsFull.set_empty_key(opt::EMPTY);
				for (vector<ID>::const_iterator j = i->begin(); j != i->end();
						j++) {
					if (tempIDsFull.find(*j) == tempIDsFull.end()) {
						google::dense_hash_map<ID, unsigned>::iterator tempItrFull =
								fullCounts.find(*j);
						if (tempItrFull == fullCounts.end()) {
							fullCounts[*j] = 1;
						} else {
							++(tempItrFull->second);
						}
						tempIDsFull.insert(*j);
					}
				}
			}
		}

		ID maxID = opt::EMPTY;
		double minVal = 1.0;
		double maxMinVal = 1.0;
		double adjustedPValThreshold = 1.0
				- pow(1.0 - opt::score, 1.0 / double(m_fullIDs.size() - 1));
		pVal = 1.0;
		for (google::dense_hash_map<ID, unsigned>::const_iterator itr =
				counts.begin(); itr != counts.end(); itr++) {
			//TODO use complement cdf? so I don't have to subtract?
			//TODO check maximum numerical precision on perFrameProb
			binomial bin(evaluatedSeeds, 1.0 - m_perFrameProb[itr->first]);
			double cumProb = cdf(bin, evaluatedSeeds - itr->second);
			if (adjustedPValThreshold > cumProb) {
				if (fullCounts[itr->first] > maxCount
						|| (fullCounts[itr->first] == maxCount
								&& maxMinVal > cumProb)) {
					maxCount = fullCounts[itr->first];
					maxID = itr->first;
					maxMinVal = cumProb;
				}
				fullSignifCounts.push_back(fullCounts[itr->first]);
				signifCounts.push_back(itr->second);
				signifResults.push_back(itr->first);
				pVals.push_back(log10(cumProb) * -10);
			}
			if (minVal > cumProb) {
				minVal = cumProb;
			}
		}
		//sidak method - leads to precision errors
//		pVal = 1.0 - pow(1.0-minVal, m_fullIDs.size() - 1);
		//bonferroni
		pVal = minVal * (m_fullIDs.size() - 1);
		pVal = pVal > 1 ? 1.0 : pVal;

		if (pVal <= opt::score) {
			return maxID;
		}
		return opt::EMPTY;
	}

	//TODO not optimized
	inline double calcProbSingleFrame(double freq) {
		double occupancy = double(m_filter.getPop()) / double(m_filter.size());
		double probTotal = 0.0;
		for (unsigned i = 0; i <= opt::allowMisses; i++) {
			double prob = nChoosek(m_filter.getHashNum(), i);
			prob *= pow(occupancy, m_filter.getHashNum() - i);
			prob *= pow(1.0 - occupancy, i);
			prob *= (1.0 - pow(1.0 - freq, m_filter.getHashNum() - i));
			probTotal += prob;
		}
		return probTotal;
	}

	inline unsigned getMinCount(unsigned length, double eventProb) {
		binomial bin(length, 1.0 - eventProb);
		double criticalScore = 1.0
				- pow(1.0 - opt::score, 1.0 / double(m_fullIDs.size() - 1));
		unsigned i = 0;
		for (; i < length; ++i) {
			double cumProb = cdf(bin, length - i);
			if (criticalScore >= cumProb) {
				break;
			}
		}
		return (i);
	}

	inline void printVerbose(const string &header, const string &comment,
			const string &seq, const vector<vector<ID> > &hitsPattern,
			const vector<unsigned> &sig, const vector<unsigned> &saturation,
			ID value) {
		unsigned evaluatedSeeds = 0;
		cout << header << ' ' << comment << ' ' << evaluatedSeeds << ' '
				<< base64_chars[value % 64] << ' ' << value << endl;
		cout << seq << endl;
		cout << vectToStr(sig, seq) << endl;
		cout << vectToStr(saturation, seq) << endl;
		cout << vectToStr(sig, hitsPattern, seq);
	}

	inline string vectToStr(const vector<unsigned> &hitsVector,
			const vector<vector<ID> > &hitsPattern, const string &seq) {
		stringstream ss;
		{
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
						ss << " ";
					}
					unsigned count = 0;
					for (unsigned hVal = 0; hVal < m_filter.getHashNum();
							++hVal) {
						if (hitsPattern[index][hVal] != opt::EMPTY) {
							++count;
						}
					}
					ss << count;
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
						ss << " ";
					}
					if (hitsPattern[index][hVal] == opt::EMPTY) {
						ss << " ";
					} else {
						ss << base64_chars[hitsPattern[index][hVal] % 64];
					}
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
		for (unsigned i = prevIndex;
				i <= (seq.size() - m_filter.getKmerSize() + 1); ++i) {
			ss << 0;
		}
		return ss.str();
	}

	unsigned getThreshold(double confidence, unsigned trials, double p) {
		unsigned threshold = 0;
		binomial quiz(trials, p);
		double probFP = pdf(quiz, threshold);
		while (probFP < confidence) {
			threshold++;
			probFP += pdf(quiz, threshold);
			cerr << probFP << " " << threshold << endl;
		}
		return threshold;
	}
};

#endif /* CHROMIUMMAP_MIBFCLASSIFIER_HPP_ */
