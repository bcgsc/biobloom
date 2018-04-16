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
#include <BioBloomClassifier.h>
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
		m_minCount.set_empty_key(0);
//		m_fprCount.set_empty_key(0);
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
					vector<double> multiMatchPVal;
					ID idIndex = evalRead(hitsPattern, evaluatedSeeds, pVal,
							maxCount, signifResults, signifValues,
							fullSignifCounts, pVals);
					resSummary.updateSummaryData(idIndex);

					const string &fullID =
							idIndex == opt::EMPTY ? NO_MATCH :
							idIndex == resSummary.getMultiMatchIndex() ?
									MULTI_MATCH : m_fullIDs.at(idIndex);
					//debugging
//#pragma omp critical(cout)
//					if (name == "noMatch") {
//						cout << pVal << "\t" << "F" << "\t" << name << "\t"
//								<< fullID << endl;
//					} else {
//						cout << pVal << "\t" << "T" << "\t" << name << "\t"
//								<< fullID << endl;
//					}
					bool match = false;
					unsigned tempCount = 0;
					for (unsigned i = 0; i < signifResults.size(); ++i) {
						//TODO resolve precision error better
						if (fullSignifCounts[i] >= maxCount) {
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

					if (!match) {
#pragma omp atomic
						++incorrectCount;
					}

#pragma omp critical(cout)
					if (!match || tempCount > 1) {
						cout << m_numRead << "\tCorrectID:" << m_idToIndex[name]
								<< "\tCorrectName:" << name << "\tPredictedID:"
								<< m_idToIndex[fullID] << "\tPredictedName:"
								<< fullID << "\tpVal:" << log10(pVal) * (-10.0)
								<< "\t" << endl;

						for (unsigned i = 0; i < signifResults.size(); ++i) {
							cout << signifResults[i] << ","
									<< m_fullIDs[signifResults[i]] << ","
									<< base64_chars[signifResults[i] % 64]
									<< "," << signifValues[i] << ","
									<< fullSignifCounts[i] << "\t";
						}
						cout << endl;
						printVerbose(name, comment, sequence, hitsPattern, sig,
								saturation, idIndex);
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

		cerr << "Multiple Map Count: " << multiCount << endl;
		cerr << "Incorrect: " << incorrectCount << endl;

		cerr << "Total Reads: " << m_numRead << endl;
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

		double startTime = omp_get_wtime();

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

					vector<pair<ID, double>> signifResults = classify2(sequence);

					if (signifResults.size() == 0) {
						resSummary.updateSummaryData(opt::EMPTY);
					} else if (signifResults.size() == 1) {
						resSummary.updateSummaryData(signifResults[0].first);
					} else {
						resSummary.updateSummaryData(
								resSummary.getMultiMatchIndex());
					}

#pragma omp critical(outputFiles)
					if (opt::outputType == "fq") {
						readsOutput << "@" << name << " ";
						for (unsigned i = 0; i < signifResults.size(); ++i) {
							readsOutput << " "
									<< m_fullIDs[signifResults[i].first];
						}
						readsOutput << "\n" << sequence << "\n+\n" << qual
								<< "\n";
					} else {
						if (signifResults.empty()) {
							readsOutput << "noMatch\t" << name << "\t"
									<< comment << "\t0\n";
						} else {
							for (unsigned i = 0; i < signifResults.size();
									++i) {
								readsOutput << m_fullIDs[signifResults[i].first]
										<< "\t" << name << "\t" << comment
										<< "\t" << signifResults.size() << "\n";
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

		cerr << "Classification time (s): " << (omp_get_wtime() - startTime)
						<< endl;

		cerr << "Total Reads:" << m_numRead << endl;
		cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv"
				<< endl;

		Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
		summaryOutput << resSummary.getResultsSummary(m_numRead);
		summaryOutput.close();
		cout.flush();
	}

	/*
	 * Filters reads -> uses paired end information
	 * Assumes only one hash signature exists (load only filters with same
	 * hash functions)
	 * prints reads
	 */
	void filterPairPrint(const string &file1, const string &file2) {

		gzFile fp1, fp2;
		fp1 = gzopen(file1.c_str(), "r");
		if (fp1 == NULL) {
			cerr << "file " << file1.c_str() << " cannot be opened" << endl;
			exit(1);
		}
		fp2 = gzopen(file2.c_str(), "r");
		if (fp2 == NULL) {
			cerr << "file " << file2.c_str() << " cannot be opened" << endl;
			exit(1);
		}
		kseq_t *kseq1 = kseq_init(fp1);
		kseq_t *kseq2 = kseq_init(fp2);
		FaRec rec1;
		FaRec rec2;

		//results summary object
		ResultsManager<unsigned> resSummary(m_fullIDs, false);

		size_t totalReads = 0;

		vector<Dynamicofstream*> outputFiles1(m_fullIDs.size(), 0);
		vector<Dynamicofstream*> outputFiles2(m_fullIDs.size(), 0);
		//initialize variables
		unsigned index = 0;

		for (vector<string>::const_iterator i = m_fullIDs.begin(); i != m_fullIDs.end(); ++i) {
			outputFiles1[index] = new Dynamicofstream(
					opt::outputPrefix + "/" + *i + "_1." + opt::outputType
							+ opt::filePostfix);
			outputFiles2[index++] = new Dynamicofstream(
					opt::outputPrefix + "/" + *i + "_2." + opt::outputType
							+ opt::filePostfix);
		}

		if (opt::verbose) {
			cerr << "Filtering Start" << "\n";
		}

		double startTime = omp_get_wtime();

#pragma omp parallel private(rec1, rec2)
		for (int l1, l2;;) {
#pragma omp critical(kseq)
			{
				l1 = kseq_read(kseq1);
				if (l1 >= 0) {
					rec1.seq = string(kseq1->seq.s, l1);
					rec1.header = string(kseq1->name.s, kseq1->name.l);
					rec1.qual = string(kseq1->qual.s, kseq1->qual.l);
					rec1.comment = string(kseq1->comment.s, kseq1->comment.l);
				}
				l2 = kseq_read(kseq2);
				if (l2 >= 0) {
					rec2.seq = string(kseq2->seq.s, l2);
					rec2.header = string(kseq2->name.s, kseq2->name.l);
					rec2.qual = string(kseq2->qual.s, kseq2->qual.l);
					rec2.comment = string(kseq1->comment.s, kseq1->comment.l);
				}
			}
			if (l1 >= 0 && l2 >= 0) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				vector<pair<ID, double>> signifResults = classifyPair(rec1.seq,
						rec2.seq);

				if (signifResults.size() == 0) {
					resSummary.updateSummaryData(opt::EMPTY);
				} else if (signifResults.size() == 1) {
					resSummary.updateSummaryData(signifResults[0].first);
					printPairToFile(signifResults[0].first, rec1, rec2,
							outputFiles1, outputFiles2);
				} else {
					resSummary.updateSummaryData(
							resSummary.getMultiMatchIndex());
					for (unsigned i = 0; i < signifResults.size(); ++i) {
						printPairToFile(signifResults[i].first, rec1, rec2,
								outputFiles1, outputFiles2);
					}
				}

			} else
				break;
		}
		cerr << "Classification time (s): " << (omp_get_wtime() - startTime)
				<< endl;
		kseq_destroy(kseq1);
		kseq_destroy(kseq2);

		//close sorting files
		for (unsigned i = 0; i < m_fullIDs.size(); ++i) {
			delete (outputFiles1[i]);
			if (opt::verbose) {
				cerr << "File written to: "
						<< opt::outputPrefix + "_" + m_fullIDs[i] + "_1."
								+ opt::outputType + opt::filePostfix << endl;
			}
			delete (outputFiles2[i]);
			if (opt::verbose) {
				cerr << "File written to: "
						<< opt::outputPrefix + "_" + m_fullIDs[i] + "_2."
								+ opt::outputType + opt::filePostfix << endl;
			}
		}

		Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
		summaryOutput << resSummary.getResultsSummary(totalReads);
		summaryOutput.close();
		cout.flush();
	}

	void filterPair(const string &file1, const string &file2) {

		gzFile fp1, fp2;
		fp1 = gzopen(file1.c_str(), "r");
		if (fp1 == NULL) {
			cerr << "file " << file1.c_str() << " cannot be opened" << endl;
			exit(1);
		}
		fp2 = gzopen(file2.c_str(), "r");
		if (fp2 == NULL) {
			cerr << "file " << file2.c_str() << " cannot be opened" << endl;
			exit(1);
		}
		kseq_t *kseq1 = kseq_init(fp1);
		kseq_t *kseq2 = kseq_init(fp2);
		FaRec rec1;
		FaRec rec2;

		//results summary object
		ResultsManager<unsigned> resSummary(m_fullIDs, false);

		size_t totalReads = 0;

		string outputName = opt::outputPrefix + "_reads.tsv";

		//TODO output fasta?
		if (opt::outputType == "fq") {
			outputName = opt::outputPrefix + "_reads.fq";
		}
		else if (opt::outputType == "fa") {
			outputName = opt::outputPrefix + "_reads.fa";
		}

		Dynamicofstream readsOutput(outputName + opt::filePostfix);

		if (opt::verbose) {
			cerr << "Filtering Start" << "\n";
		}
		double startTime = omp_get_wtime();

#pragma omp parallel private(rec1, rec2)
		for (int l1, l2;;) {
#pragma omp critical(kseq)
			{
				l1 = kseq_read(kseq1);
				if (l1 >= 0) {
					rec1.seq = string(kseq1->seq.s, l1);
					rec1.header = string(kseq1->name.s, kseq1->name.l);
					rec1.qual = string(kseq1->qual.s, kseq1->qual.l);
					rec1.comment = string(kseq1->comment.s, kseq1->comment.l);
				}
				l2 = kseq_read(kseq2);
				if (l2 >= 0) {
					rec2.seq = string(kseq2->seq.s, l2);
					rec2.header = string(kseq2->name.s, kseq2->name.l);
					rec2.qual = string(kseq2->qual.s, kseq2->qual.l);
					rec2.comment = string(kseq1->comment.s, kseq1->comment.l);
				}
			}
			if (l1 >= 0 && l2 >= 0) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				vector<pair<ID, double>> signifResults = classifyPair(rec1.seq,
						rec2.seq);

				if (signifResults.size() == 0) {
					resSummary.updateSummaryData(opt::EMPTY);
				} else if (signifResults.size() == 1) {
					resSummary.updateSummaryData(signifResults[0].first);
					printPairToFile(rec1, rec2, readsOutput);
				} else {
					resSummary.updateSummaryData(
							resSummary.getMultiMatchIndex());
					for (unsigned i = 0; i < signifResults.size(); ++i) {
						printPairToFile(rec1, rec2, readsOutput);
					}
				}
#pragma omp critical(outputFiles)
				if (opt::outputType == "") {
					if (signifResults.empty()) {
						readsOutput << "noMatch\t" << rec1.header << "\t"
								<< rec1.comment << "\t0\n";
					} else {
						for (unsigned i = 0; i < signifResults.size(); ++i) {
							readsOutput << m_fullIDs[signifResults[i].first]
									<< "\t" << rec1.header << "\t"
									<< rec1.comment << "\t"
									<< signifResults.size() << "\n";
						}
					}
				}
			} else
				break;
		}
		cerr << "Classification time (s): " << (omp_get_wtime() - startTime)
				<< endl;
		kseq_destroy(kseq1);
		kseq_destroy(kseq2);
		readsOutput.close();

		Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
		summaryOutput << resSummary.getResultsSummary(totalReads);
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
	google::dense_hash_map<unsigned, boost::shared_ptr<vector<unsigned>>> m_minCount;
//	google::dense_hash_map<unsigned, boost::shared_ptr<vector<unsigned>>> m_fprCount;

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
	inline vector<pair<ID, double>> classify(const string &seq) {
		if (m_filter.getSeedValues().size() > 0) {
			RollingHashIterator itr(seq, m_filter.getKmerSize(),
					m_filter.getSeedValues());
			return m_filter.query(itr, m_perFrameProb, m_perFrameProbMulti,
					opt::score, opt::multiThresh, opt::allowMisses);

		} else {
			ntHashIterator itr(seq, m_filter.getHashNum(),
					m_filter.getKmerSize());
			return m_filter.query(itr, m_perFrameProb, m_perFrameProbMulti,
					opt::score, opt::multiThresh, opt::allowMisses);
		}
	}

	/*
	 * testing heuristic faster code
	 */
	inline vector<pair<ID, double>> classify2(const string &seq) {
		unsigned frameCount = seq.size() - m_filter.getKmerSize() + 1;
		if (m_minCount.find(frameCount) == m_minCount.end()) {
			m_minCount[frameCount] = boost::shared_ptr<vector<unsigned>>(
					new vector<unsigned>(m_fullIDs.size()));
			for (size_t i = 0; i < m_minCount.size(); ++i) {
				//TODO remove hard coded values
				(*m_minCount[frameCount])[i] = getMinCount(frameCount,
						m_perFrameProb[i]);
			}
//			m_fprCount[frameCount] = boost::shared_ptr<vector<unsigned>>(
//					new vector<unsigned>(m_filter.getHashNum()));
//			for (unsigned i = 0; i < m_filter.getHashNum(); ++i) {
//				(*m_fprCount[frameCount])[i] = getFPRCount(frameCount,
//						pow(double(m_filter.getPop()) / double(m_filter.size()),
//								i + 1));
////				cerr << i << " " << m_fprCount[seq.size()]->at(i) << endl;
//			}
		}
		if (m_filter.getSeedValues().size() > 0) {
			RollingHashIterator itr(seq, m_filter.getKmerSize(),
					m_filter.getSeedValues());
			return m_filter.query(itr, *m_minCount[frameCount], frameCount);

		} else {
			ntHashIterator itr(seq, m_filter.getHashNum(),
					m_filter.getKmerSize());
			return m_filter.query(itr, *m_minCount[frameCount], frameCount);
		}
	}

	/*
	 * heuristic code
	 *
	 */
	inline vector<pair<ID, double>> classifyPair2(const string &seq1,
			const string &seq2) {
		unsigned frameCount = seq1.size() + seq2.size() - (m_filter.getKmerSize() + 1)*2;
		if (m_minCount.find(frameCount) == m_minCount.end()) {
			m_minCount[frameCount] = boost::shared_ptr<vector<unsigned>>(
					new vector<unsigned>(m_fullIDs.size()));
			for (size_t i = 0; i < m_minCount.size(); ++i) {
				//TODO remove hard coded values
				(*m_minCount[frameCount])[i] = getMinCount(frameCount,
						m_perFrameProb[i]);
			}
		}
		if (m_filter.getSeedValues().size() > 0) {
			RollingHashIterator itr1(seq1, m_filter.getKmerSize(),
					m_filter.getSeedValues());
			RollingHashIterator itr2(seq2, m_filter.getKmerSize(),
					m_filter.getSeedValues());
			return m_filter.query(itr1, itr2, *m_minCount[frameCount]);

		} else {
			ntHashIterator itr1(seq1, m_filter.getHashNum(),
					m_filter.getKmerSize());
			ntHashIterator itr2(seq2, m_filter.getHashNum(),
					m_filter.getKmerSize());
			return m_filter.query(itr1, itr2, *m_minCount[frameCount]);
		}
	}

	inline vector<pair<ID, double>> classifyPair(const string &seq1,
			const string &seq2) {
		if (m_filter.getSeedValues().size() > 0) {
			RollingHashIterator itr1(seq1, m_filter.getKmerSize(),
					m_filter.getSeedValues());
			RollingHashIterator itr2(seq2, m_filter.getKmerSize(),
					m_filter.getSeedValues());
			return m_filter.query(itr1, itr2, m_perFrameProb, opt::score,
					opt::allowMisses);

		} else {
			ntHashIterator itr1(seq1, m_filter.getHashNum(),
					m_filter.getKmerSize());
			ntHashIterator itr2(seq2, m_filter.getHashNum(),
					m_filter.getKmerSize());
			return m_filter.query(itr1, itr2, m_perFrameProb, opt::score,
					opt::allowMisses);
		}
	}

	inline void printPairToFile(const FaRec rec1, const FaRec rec2,
			Dynamicofstream &outputFile) {
		if (opt::outputType == "fa") {
#pragma omp critical(outputFiles)
			{
				outputFile << ">" << rec1.header << "\n" << rec1.seq << "\n";
				outputFile << ">" << rec2.header << "\n" << rec2.seq << "\n";
			}
		} else {
#pragma omp critical(outputFiles)
			{
				outputFile << "@" << rec1.header << "\n" << rec1.seq << "\n+\n"
						<< rec1.qual << "\n";
				outputFile << "@" << rec2.header << "\n" << rec2.seq << "\n+\n"
						<< rec2.qual << "\n";
			}
		}
	}

	inline void printPairToFile(unsigned outputFileIndex, const FaRec rec1,
			const FaRec rec2, vector<Dynamicofstream*> &outputFiles1,
			vector<Dynamicofstream*> &outputFiles2) {
		if (opt::outputType == "fa") {
#pragma omp critical(outputFiles)
			{
				(*outputFiles1[outputFileIndex]) << ">" << rec1.header << "\n"
						<< rec1.seq << "\n";
				(*outputFiles2[outputFileIndex]) << ">" << rec2.header << "\n"
						<< rec2.seq << "\n";
			}
		} else {
#pragma omp critical(outputFiles)
			{
				(*outputFiles1[outputFileIndex]) << "@" << rec1.header << "\n"
						<< rec1.seq << "\n+\n" << rec1.qual << "\n";
				(*outputFiles2[outputFileIndex]) << "@" << rec2.header << "\n"
						<< rec2.seq << "\n+\n" << rec2.qual << "\n";
			}
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

	/*
	 * Old function for evaluating pValues
	 * For debugging now
	 */
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
			if (criticalScore > cumProb) {
				break;
			}
		}
		return (i);
	}

	inline unsigned getFPRCount(unsigned length, double eventProb) {
		binomial bin(length, 1.0 - eventProb);
		unsigned i = 0;
		for (; i < length; ++i) {
			double cumProb = cdf(bin, length - i);
			if (opt::score > cumProb) {
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
