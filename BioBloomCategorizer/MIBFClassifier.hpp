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
#include "btl_bloomfilter/stHashIterator.hpp"
#include "btl_bloomfilter/ntHashIterator.hpp"
#include "bloomfilter/MIBFQuerySupport.hpp"
#include <iostream>
#include "Common/Options.h"
#include "ResultsManager.hpp"
//#include <BioBloomClassifier.h>
#include <tuple>
#include "Common/concurrentqueue.h"

#include <zlib.h>
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "Common/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/
#include "Common/kseq_util.h"

#include <boost/math/distributions/binomial.hpp>

using namespace std;
using namespace boost::math::policies;
using namespace boost::math;

static const std::string base64_chars =
		"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

const static size_t s_bulkSize = 256;

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
		} else if (opt::verbose) {
			cerr << "Loading ID file: " << idFile.c_str() << endl;
		}
		//first element is considered repetitive
		m_fullIDs.push_back("repeat");
		ifstream idFH(idFile.c_str(), ifstream::in);
		string line;
		getline(idFH, line);
		while (idFH.good()) {
			ID id;
			string name;
			stringstream converter(line);
			converter >> id;
			converter >> name;
			m_fullIDs.push_back(name);
			getline(idFH, line);
		}
		idFH.close();
		m_allowedMiss =
				(opt::frameMatches > m_filter.getHashNum())
						|| (m_filter.getSeedValues().size() == 0) ?
						0 : m_filter.getHashNum() - opt::frameMatches;
		if (opt::verbose) {
			cerr << "Calculating frame probabilities" << endl;
		}

		m_perFrameProb = vector<double>(m_fullIDs.size());
		m_rateSaturated = m_filter.calcFrameProbs(m_perFrameProb,
				m_allowedMiss);
		m_minCount.set_empty_key(0);
		if (opt::verbose) {
			cerr << "Filter loading complete" << endl;
		}
	}

	void filterDebug(const vector<string> &inputFiles) {
		cerr << "Filtering Start" << endl;

		//debugging
		size_t multiCount = 0;

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
			MIBFQuerySupport<ID> support = MIBFQuerySupport<ID>(m_filter,
					m_perFrameProb, opt::multiThresh, opt::streakThreshold,
					m_allowedMiss, opt::minCountNonSatCount,
					opt::bestHitCountAgree);
#pragma omp parallel private(l) firstprivate(support)
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
					vector<vector<pair<ID, bool>> > hitsPattern;
					unsigned evaluatedSeeds = 0;
					vector<unsigned> sig = support.getMatchSignature(sequence,
							evaluatedSeeds, hitsPattern);
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
					classify(support, sequence);

					const string &fullID = m_fullIDs.at(idIndex);
					unsigned tempCount = 0;
					for (unsigned i = 0; i < signifResults.size(); ++i) {
						//TODO resolve precision error better
						if (fullSignifCounts[i] >= maxCount) {
							++tempCount;
						}
					}
					if (tempCount > 1) {
#pragma omp atomic
						++multiCount;
					}

#pragma omp critical(cout)
					if (tempCount > 0) {
						cout << m_numRead << "\tCorrectName:" << name
								<< "\tPredictedName:" << fullID << "\tpVal:"
								<< log10(pVal) * (-10.0) << "\t" << endl;

						for (unsigned i = 0; i < signifResults.size(); ++i) {
							cout << signifResults[i] << ","
									<< m_fullIDs[signifResults[i]] << ","
									<< base64_chars[signifResults[i] % 64]
									<< "," << signifValues[i] << ","
									<< fullSignifCounts[i] << "," << pVals[i]
									<< "\t";
						}
						cout << endl;
						printVerbose(name, comment, sequence, hitsPattern, sig,
								idIndex);
//						printCountHistogram(hitsPattern);
						support.printAllCounts(m_fullIDs);
					}
				} else {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}

		cerr << "Multiple Map Count: " << multiCount << endl;
		cerr << "Total Reads: " << m_numRead << endl;
		cout.flush();
	}

	void filter(const vector<string> &inputFiles) {
		string outputName = opt::outputPrefix + "_reads.tsv";

		if (opt::outputType == opt::FASTQ) {
			outputName = opt::outputPrefix + "_reads.fq";
		} else if (opt::outputType == opt::FASTA) {
			outputName = opt::outputPrefix + "_reads.fa";
		}

		//print out header info and initialize variables
		ResultsManager<ID> resSummary(m_fullIDs, false);

		cerr << "Filtering Start" << endl;

		double startTime = omp_get_wtime();
		for (vector<string>::const_iterator it = inputFiles.begin();
				it != inputFiles.end(); ++it) {
			if (opt::threads == 1) {
				gzFile fp;
				fp = gzopen(it->c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				string outBuffer;
				MIBFQuerySupport<ID> support = MIBFQuerySupport<ID>(m_filter,
						m_perFrameProb, opt::multiThresh, opt::streakThreshold,
						m_allowedMiss, opt::minCountNonSatCount,
						opt::bestHitCountAgree);
				while (kseq_read(seq) >= 0) {
					filterSingleRead(*seq, support, resSummary,
							outBuffer);
				}
			} else {
				moodycamel::ConcurrentQueue<kseq_t> workQueue(
						opt::threads * s_bulkSize);
				bool good = true;
				typedef std::vector<kseq_t>::iterator iter_t;
#pragma omp parallel
				{
					std::vector<kseq_t> readBuffer;
					readBuffer.reserve(s_bulkSize);
					string outBuffer;
					MIBFQuerySupport<ID> support = MIBFQuerySupport<ID>(m_filter,
							m_perFrameProb, opt::multiThresh, opt::streakThreshold,
							m_allowedMiss, opt::minCountNonSatCount,
							opt::bestHitCountAgree);
					if (omp_get_thread_num() == 0) {
						//file reading init
						gzFile fp;
						fp = gzopen(it->c_str(), "r");
						kseq_t *seq = kseq_init(fp);

						//per thread token
						moodycamel::ProducerToken ptok(workQueue);

						unsigned size = 0;
						while (kseq_read(seq) >= 0) {
							readBuffer.push_back(kseq_t()); // Don't like this, need to reallocate memory twice
							cpy_kseq(&readBuffer[size++], seq); //TODO Make proper copy constructor for kseq?
							if (++m_numRead % opt::fileInterval == 0) {
								cerr << "Currently Reading Read Number: "
										<< m_numRead << endl;
							}
							if (s_bulkSize == size) {
								//try to insert, if cannot queue is full
								while (!workQueue.try_enqueue_bulk(ptok,
										std::move_iterator < iter_t
												> (readBuffer.begin()), size)) {
									//try to work
									if (kseq_read(seq) >= 0) {
										if (++m_numRead % opt::fileInterval
												== 0) {
											cerr
													<< "Currently Reading Read Number: "
													<< m_numRead << endl;
										}
										//------------------------WORK CODE START---------------------------------------
										filterSingleRead(*seq, support,
												resSummary, outBuffer);
										//------------------------WORK CODE END-----------------------------------------
									} else {
										break;
									}
								}
								//reset buffer
								readBuffer.clear();
								size = 0;
							}
						}
						//finish off remaining work
						for (unsigned i = 0; i < size; ++i) {
							//------------------------WORK CODE START---------------------------------------
							filterSingleRead(readBuffer[i], support, resSummary,
									outBuffer);
							//------------------------WORK CODE END-----------------------------------------
						}
						readBuffer.clear();
						if (workQueue.size_approx()) {
							moodycamel::ConsumerToken ctok(workQueue);
							//join in if others are still not finished
							while (workQueue.size_approx()) {
								size_t num = workQueue.try_dequeue_bulk(ctok, readBuffer.begin(),
										s_bulkSize);
								if (num) {
									for (unsigned i = 0; i < num; ++i) {
										//------------------------WORK CODE START---------------------------------------
										filterSingleRead(readBuffer[i], support,
												resSummary, outBuffer);
										//------------------------WORK CODE END-----------------------------------------
									}
								}
							}
						}
						good = false;
						kseq_destroy(seq);
						gzclose(fp);
					}
					else {
						moodycamel::ConsumerToken ctok(workQueue);
						while (good) {
							if (workQueue.size_approx() >= s_bulkSize) {
								size_t num = workQueue.try_dequeue_bulk(ctok, readBuffer.begin(),
										s_bulkSize);
								if (num) {
									for (unsigned i = 0; i < num; ++i) {
										//------------------------WORK CODE START---------------------------------------
										filterSingleRead(readBuffer[i], support,
												resSummary, outBuffer);
										//------------------------WORK CODE END-----------------------------------------
									}
								}
							}
						}
					}
				}
			}
		}

		cerr << "Classification time (s): " << (omp_get_wtime() - startTime)
				<< endl;

		cerr << "Total Reads:" << m_numRead << endl;
		cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv"
				<< endl;

		ofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
		summaryOutput << resSummary.getResultsSummary(m_numRead);
		summaryOutput.close();
		cout.flush();
	}

	void filterPair(const string &file1, const string &file2) {

		//results summary object
		ResultsManager<ID> resSummary(m_fullIDs, false);
		if (opt::verbose) {
			cerr << "Filtering Start" << "\n";
		}

		moodycamel::ConcurrentQueue<pair<kseq_t, kseq_t>> workQueue(
				opt::threads * s_bulkSize);

		typedef std::vector<pair<kseq_t, kseq_t>>::iterator iter_t;

		double startTime = omp_get_wtime();
		bool good = true;

#pragma omp parallel
		{
			vector<pair<kseq_t, kseq_t>> readBuffer;
			readBuffer.reserve(s_bulkSize);
			MIBFQuerySupport<ID> support = MIBFQuerySupport<ID>(m_filter,
					m_perFrameProb, opt::multiThresh, opt::streakThreshold,
					m_allowedMiss, opt::minCountNonSatCount,
					opt::bestHitCountAgree);
			string outBuffer;
			if (omp_get_thread_num() == 0) {
				//file reading init
				gzFile fp1 = gzopen(file1.c_str(), "r");
				gzFile fp2 = gzopen(file2.c_str(), "r");
				kseq_t *seq1 = kseq_init(fp1);
				kseq_t *seq2 = kseq_init(fp2);

				//per thread token
				moodycamel::ProducerToken ptok(workQueue);

				unsigned size = 0;
				while ((kseq_read(seq1) >= 0) && (kseq_read(seq2) >= 0)) {
					readBuffer.push_back(pair<kseq_t,kseq_t>(kseq_t(),kseq_t()));
					cpy_kseq(&(readBuffer[size].first), seq1);
					cpy_kseq(&(readBuffer[size++].second), seq2);
					if (++m_numRead % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: "
								<< m_numRead << endl;
					}
					if (s_bulkSize == size) {
						//try to insert, if cannot queue is full
						while (!workQueue.try_enqueue_bulk(ptok,
								std::move_iterator < iter_t > (readBuffer.begin()),
								size)) {
							if (++m_numRead % opt::fileInterval == 0) {
								cerr << "Currently Reading Read Number: "
										<< m_numRead << endl;
							}
							//try to work
							if ((kseq_read(seq1) >= 0) && (kseq_read(seq2) >= 0)) {
								filterPairedRead(*seq1, *seq2, support,
										resSummary, outBuffer);
							} else {
								break;
							}
						}
						readBuffer.clear();
						size = 0;
					}
				}
				//finish off remaining work
				for (unsigned i = 0; i < size; ++i) {
					filterPairedRead(readBuffer[i].first,
							readBuffer[i].second, support,
							resSummary, outBuffer);
				}
				if (workQueue.size_approx()) {
					moodycamel::ConsumerToken ctok(workQueue);
					//join in if others are still not finished
					while (workQueue.size_approx()) {
						size_t num = workQueue.try_dequeue_bulk(ctok,
								std::move_iterator < iter_t > (readBuffer.begin()),
								s_bulkSize);
						if (num) {
							for (unsigned i = 0; i < num; ++i) {
								filterPairedRead(readBuffer[i].first,
										readBuffer[i].second, support,
										resSummary, outBuffer);
							}
						}
					}
				}
#pragma omp atomic update
				good &= false;
				kseq_destroy(seq1);
				kseq_destroy(seq2);
				gzclose(fp1);
				gzclose(fp2);
			} else {
				moodycamel::ConsumerToken ctok(workQueue);
				while (good) {
					if (workQueue.size_approx() >= s_bulkSize) {
						size_t num = workQueue.try_dequeue_bulk(ctok,
								std::move_iterator < iter_t > (readBuffer.begin()),
								s_bulkSize);
						if (num) {
							for (unsigned i = 0; i < num; ++i) {
								filterPairedRead(readBuffer[i].first,
										readBuffer[i].second, support,
										resSummary, outBuffer);
							}
						}
					}
				}
			}
		}
		cerr << "Classification time (s): " << (omp_get_wtime() - startTime)
				<< endl;

		ofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
		summaryOutput << resSummary.getResultsSummary(m_numRead);
		summaryOutput.close();
		cout.flush();
	}

private:
	MIBloomFilter<ID> m_filter;
	size_t m_numRead;
	vector<string> m_fullIDs;
	vector<double> m_perFrameProb;
	google::dense_hash_map<unsigned, boost::shared_ptr<vector<unsigned>>>m_minCount;
	double m_rateSaturated;
	unsigned m_allowedMiss;
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

	//Relies on short string optimization
	void appendResults(const vector<MIBFQuerySupport<ID>::QueryResult> &signifResults,
			string &outStr) {
		outStr += "\t";
		for (unsigned i = 0; i < signifResults.size(); ++i) {
			if (i != 0) {
				outStr += ";";
			}
			outStr += m_fullIDs[signifResults[i].id];
			outStr += ",";
			outStr += std::to_string((unsigned) signifResults[i].nonSatFrameCount);
			outStr += ",";
			outStr += std::to_string((unsigned) signifResults[i].solidCount);
			outStr += ",";
			outStr += std::to_string((unsigned) signifResults[i].count);
			outStr += ",";
			outStr += std::to_string((unsigned) signifResults[i].nonSatCount);
			outStr += ",";
			outStr += std::to_string((unsigned) signifResults[i].totalNonSatCount);
			outStr += ",";
			outStr += std::to_string((unsigned) signifResults[i].totalCount);
		}
	}

	void formatOutStr(const kseq_t &read, string &outStr, MIBFQuerySupport<ID> &support,
			const vector<MIBFQuerySupport<ID>::QueryResult> &signifResults) {
		if (opt::outputType == opt::FASTA) {
			outStr += ">";
			outStr += read.name.s;
			if(read.comment.l) {
				outStr += " ";
				outStr += read.comment.s;
			}
			if (!signifResults.empty()) {
				unsigned i = 0;
				outStr += "\t";
				outStr += m_fullIDs[signifResults[i].id];
				if(signifResults.size() > 1) {
					appendResults(signifResults, outStr);
				}
			}
			outStr += "\n";
			outStr += read.seq.s;
		} else if (opt::outputType == opt::FASTQ) {
			outStr += "@";
			outStr += read.name.s;
			if(read.comment.l) {
				outStr += " ";
				outStr += read.comment.s;
			}
			if (!signifResults.empty()) {
				unsigned i = 0;
				outStr += "\t";
				outStr += m_fullIDs[signifResults[i].id];
				if(signifResults.size() > 1) {
					appendResults(signifResults, outStr);
				}
			}
			outStr += "\n";
			outStr += read.seq.s;
			outStr += "\n+";
			outStr += read.qual.s;
		}
		else {
			if (signifResults.empty()) {
				outStr += "*\t";
				outStr += read.name.s;
				if(read.comment.l) {
					outStr += " ";
					outStr += read.comment.s;
				}
				outStr += "\t";
				outStr += std::to_string(support.getSatCount());
				outStr += "\t*";
			} else {
				unsigned i = 0;
				outStr += m_fullIDs[signifResults[i].id];
				outStr += "\t";
				outStr += read.name.s;
				if(read.comment.l) {
					outStr += " ";
					outStr += read.comment.s;
				}
				outStr += "\t";
				outStr += std::to_string(support.getSatCount());
				if (signifResults.size() == 1) {
					outStr += "\t*";
				} else {
					appendResults(signifResults, outStr);
				}
			}
		}
		outStr += "\n";
	}

	void filterSingleRead(const kseq_t &read, MIBFQuerySupport<ID> &support,
			ResultsManager<ID> &resSummary, string &outStr) {
		const vector<MIBFQuerySupport<ID>::QueryResult> &signifResults =
				classify(support, read.seq.s);
		resSummary.updateSummaryData(signifResults);
		outStr.clear();
		formatOutStr(read, outStr, support, signifResults);
		cout << outStr;
	}

	void filterPairedRead(const kseq_t &read1, const kseq_t &read2, MIBFQuerySupport<ID> &support,
			ResultsManager<ID> &resSummary, string &outStr) {
		const vector<MIBFQuerySupport<ID>::QueryResult> &signifResults =
		classify(support, read1.seq.s, read2.seq.s);
		outStr.clear();
		resSummary.updateSummaryData(signifResults);
		formatOutStr(read1, outStr, support, signifResults);
		formatOutStr(read2, outStr, support, signifResults);
		cout << outStr;
	}

	/*
	 * testing heuristic faster code
	 */
	//TODO: Reuse itr object
	inline const vector<MIBFQuerySupport<ID>::QueryResult> &classify(MIBFQuerySupport<ID> &support, const string &seq) {
		unsigned frameCount = seq.size() - m_filter.getKmerSize() + 1;
#pragma omp critical(minCount)
		if (m_minCount.find(frameCount) == m_minCount.end()) {
			m_minCount[frameCount] = boost::shared_ptr<vector<unsigned>> (
					new vector<unsigned>(m_fullIDs.size()));
			for (size_t i = 1; i < m_fullIDs.size(); ++i) {
				(*m_minCount[frameCount])[i] = getMinCount(frameCount,
						m_perFrameProb[i]);
			}
		}
		if (m_filter.getSeedValues().size() > 0) {
			stHashIterator itr(seq, m_filter.getSeedValues(),
					m_filter.getHashNum(), m_filter.getKmerSize());
			return support.query(itr, *m_minCount[frameCount]);
		} else {
			ntHashIterator itr(seq, m_filter.getHashNum(),
					m_filter.getKmerSize());
			return support.query(itr, *m_minCount[frameCount]);
		}
	}

	/*
	 * heuristic code
	 *
	 */
	inline const vector<MIBFQuerySupport<ID>::QueryResult> &classify(MIBFQuerySupport<ID> &support, const string &seq1,
			const string &seq2) {
		unsigned frameCount = seq1.size() + seq2.size() - (m_filter.getKmerSize() + 1)*2;
#pragma omp critical(m_minCount)
		if (m_minCount.find(frameCount) == m_minCount.end()) {
			m_minCount[frameCount] = boost::shared_ptr<vector<unsigned>>(
					new vector<unsigned>(m_fullIDs.size()));
			for (size_t i = 1; i < m_fullIDs.size(); ++i) {
				(*m_minCount[frameCount])[i] = getMinCount(frameCount,
						m_perFrameProb[i]);
			}
		}
		if (m_filter.getSeedValues().size() > 0) {
			stHashIterator itr1(seq1, m_filter.getSeedValues(), m_filter.getHashNum(),
					m_filter.getKmerSize());
			stHashIterator itr2(seq2, m_filter.getSeedValues(), m_filter.getHashNum(),
					m_filter.getKmerSize());
			return support.query(itr1, itr2, *m_minCount[frameCount]);
		} else {
			ntHashIterator itr1(seq1, m_filter.getHashNum(),
					m_filter.getKmerSize());
			ntHashIterator itr2(seq2, m_filter.getHashNum(),
					m_filter.getKmerSize());
			return support.query(itr1, itr2, *m_minCount[frameCount]);
		}
	}

	/*
	 * Old function for evaluating pValues
	 * For debugging now
	 */
	inline ID evalRead(const vector<vector<pair<ID, bool> > > &hitsPattern,
			unsigned evaluatedSeeds, double &pVal, unsigned &maxCount,
			vector<ID> &signifResults, vector<unsigned> &signifCounts,
			vector<unsigned> &fullSignifCounts, vector<double> &pVals) {
		google::dense_hash_map<ID, unsigned> counts;
		google::dense_hash_map<ID, unsigned> fullCounts;
		counts.set_empty_key(opt::EMPTY);
		fullCounts.set_empty_key(opt::EMPTY);
		for (unsigned i = 0;
				i < hitsPattern.size(); i++) {
			//to determine if already added for this frame
			google::dense_hash_set<ID> tempIDs;
			tempIDs.set_empty_key(opt::EMPTY);
			unsigned count = 0;
			for (vector<pair<ID, bool>>::const_iterator j = hitsPattern[i].begin(); j != hitsPattern[i].end();
					j++) {
				if (j->first != opt::EMPTY) {
					if (tempIDs.find(j->first) == tempIDs.end()) {
						google::dense_hash_map<ID, unsigned>::iterator tempItr =
						counts.find(j->first);
						if (tempItr == counts.end()) {
							counts[j->first] = 1;
						} else {
							++(tempItr->second);
						}
						tempIDs.insert(j->first);
					}
					++count;
				}
			}
			//make full counts consider only non-saturation
			if (count == m_filter.getHashNum()) {
				google::dense_hash_set<ID> tempIDsFull;
				tempIDsFull.set_empty_key(opt::EMPTY);
				for (vector<pair<ID, bool>>::const_iterator j = hitsPattern[i].begin(); j != hitsPattern[i].end();
						j++) {
					if (tempIDsFull.find(j->first) == tempIDsFull.end() && !j->second) {
						google::dense_hash_map<ID, unsigned>::iterator tempItrFull =
						fullCounts.find(j->first);
						if (tempItrFull == fullCounts.end()) {
							fullCounts[j->first] = 1;
						} else {
							++(tempItrFull->second);
						}
						tempIDsFull.insert(j->first);
					}
				}
			}
		}

		ID maxID = opt::EMPTY;
		double minVal = 1.0;
		double maxMinVal = 1.0;
		double adjustedPValThreshold = opt::score / double(m_fullIDs.size() - 1);
		pVal = 1.0;
		for (google::dense_hash_map<ID, unsigned>::const_iterator itr =
				counts.begin(); itr != counts.end(); itr++) {
			typedef boost::math::binomial_distribution<
			double,
			policy<discrete_quantile<integer_round_up> > >
			binom_round_up;
			binom_round_up bin(evaluatedSeeds, m_perFrameProb[itr->first]);
			double cumProb = cdf(complement(bin, itr->second));
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
		//bonferroni
		pVal = minVal * (m_fullIDs.size() - 1);
		pVal = pVal > 1 ? 1.0 : pVal;

		if (pVal <= opt::score) {
			return maxID;
		}
		return opt::EMPTY;
	}

	inline unsigned getMinCount(unsigned length, double eventProb) {
		typedef boost::math::binomial_distribution<
		double,
		policy<discrete_quantile<integer_round_up> > >
		binom_round_up;
		binom_round_up bin(length, eventProb);
		double criticalScore = opt::score / double(m_fullIDs.size() - 1);
		unsigned i = quantile(complement(bin, criticalScore));
		return i > 1 ? i : 1;
	}

	//debug helper methods
	/*
	 * debugging function
	 * prints this count hist for all elements
	 * allCounts, counts
	 * TODO satuCounts, allSatuCounts
	 */
	void printCountHistogram(const vector<vector<pair<ID,bool>> > &hitsPattern) {
		google::dense_hash_map< ID, unsigned> allCounts;
		allCounts.set_empty_key(0);
		google::dense_hash_map< ID, unsigned> counts;
		counts.set_empty_key(0);
		google::dense_hash_map< ID, unsigned> satCounts;
		satCounts.set_empty_key(0);
		for(ID i = 1; i < m_fullIDs.size(); ++i) {
			allCounts[i] = 0;
			counts[i] = 0;
			satCounts[i] = 0;
		}
		for (vector<vector<pair<ID,bool>> >::const_iterator i = hitsPattern.begin();
				i != hitsPattern.end(); i++) {
			//to determine if already added for this frame
			google::dense_hash_set<ID> tempIDs;
			tempIDs.set_empty_key(opt::EMPTY);
			for (vector<pair<ID,bool>>::const_iterator j = i->begin(); j != i->end();
					j++) {
				if (j->first != opt::EMPTY) {
					if (tempIDs.find(j->first) == tempIDs.end()) {
						tempIDs.insert(j->first);
						++counts[j->first];
					}
					if(!j->second) {
						++satCounts[j->first];
					}
					++allCounts[j->first];
				}
			}
		}

		for (ID i = 1; i < m_fullIDs.size(); ++i) {
			if(allCounts[i] > 0) {
				cout << i << "\t" << m_fullIDs[i] << "\t" << counts[i] << "\t"
				<< allCounts[i] << "\t"
				<< satCounts[i] << endl;
			}
		}
	}

	inline void printVerbose(const string &header, const string &comment,
			const string &seq, const vector<vector<pair<ID,bool>> > &hitsPattern,
			const vector<unsigned> &sig, ID value) {
		unsigned evaluatedSeeds = 0;
		cout << header << ' ' << comment << ' ' << evaluatedSeeds << ' '
		<< base64_chars[value % 64] << ' ' << value << endl;
		cout << seq << endl;
		cout << vectToStr(sig, seq) << endl;
		cout << vectToStr(sig, hitsPattern, seq);
	}

	inline string vectToStr(const vector<unsigned> &hitsVector,
			const vector<vector<pair<ID,bool>> > &hitsPattern, const string &seq) {
		stringstream ss;
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
					if (hitsPattern[index][hVal].first != opt::EMPTY) {
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
					if (hitsPattern[index][hVal].first == opt::EMPTY) {
						ss << " ";
					} else {
						ss << base64_chars[hitsPattern[index][hVal].first % 64];
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
					if (hitsPattern[index][hVal].second == opt::EMPTY) {
						ss << 0;
					} else {
						ss << 1;
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
