/*
 * BioBloomClassifier.h
 *
 *  Created on: Oct 17, 2012
 *      Author: cjustin
 */

#ifndef BIOBLOOMCLASSIFIER_H_
#define BIOBLOOMCLASSIFIER_H_
#include <vector>
#include <string>
#include "boost/unordered/unordered_map.hpp"
#include "boost/shared_ptr.hpp"
#include "Common/BloomFilterInfo.h"
#include "MultiFilter.h"
#include "Common/ReadsProcessor.h"
//#include "Common/Uncompress.h"
#include "Common/BloomFilter.h"
#include "ResultsManager.h"
#include "Common/Dynamicofstream.h"
#include "Common/SeqEval.h"
#include <zlib.h>
#include <iostream>

using namespace std;
using namespace boost;

static const string NO_MATCH = "noMatch";
static const string MULTI_MATCH = "multiMatch";

/** for modes of filtering */
enum mode {
	COLLAB, MINHITONLY, BESTHIT, STD, SCORES
};

///** for modes of printing out files */
//enum printMode {FASTA, FASTQ, BEST_FASTA, BEST_FASTQ};
//enum printMode {NORMAL, WITH_SCORE};

class BioBloomClassifier {
public:
	explicit BioBloomClassifier(const vector<string> &filterFilePaths,
			double scoreThreshold, const string &outputPrefix,
			const string &outputPostFix, unsigned minHit, bool minHitOnly,
			bool withScore);
	void filter(const vector<string> &inputFiles);
	void filterPrint(const vector<string> &inputFiles,
			const string &outputType);
	void filterPair(const string &file1, const string &file2);
	void filterPairPrint(const string &file1, const string &file2,
			const string &outputType);
	void filterPair(const string &file);
	void filterPairPrint(const string &file, const string &outputType);

	void setCollabFilter() {
		m_mode = COLLAB;
		if (m_hashSigs.size() != 1) {
			cerr
					<< "To use collaborative filtering all filters must use the same k and same number of hash functions."
					<< endl;
			exit(1);
		}
	}

	void setInclusive() {
		m_inclusive = true;
	}

	void setEvalMode(SeqEval::EvalMode mode) {
		m_evalMode = mode;
	}

	SeqEval::EvalMode getEvalMode() {
		return m_evalMode;
	}

	void setMainFilter(const string &filtername);

	virtual ~BioBloomClassifier();

private:
	//group filters with same hash number
	unordered_map<string, vector<boost::shared_ptr<BloomFilterInfo> > > m_infoFiles;
	unordered_map<string, boost::shared_ptr<MultiFilter> > m_filters;
	unordered_map<string, boost::shared_ptr<BloomFilter> > m_filtersSingle;
	vector<string> m_filterOrder;
	vector<string> m_hashSigs;
	double m_scoreThreshold;
	unsigned m_filterNum;
	const string &m_prefix;
	const string &m_postfix;
	const unsigned m_minHit;

	// modes of filtering
	mode m_mode;
	// Match scoring method. Possible values:
	// i) EVAL_STANDARD => score in range (0,1)
	// ii) EVAL_MIN_MATCH_LEN => minimum match length (in bases)
	SeqEval::EvalMode m_evalMode;

	string m_mainFilter;
	bool m_inclusive;

	struct FaRec {
		string header;
		string seq;
		string qual;
	};

	void loadFilters(const vector<string> &filterFilePaths);
	bool fexists(const string &filename) const;
	void evaluateReadStd(const string &rec, const string &hashSig,
			unordered_map<string, bool> &hits);
	void evaluateReadMin(const string &rec, const string &hashSig,
			unordered_map<string, bool> &hits);
	void evaluateReadCollab(const string &rec, const string &hashSig,
			unordered_map<string, bool> &hits);
	double evaluateReadBestHit(const string &rec, const string &hashSig,
			unordered_map<string, bool> &hits, vector<double> &scores);
	void evaluateReadScore(const string &rec, const string &hashSig,
			unordered_map<string, bool> &hits, vector<double> &scores);

	inline void printSingle(const FaRec &rec, double score,
			const string &filterID) {
		if (m_mainFilter == filterID) {
			if (m_mode == BESTHIT) {
#pragma omp critical(cout)
				{
					cout << "@" << rec.header << " " << score << "\n" << rec.seq
							<< "\n+\n" << rec.qual << "\n";
				}
			} else {
#pragma omp critical(cout)
				{
					cout << "@" << rec.header << "\n" << rec.seq << "\n+\n"
							<< rec.qual << "\n";
				}
			}
		}
	}

	inline void printSingleToFile(const string &outputFileName,
			const FaRec &rec,
			unordered_map<string, boost::shared_ptr<Dynamicofstream> > &outputFiles,
			string const &outputType, double score, vector<double> &scores) {
		if (outputType == "fa") {
			if (m_mode == SCORES && outputFileName == MULTI_MATCH) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << ">" << rec.header;
					for (vector<double>::iterator i = scores.begin();
							i != scores.end(); ++i) {
						(*outputFiles[outputFileName]) << " " << *i;
					}
					(*outputFiles[outputFileName]) << "\n" << rec.seq << "\n";
				}
			} else if (m_mode == BESTHIT) {
				if (outputFileName == MULTI_MATCH)
#pragma omp critical(outputFiles)
						{
					(*outputFiles[outputFileName]) << ">" << rec.header;
					for (vector<double>::iterator i = scores.begin();
							i != scores.end(); ++i) {
						(*outputFiles[outputFileName]) << " " << *i;
					}
					(*outputFiles[outputFileName]) << "\n" << rec.seq << "\n";
				} else
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << ">" << rec.header << " "
							<< score << "\n" << rec.seq << "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << ">" << rec.header << "\n"
							<< rec.seq << "\n";
				}
			}
		} else {
			if (m_mode == SCORES && outputFileName == MULTI_MATCH) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << "@" << rec.header;
					for (vector<double>::iterator i = scores.begin();
							i != scores.end(); ++i) {
						(*outputFiles[outputFileName]) << " " << *i;
					}
					(*outputFiles[outputFileName]) << "\n" << rec.seq << "\n+\n"
							<< rec.qual << "\n";
				}
			} else if (m_mode == BESTHIT) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << "@" << rec.header << " "
							<< score << "\n" << rec.seq << "\n+\n" << rec.qual
							<< "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << "@" << rec.header << "\n"
							<< rec.seq << "\n+\n" << rec.qual << "\n";
				}
			}
		}
	}

	inline void printPair(const FaRec &rec1, const FaRec &rec2, double score1,
			double score2, const string &filterID) {
		if (m_mainFilter == filterID) {
			if (m_mode == BESTHIT) {
#pragma omp critical(cout)
				{
					cout << "@" << rec1.header << " " << score1 << "\n"
							<< rec1.seq << "\n+\n" << rec1.qual << "\n";
					cout << "@" << rec2.seq << " " << score2 << "\n" << rec2.seq
							<< "\n+\n" << rec2.seq << "\n";
				}
			} else {
#pragma omp critical(cout)
				{
					cout << "@" << rec1.header << "\n" << rec1.seq << "\n+\n"
							<< rec1.qual << "\n";
					cout << "@" << rec2.seq << "\n" << rec2.seq << "\n+\n"
							<< rec2.seq << "\n";
				}
			}
		}
	}

	inline void printPairToFile(const string &outputFileName, const FaRec &rec1,
			const FaRec &rec2,
			unordered_map<string, boost::shared_ptr<Dynamicofstream> > &outputFiles,
			string const &outputType, double score1, double score2,
			vector<double> &scores1, vector<double> &scores2) {
		if (outputType == "fa") {
			if (m_mode == SCORES && outputFileName == MULTI_MATCH) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName + "_1"]) << ">" << rec1.header;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles[outputFileName + "_1"]) << " " << *i;
					}
					(*outputFiles[outputFileName + "_1"]) << "\n" << rec1.seq
							<< "\n";
					(*outputFiles[outputFileName + "_2"]) << ">" << rec2.seq;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles[outputFileName + "_2"]) << " " << *i;
					}
					(*outputFiles[outputFileName + "_2"]) << "\n" << rec2.seq
							<< "\n";
				}
			} else if (m_mode == BESTHIT) {
				if (outputFileName == MULTI_MATCH)
#pragma omp critical(outputFiles)
						{
					(*outputFiles[outputFileName + "_1"]) << ">" << rec1.header;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles[outputFileName + "_1"]) << " " << *i;
					}
					(*outputFiles[outputFileName + "_1"]) << "\n" << rec1.seq
							<< "\n";
					(*outputFiles[outputFileName + "_2"]) << ">" << rec2.seq;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles[outputFileName + "_2"]) << " " << *i;
					}
					(*outputFiles[outputFileName + "_2"]) << "\n" << rec2.seq
							<< "\n";
				} else
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName + "_1"]) << ">" << rec1.header
							<< " " << score1 << "\n" << rec1.seq << "\n";
					(*outputFiles[outputFileName + "_2"]) << ">" << rec2.seq
							<< " " << score2 << "\n" << rec2.seq << "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName + "_1"]) << ">" << rec1.header
							<< "\n" << rec1.seq << "\n";
					(*outputFiles[outputFileName + "_2"]) << ">" << rec2.seq
							<< "\n" << rec2.seq << "\n";
				}
			}
		} else {
			if (m_mode == SCORES && outputFileName == MULTI_MATCH) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName + "_1"]) << "@" << rec1.header;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles[outputFileName + "_1"]) << " " << *i;
					}
					(*outputFiles[outputFileName + "_1"]) << "\n" << rec1.seq
							<< "\n+\n" << rec1.qual << "\n";
					(*outputFiles[outputFileName + "_2"]) << "@" << rec2.seq;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles[outputFileName + "_2"]) << " " << *i;
					}
					(*outputFiles[outputFileName + "_2"]) << "\n" << rec2.seq
							<< "\n+\n" << rec2.seq << "\n";
				}
			} else if (m_mode == BESTHIT) {
				if (outputFileName == MULTI_MATCH)
#pragma omp critical(outputFiles)
						{
					(*outputFiles[outputFileName + "_1"]) << "@" << rec1.header;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles[outputFileName + "_1"]) << " " << *i;
					}
					(*outputFiles[outputFileName + "_1"]) << "\n" << rec1.seq
							<< "\n+\n" << rec1.qual << "\n";
					(*outputFiles[outputFileName + "_2"]) << "@" << rec2.seq;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles[outputFileName + "_2"]) << " " << *i;
					}
					(*outputFiles[outputFileName + "_2"]) << "\n" << rec2.seq
							<< "\n+\n" << rec2.seq << "\n";
				} else
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName + "_1"]) << "@" << rec1.header
							<< " " << score1 << "\n" << rec1.seq << "\n+\n"
							<< rec1.qual << "\n";
					(*outputFiles[outputFileName + "_2"]) << "@" << rec2.seq
							<< " " << score2 << "\n" << rec2.seq << "\n+\n"
							<< rec2.seq << "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName + "_1"]) << "@" << rec1.header
							<< "\n" << rec1.seq << "\n+\n" << rec1.qual << "\n";
					(*outputFiles[outputFileName + "_2"]) << "@" << rec2.seq
							<< "\n" << rec2.seq << "\n+\n" << rec2.seq << "\n";
				}
			}
		}
	}

	inline void evaluateRead(const string &rec, const string &hashSig,
			unordered_map<string, bool> &hits, double &score,
			vector<double> &scores) {
		switch (m_mode) {
		case COLLAB: {
			evaluateReadCollab(rec, hashSig, hits);
			break;
		}
		case MINHITONLY: {
			evaluateReadMin(rec, hashSig, hits);
			break;
		}
		case BESTHIT: {
			score = evaluateReadBestHit(rec, hashSig, hits, scores);
			break;
		}
		case SCORES: {
			evaluateReadScore(rec, hashSig, hits, scores);
			break;
		}
		default: {
			evaluateReadStd(rec, hashSig, hits);
			break;
		}
		}
	}
};

#endif /* BIOBLOOMCLASSIFIER_H_ */
