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
#include "google/dense_hash_map"
#include "Common/BloomFilterInfo.h"
#include "Common/Dynamicofstream.h"
#include "Common/SeqEval.h"
#include <zlib.h>
#include <iostream>
#include "ResultsManager.hpp"
#include "BioBloomCategorizer/Options.h"
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "Common/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

/** for modes of filtering */
enum mode {
	ORDERED, BESTHIT, STD, SCORES
};

struct FaRec {
	string header;
	string seq;
	string qual;
	string comment;
};

///** for modes of printing out files */
//enum printMode {FASTA, FASTQ, BEST_FASTA, BEST_FASTQ};
//enum printMode {NORMAL, WITH_SCORE};

//To prevent unused variable warning
static void __attribute__((unused)) wno_unused_kseq(void) {
	(void) &kseq_init;
	(void) &kseq_read;
	(void) &kseq_destroy;
	return;
}

class BioBloomClassifier {
public:
	explicit BioBloomClassifier(const vector<string> &filterFilePaths,
			double scoreThreshold, const string &outputPrefix,
			const string &outputPostFix, bool withScore);
	void filter(const vector<string> &inputFiles);
	void filterPrint(const vector<string> &inputFiles,
			const string &outputType);
	void filterPair(const string &file1, const string &file2);
	void filterPairPrint(const string &file1, const string &file2,
			const string &outputType);
	void filterPair(const string &file);
	void filterPairPrint(const string &file, const string &outputType);

	//file list functions
	void filterPair(const vector<string> &inputFiles1,
			const vector<string> &inputFiles2);
	void filterPairPrint(const vector<string> &inputFiles1,
			const vector<string> &inputFiles2, const string &outputType);

	void setOrderedFilter() {
		if (m_mode == BESTHIT) {
			cerr
					<< "Best Hit mode and Ordered mode detected. Not yet supported."
					<< endl;
			exit(1);
		}
		m_mode = ORDERED;
	}

	void setInclusive() {
		m_inclusive = true;
	}

	void setStdout() {
		m_stdout = true;
	}

	virtual ~BioBloomClassifier();

private:
	//TODO: change some of these variable to static global variable in option namespace
	vector<BloomFilterInfo*> m_infoFiles;
	vector<BloomFilter*> m_filters;
	vector<string> m_filterOrder;
	double m_scoreThreshold;
	const unsigned m_filterNum;
	const string &m_prefix;
	const string &m_postfix;

	// modes of filtering
	mode m_mode;

	bool m_stdout;
	bool m_inclusive;

	void loadFilters(const vector<string> &filterFilePaths);
	void evaluateReadStd(const string &rec, vector<unsigned> &hits);
//	void evaluateReadMin(const string &rec, vector<unsigned> &hits);
//	void evaluateReadCollab(const string &rec, vector<unsigned> &hits);
	void evaluateReadOrdered(const string &rec, vector<unsigned> &hits);
	double evaluateReadBestHit(const string &rec, vector<unsigned> &hits,
			vector<double> &scores);
	void evaluateReadScore(const string &rec, vector<unsigned> &hits,
			vector<double> &scores);

//	void evaluateReadCollabPair(const string &rec1, const string &rec2,
//			vector<unsigned> &hits1, vector<unsigned> &hits2);
	void evaluateReadOrderedPair(const string &rec1, const string &rec2,
			vector<unsigned> &hits1, vector<unsigned> &hits2);

	inline void printSingle(const FaRec &rec, double score, unsigned filterID) {
		if (m_stdout) {
			if (filterID == 0 && !opt::inverse) {
				if (m_mode == BESTHIT) {
#pragma omp critical(cout)
					{
						cout << "@" << rec.header << " " << rec.comment << " " << score << "\n"
								<< rec.seq << "\n+\n" << rec.qual << "\n";
					}
				} else {
#pragma omp critical(cout)
					{
						cout << "@" << rec.header << " " << rec.comment << "\n" << rec.seq << "\n+\n"
								<< rec.qual << "\n";
					}
				}
			} else if (filterID != 0) {
#pragma omp critical(cout)
				{
					cout << "@" << rec.header << " " << rec.comment << "\n" << rec.seq << "\n+\n"
							<< rec.qual << "\n";
				}
			}
		}
	}

	inline void printSingleToFile(unsigned outputFileName, const FaRec &rec,
			vector<Dynamicofstream*> &outputFiles, string const &outputType,
			double score, vector<double> &scores,
			const ResultsManager<unsigned> &rm) {
		if (outputType == "fa") {
			if (m_mode == SCORES && outputFileName == rm.getMultiMatchIndex()) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << ">" << rec.header << " " << rec.comment;
					for (vector<double>::iterator i = scores.begin();
							i != scores.end(); ++i) {
						(*outputFiles[outputFileName]) << " " << *i;
					}
					(*outputFiles[outputFileName]) << "\n" << rec.seq << "\n";
				}
			} else if (m_mode == BESTHIT) {
				if (outputFileName == rm.getMultiMatchIndex())
#pragma omp critical(outputFiles)
						{
					(*outputFiles[outputFileName]) << ">" << rec.header << " " << rec.comment;
					for (vector<double>::iterator i = scores.begin();
							i != scores.end(); ++i) {
						(*outputFiles[outputFileName]) << " " << *i;
					}
					(*outputFiles[outputFileName]) << "\n" << rec.seq << "\n";
				} else
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << ">" << rec.header << " " << rec.comment << " "
							<< score << "\n" << rec.seq << "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << ">" << rec.header << " " << rec.comment << "\n"
							<< rec.seq << "\n";
				}
			}
		} else {
			if (m_mode == SCORES && outputFileName == rm.getMultiMatchIndex()) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << "@" << rec.header << " " << rec.comment;
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
					(*outputFiles[outputFileName]) << "@" << rec.header << " " << rec.comment << " "
							<< score << "\n" << rec.seq << "\n+\n" << rec.qual
							<< "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << "@" << rec.header << " " << rec.comment << "\n"
							<< rec.seq << "\n+\n" << rec.qual << "\n";
				}
			}
		}
	}

	inline void printPair(const FaRec &rec1, const FaRec &rec2, double score1,
			double score2, unsigned filterID) {
		if (m_stdout) {
			if (filterID == 0 && !opt::inverse) {
				if (m_mode == BESTHIT) {
#pragma omp critical(cout)
					{
						cout << "@" << rec1.header << " " << rec1.comment << " " << score1 << "\n"
								<< rec1.seq << "\n+\n" << rec1.qual << "\n";
						cout << "@" << rec2.header << " " << rec2.comment << " " << score2 << "\n"
								<< rec2.seq << "\n+\n" << rec2.qual << "\n";
					}
				} else {
#pragma omp critical(cout)
					{
						cout << "@" << rec1.header << " " << rec1.comment << "\n" << rec1.seq
								<< "\n+\n" << rec1.qual << "\n";
						cout << "@" << rec2.header << " " << rec2.comment << "\n" << rec2.seq
								<< "\n+\n" << rec2.qual << "\n";
					}
				}
			} else if (filterID != 0) {
#pragma omp critical(cout)
				{
					cout << "@" << rec1.header << " " << rec1.comment << "\n" << rec1.seq
							<< "\n+\n" << rec1.qual << "\n";
					cout << "@" << rec2.header << " " << rec2.comment << "\n" << rec2.seq
							<< "\n+\n" << rec2.qual << "\n";
				}
			}
		}
	}

	inline void printPair(const kseq_t * rec1, const kseq_t * rec2,
			double score1, double score2, unsigned filterID) {
		if (m_stdout) {
			if (filterID == 0 && !opt::inverse) {
				if (m_mode == BESTHIT) {
#pragma omp critical(cout)
					{
						cout << "@" << rec1->name.s << " " << rec1->comment.s << " " << score1 << "\n"
								<< rec1->seq.s << "\n+\n" << rec1->qual.s
								<< "\n";
						cout << "@" << rec2->name.s << " " << rec2->comment.s << " " << score2 << "\n"
								<< rec2->seq.s << "\n+\n" << rec2->qual.s
								<< "\n";
					}
				} else {
#pragma omp critical(cout)
					{
						cout << "@" << rec1->name.s << " " << rec1->comment.s << "\n" << rec1->seq.s
								<< "\n+\n" << rec1->qual.s << "\n";
						cout << "@" << rec2->name.s << " " << rec2->comment.s << "\n" << rec2->seq.s
								<< "\n+\n" << rec2->qual.s << "\n";
					}
				}
			} else if (filterID != 0) {
#pragma omp critical(cout)
				{
					cout << "@" << rec1->name.s << " " << rec1->comment.s << "\n" << rec1->seq.s
							<< "\n+\n" << rec1->qual.s << "\n";
					cout << "@" << rec2->name.s << " " << rec2->comment.s << "\n" << rec2->seq.s
							<< "\n+\n" << rec2->qual.s << "\n";
				}
			}
		}
	}

	inline void printPairToFile(unsigned outputFileIndex, const FaRec &rec1,
			const FaRec &rec2, vector<Dynamicofstream*> &outputFiles1,
			vector<Dynamicofstream*> &outputFiles2, string const &outputType,
			double score1, double score2, vector<double> &scores1,
			vector<double> &scores2, const ResultsManager<unsigned> &rm) {
		if (outputType == "fa") {
			if (m_mode == SCORES
					&& outputFileIndex == rm.getMultiMatchIndex()) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << ">" << rec1.header << " " << rec1.comment;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles1[outputFileIndex]) << " " << *i;
					}
					(*outputFiles1[outputFileIndex]) << "\n" << rec1.seq
							<< "\n";
					(*outputFiles2[outputFileIndex]) << ">" << rec2.header << " " << rec2.comment;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles2[outputFileIndex]) << " " << *i;
					}
					(*outputFiles2[outputFileIndex]) << "\n" << rec2.seq
							<< "\n";
				}
			} else if (m_mode == BESTHIT) {
				if (outputFileIndex == rm.getMultiMatchIndex())
#pragma omp critical(outputFiles)
						{
					(*outputFiles1[outputFileIndex]) << ">" << rec1.header << " " << rec1.comment;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles1[outputFileIndex]) << " " << *i;
					}
					(*outputFiles1[outputFileIndex]) << "\n" << rec1.seq
							<< "\n";
					(*outputFiles2[outputFileIndex]) << ">" << rec2.header << " " << rec2.comment;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles2[outputFileIndex]) << " " << *i;
					}
					(*outputFiles2[outputFileIndex]) << "\n" << rec2.seq
							<< "\n";
				} else
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << ">" << rec1.header << " " << rec1.comment
							<< " " << score1 << "\n" << rec1.seq << "\n";
					(*outputFiles2[outputFileIndex]) << ">" << rec2.header << " " << rec2.comment
							<< " " << score2 << "\n" << rec2.seq << "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << ">" << rec1.header << " " << rec1.comment
							<< "\n" << rec1.seq << "\n";
					(*outputFiles2[outputFileIndex]) << ">" << rec2.header << " " << rec2.comment
							<< "\n" << rec2.seq << "\n";
				}
			}
		} else {
			if (m_mode == SCORES
					&& outputFileIndex == rm.getMultiMatchIndex()) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << "@" << rec1.header << " " << rec1.comment;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles1[outputFileIndex]) << " " << *i;
					}
					(*outputFiles1[outputFileIndex]) << "\n" << rec1.seq
							<< "\n+\n" << rec1.qual << "\n";
					(*outputFiles2[outputFileIndex]) << "@" << rec2.header << " " << rec2.comment;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles2[outputFileIndex]) << " " << *i;
					}
					(*outputFiles2[outputFileIndex]) << "\n" << rec2.seq
							<< "\n+\n" << rec2.qual << "\n";
				}
			} else if (m_mode == BESTHIT) {
				if (outputFileIndex == rm.getMultiMatchIndex())
#pragma omp critical(outputFiles)
						{
					(*outputFiles1[outputFileIndex]) << "@" << rec1.header << " " << rec1.comment;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles1[outputFileIndex]) << " " << *i;
					}
					(*outputFiles1[outputFileIndex]) << "\n" << rec1.seq
							<< "\n+\n" << rec1.qual << "\n";
					(*outputFiles2[outputFileIndex]) << "@" << rec2.header << " " << rec2.comment;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles2[outputFileIndex]) << " " << *i;
					}
					(*outputFiles2[outputFileIndex]) << "\n" << rec2.seq
							<< "\n+\n" << rec2.qual << "\n";
				} else
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << "@" << rec1.header << " " << rec1.comment
							<< " " << score1 << "\n" << rec1.seq << "\n+\n"
							<< rec1.qual << "\n";
					(*outputFiles2[outputFileIndex]) << "@" << rec2.header << " " << rec2.comment
							<< " " << score2 << "\n" << rec2.seq << "\n+\n"
							<< rec2.qual << "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << "@" << rec1.header << " " << rec1.comment
							<< "\n" << rec1.seq << "\n+\n" << rec1.qual << "\n";
					(*outputFiles2[outputFileIndex]) << "@" << rec2.header << " " << rec2.comment
							<< "\n" << rec2.seq << "\n+\n" << rec2.qual << "\n";
				}
			}
		}
	}

	inline void printPairToFile(unsigned outputFileIndex, const kseq_t * rec1,
			const kseq_t * rec2, vector<Dynamicofstream*> &outputFiles1,
			vector<Dynamicofstream*> &outputFiles2, string const &outputType,
			double score1, double score2, vector<double> &scores1,
			vector<double> &scores2, const ResultsManager<unsigned> &rm) {
		if (outputType == "fa") {
			if (m_mode == SCORES
					&& outputFileIndex == rm.getMultiMatchIndex()) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << ">" << rec1->name.s << " " << rec1->comment.s;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles1[outputFileIndex]) << " " << *i;
					}
					(*outputFiles1[outputFileIndex]) << "\n" << rec1->seq.s
							<< "\n";
					(*outputFiles2[outputFileIndex]) << ">" << rec2->name.s << " " << rec2->comment.s;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles2[outputFileIndex]) << " " << *i;
					}
					(*outputFiles2[outputFileIndex]) << "\n" << rec2->seq.s
							<< "\n";
				}
			} else if (m_mode == BESTHIT) {
				if (outputFileIndex == rm.getMultiMatchIndex())
#pragma omp critical(outputFiles)
						{
					(*outputFiles1[outputFileIndex]) << ">" << rec1->name.s << " " << rec1->comment.s;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles1[outputFileIndex]) << " " << *i;
					}
					(*outputFiles1[outputFileIndex]) << "\n" << rec1->seq.s
							<< "\n";
					(*outputFiles2[outputFileIndex]) << ">" << rec2->name.s << " " << rec2->comment.s;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles2[outputFileIndex]) << " " << *i;
					}
					(*outputFiles2[outputFileIndex]) << "\n" << rec2->seq.s
							<< "\n";
				} else
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << ">" << rec1->name.s << " " << rec1->comment.s
							<< " " << score1 << "\n" << rec1->seq.s << "\n";
					(*outputFiles2[outputFileIndex]) << ">" << rec2->name.s << " " << rec2->comment.s
							<< " " << score2 << "\n" << rec2->seq.s << "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << ">" << rec1->name.s << " " << rec1->comment.s
							<< "\n" << rec1->seq.s << "\n";
					(*outputFiles2[outputFileIndex]) << ">" << rec2->name.s << " " << rec2->comment.s
							<< "\n" << rec2->seq.s << "\n";
				}
			}
		} else {
			if (m_mode == SCORES
					&& outputFileIndex == rm.getMultiMatchIndex()) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << "@" << rec1->name.s << " " << rec1->comment.s;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles1[outputFileIndex]) << " " << *i;
					}
					(*outputFiles1[outputFileIndex]) << "\n" << rec1->seq.s
							<< "\n+\n" << rec1->qual.s << "\n";
					(*outputFiles2[outputFileIndex]) << "@" << rec2->name.s << " " << rec2->comment.s;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles2[outputFileIndex]) << " " << *i;
					}
					(*outputFiles2[outputFileIndex]) << "\n" << rec2->seq.s
							<< "\n+\n" << rec2->qual.s << "\n";
				}
			} else if (m_mode == BESTHIT) {
				if (outputFileIndex == rm.getMultiMatchIndex())
#pragma omp critical(outputFiles)
						{
					(*outputFiles1[outputFileIndex]) << "@" << rec1->name.s << " " << rec1->comment.s;
					for (vector<double>::iterator i = scores1.begin();
							i != scores1.end(); ++i) {
						(*outputFiles1[outputFileIndex]) << " " << *i;
					}
					(*outputFiles1[outputFileIndex]) << "\n" << rec1->seq.s
							<< "\n+\n" << rec1->qual.s << "\n";
					(*outputFiles2[outputFileIndex]) << "@" << rec2->name.s << " " << rec2->comment.s;
					for (vector<double>::iterator i = scores2.begin();
							i != scores2.end(); ++i) {
						(*outputFiles2[outputFileIndex]) << " " << *i;
					}
					(*outputFiles2[outputFileIndex]) << "\n" << rec2->seq.s
							<< "\n+\n" << rec2->qual.s << "\n";
				} else
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << "@" << rec1->name.s << " " << rec1->comment.s
							<< " " << score1 << "\n" << rec1->seq.s << "\n+\n"
							<< rec1->qual.s << "\n";
					(*outputFiles2[outputFileIndex]) << "@" << rec2->name.s << " " << rec2->comment.s
							<< " " << score2 << "\n" << rec2->seq.s << "\n+\n"
							<< rec2->qual.s << "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles1[outputFileIndex]) << "@" << rec1->name.s << " " << rec1->comment.s
							<< "\n" << rec1->seq.s << "\n+\n" << rec1->qual.s
							<< "\n";
					(*outputFiles2[outputFileIndex]) << "@" << rec2->name.s << " " << rec2->comment.s
							<< "\n" << rec2->seq.s << "\n+\n" << rec2->qual.s
							<< "\n";
				}
			}
		}
	}

	inline void evaluateRead(const string &rec, vector<unsigned> &hits,
			double &score, vector<double> &scores) {
		switch (m_mode) {
		case ORDERED: {
			evaluateReadOrdered(rec, hits);
			break;
		}
//		case MINHITONLY: {
//			evaluateReadMin(rec, hits);
//			break;
//		}
		case BESTHIT: {
			score = evaluateReadBestHit(rec, hits, scores);
			break;
		}
		case SCORES: {
			evaluateReadScore(rec, hits, scores);
			break;
		}
		default: {
			evaluateReadStd(rec, hits);
			break;
		}
		}
	}

	inline void evaluateReadPair(const string &rec1, const string &rec2,
			vector<unsigned> &hits1, vector<unsigned> &hits2, double &score1,
			double &score2, vector<double> &scores1, vector<double> &scores2) {
		switch (m_mode) {
		case ORDERED: {
			evaluateReadOrderedPair(rec1, rec2, hits1, hits2);
			break;
		}
		default: {
			evaluateRead(rec1, hits1, score1, scores1);
			evaluateRead(rec2, hits2, score2, scores2);
			break;
		}
		}
	}
};

#endif /* BIOBLOOMCLASSIFIER_H_ */
