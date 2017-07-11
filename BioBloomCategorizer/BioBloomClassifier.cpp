/*
 * BioBloomClassifier.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: cjustin
 */

#include "BioBloomClassifier.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "ResultsManager.h"
#include "Common/Options.h"
#include <map>
#if _OPENMP
# include <omp.h>
#endif

BioBloomClassifier::BioBloomClassifier(const vector<string> &filterFilePaths,
		double scoreThreshold, const string &prefix,
		const string &outputPostFix, unsigned minHit, bool minHitOnly,
		bool withScore) :
		m_scoreThreshold(scoreThreshold), m_filterNum(filterFilePaths.size()), m_prefix(
				prefix), m_postfix(outputPostFix), m_minHit(minHit), m_mode(
				STD), m_mainFilter(""), m_inclusive(false) {
	loadFilters(filterFilePaths);
	if (minHitOnly && withScore) {
		cerr << "minHit, withScore cannot be used together" << endl;
		exit(1);
	}
	if (minHitOnly) {
		m_mode = MINHITONLY;
	} else if (withScore) {
		m_mode = SCORES;
	}
	if (m_scoreThreshold == 1) {
		m_mode = BESTHIT;
		//TODO: make it possible later?
		//best hit will not allow for more than one hash sig.
		assert(m_mode == BESTHIT && m_hashSigs.size() == 1);
	}
}

/*
 * Generic filtering function (single end, no fa or fq file outputs)
 */
void BioBloomClassifier::filter(const vector<string> &inputFiles) {

	//results summary object
	ResultsManager resSummary(m_filterOrder, m_inclusive);

	size_t totalReads = 0;

	//print out header info and initialize variables

	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it) {
		gzFile fp;
		fp = gzopen(it->c_str(), "r");
		if (fp == Z_NULL) {
			cerr << "file " << *it << " cannot be opened" << endl;
			exit(1);
		}
		kseq_t *kseq = kseq_init(fp);
		FaRec rec;

#pragma omp parallel private(rec)
		for (int l;;) {
#pragma omp critical(kseq_read)
			{
				l = kseq_read(kseq);
				if (l >= 0) {
					rec.seq = string(kseq->seq.s, l);
					rec.header = string(kseq->name.s, kseq->name.l);
					rec.qual = string(kseq->qual.s, kseq->qual.l);
				}
			}
			if (l >= 0) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				unordered_map<string, bool> hits(m_filterNum);
				double score = 0; //Todo: figure out what happens to this if multiple hashSigs are used
				vector<double> scores(m_filterNum, 0.0);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = m_hashSigs.begin();
						j != m_hashSigs.end(); ++j) {
					evaluateRead(rec.seq, *j, hits, score, scores);
				}

				//Evaluate hit data and record for summary and print if needed
				printSingle(rec, score, resSummary.updateSummaryData(hits));

			} else
				break;
		}
		kseq_destroy(kseq);
		gzclose(fp);
	}

	cerr << "Total Reads:" << totalReads << endl;

	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

/*
 * Filters reads
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * Prints reads into separate files
 */
void BioBloomClassifier::filterPrint(const vector<string> &inputFiles,
		const string &outputType) {

	//results summary object
	ResultsManager resSummary(m_filterOrder, m_inclusive);

	size_t totalReads = 0;

	unordered_map<string, boost::shared_ptr<Dynamicofstream> > outputFiles;
	boost::shared_ptr<Dynamicofstream> no_match(
			new Dynamicofstream(
					m_prefix + "_" + NO_MATCH + "." + outputType + m_postfix));
	boost::shared_ptr<Dynamicofstream> multi_match(
			new Dynamicofstream(
					m_prefix + "_" + MULTI_MATCH + "." + outputType
							+ m_postfix));
	outputFiles[NO_MATCH] = no_match;
	outputFiles[MULTI_MATCH] = multi_match;

	//initialize variables
	for (vector<string>::const_iterator j = m_hashSigs.begin();
			j != m_hashSigs.end(); ++j) {
		const vector<string> idsInFilter = (*m_filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			boost::shared_ptr<Dynamicofstream> temp(
					new Dynamicofstream(
							m_prefix + "_" + *i + "." + outputType
									+ m_postfix));
			outputFiles[*i] = temp;
		}
	}

	//print out header info and initialize variables

	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it) {
		gzFile fp;
		fp = gzopen(it->c_str(), "r");
		if (fp == Z_NULL) {
			cerr << "file " << *it << " cannot be opened" << endl;
			exit(1);
		}

		kseq_t *kseq = kseq_init(fp);
		FaRec rec;

#pragma omp parallel private(rec)
		for (int l;;) {
#pragma omp critical(kseq_read)
			{
				l = kseq_read(kseq);
				if (l >= 0) {
					rec.seq = string(kseq->seq.s, l);
					rec.header = string(kseq->name.s, kseq->name.l);
					rec.qual = string(kseq->qual.s, kseq->qual.l);
				}
			}
			if (l >= 0) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				unordered_map<string, bool> hits(m_filterNum);
				double score = 0.0;
				vector<double> scores(m_filterNum, 0.0);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = m_hashSigs.begin();
						j != m_hashSigs.end(); ++j) {
					evaluateRead(rec.seq, *j, hits, score, scores);
				}

				//Evaluate hit data and record for summary
				const string &outputFileName = resSummary.updateSummaryData(
						hits);

				printSingle(rec, score, outputFileName);

				printSingleToFile(outputFileName, rec, outputFiles, outputType,
						score, scores);

			} else
				break;
		}
		kseq_destroy(kseq);
		gzclose(fp);
	}

	//close sorting files
	for (unordered_map<string, boost::shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j) {
		j->second->close();
		cerr << "File written to: "
				<< m_prefix + "_" + j->first + "." + outputType + m_postfix
				<< endl;
	}
	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPair(const string &file1, const string &file2) {

	//results summary object
	ResultsManager resSummary(m_filterOrder, m_inclusive);

	size_t totalReads = 0;

	cerr << "Filtering Start" << "\n";

	gzFile fp1, fp2;
	fp1 = gzopen(file1.c_str(), "r");
	if (fp1 == Z_NULL) {
		cerr << "file " << file1.c_str() << " cannot be opened" << endl;
		exit(1);
	}
	fp2 = gzopen(file2.c_str(), "r");
	if (fp2 == Z_NULL) {
		cerr << "file " << file2.c_str() << " cannot be opened" << endl;
		exit(1);
	}
	kseq_t *kseq1 = kseq_init(fp1);
	kseq_t *kseq2 = kseq_init(fp2);
	FaRec rec1;
	FaRec rec2;

#pragma omp parallel private(rec1, rec2)
	for (int l1, l2;;) {
#pragma omp critical(kseq)
		{
			l1 = kseq_read(kseq1);
			if (l1 >= 0) {
				rec1.seq = string(kseq1->seq.s, l1);
				rec1.header = string(kseq1->name.s, kseq1->name.l);
				rec1.qual = string(kseq1->qual.s, kseq1->qual.l);
			}
			l2 = kseq_read(kseq2);
			if (l2 >= 0) {
				rec2.seq = string(kseq2->seq.s, l2);
				rec2.header = string(kseq2->name.s, kseq2->name.l);
				rec2.qual = string(kseq2->qual.s, kseq2->qual.l);
			}
		}
		if (l1 >= 0 && l2 >= 0) {
#pragma omp critical(totalReads)
			{
				++totalReads;
				if (totalReads % 10000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}
			}

			//hits results stored in hashmap of filter names and hits
			unordered_map<string, bool> hits1(m_filterNum);
			unordered_map<string, bool> hits2(m_filterNum);

			double score1 = 0;
			double score2 = 0;

			vector<double> scores1(m_filterNum, 0.0);
			vector<double> scores2(m_filterNum, 0.0);

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = m_hashSigs.begin();
					j != m_hashSigs.end(); ++j) {
				evaluateReadPair(rec1.seq, rec2.seq, *j, hits1, hits2, score1,
						score2, scores1, scores2);
			}

			//Evaluate hit data and record for summary
			printPair(rec1, rec2, score1, score2,
					resSummary.updateSummaryData(hits1, hits2));
		} else
			break;
	}
	kseq_destroy(kseq1);
	kseq_destroy(kseq2);
	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

/*
 * Filtering using sets of files
 */
void BioBloomClassifier::filterPair(const vector<string> &inputFiles1,
		const vector<string> &inputFiles2) {

	//results summary object
	ResultsManager resSummary(m_filterOrder, m_inclusive);

	size_t totalReads = 0;

	cerr << "Filtering Start" << "\n";

#pragma omp parallel for
	for (unsigned i = 0; i < inputFiles1.size(); ++i) {
		gzFile fp1, fp2;
		fp1 = gzopen(inputFiles1[i].c_str(), "r");
		if (fp1 == Z_NULL) {
			cerr << "file " << inputFiles1[i].c_str() << " cannot be opened"
					<< endl;
			exit(1);
		}
		fp2 = gzopen(inputFiles2[i].c_str(), "r");
		if (fp2 == Z_NULL) {
			cerr << "file " << inputFiles2[i].c_str() << " cannot be opened"
					<< endl;
			exit(1);
		}
		kseq_t *kseq1 = kseq_init(fp1);
		kseq_t *kseq2 = kseq_init(fp2);

		for (int l1, l2;;) {
			l1 = kseq_read(kseq1);
			l2 = kseq_read(kseq2);
			if (l1 >= 0 && l2 >= 0) {
#pragma omp atomic
				++totalReads;
#pragma omp critical(totalReads)
				{
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				//hits results stored in hashmap of filter names and hits
				unordered_map<string, bool> hits1(m_filterNum);
				unordered_map<string, bool> hits2(m_filterNum);

				double score1 = 0;
				double score2 = 0;

				vector<double> scores1(m_filterNum, 0.0);
				vector<double> scores2(m_filterNum, 0.0);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = m_hashSigs.begin();
						j != m_hashSigs.end(); ++j) {
					evaluateReadPair(kseq1->seq.s, kseq2->seq.s, *j, hits1,
							hits2, score1, score2, scores1, scores2);
				}

				//Evaluate hit data and record for summary
				const string &outputFileName = resSummary.updateSummaryData(
						hits1, hits2);
				printPair(kseq1, kseq2, score1, score2, outputFileName);
			} else
				break;
		}
		kseq_destroy(kseq1);
		kseq_destroy(kseq2);
	}
	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * prints reads
 */
void BioBloomClassifier::filterPairPrint(const string &file1,
		const string &file2, const string &outputType) {

	//results summary object
	ResultsManager resSummary(m_filterOrder, m_inclusive);

	size_t totalReads = 0;

	unordered_map<string, boost::shared_ptr<Dynamicofstream> > outputFiles;
	boost::shared_ptr<Dynamicofstream> noMatch1(
			new Dynamicofstream(
					m_prefix + "_" + NO_MATCH + "_1." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> noMatch2(
			new Dynamicofstream(
					m_prefix + "_" + NO_MATCH + "_2." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> multiMatch1(
			new Dynamicofstream(
					m_prefix + "_" + MULTI_MATCH + "_1." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> multiMatch2(
			new Dynamicofstream(
					m_prefix + "_" + MULTI_MATCH + "_2." + outputType
							+ m_postfix));
	outputFiles[NO_MATCH + "_1"] = noMatch1;
	outputFiles[NO_MATCH + "_2"] = noMatch2;
	outputFiles[MULTI_MATCH + "_1"] = multiMatch1;
	outputFiles[MULTI_MATCH + "_2"] = multiMatch2;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = m_hashSigs.begin();
			j != m_hashSigs.end(); ++j) {
		const vector<string> idsInFilter = (*m_filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			boost::shared_ptr<Dynamicofstream> temp1(
					new Dynamicofstream(
							m_prefix + "_" + *i + "_1." + outputType
									+ m_postfix));
			boost::shared_ptr<Dynamicofstream> temp2(
					new Dynamicofstream(
							m_prefix + "_" + *i + "_2." + outputType
									+ m_postfix));
			outputFiles[*i + "_1"] = temp1;
			outputFiles[*i + "_2"] = temp2;
		}
	}

	//for output files in consistent order

	cerr << "Filtering Start" << "\n";

	gzFile fp1, fp2;
	fp1 = gzopen(file1.c_str(), "r");
	if (fp1 == Z_NULL) {
		cerr << "file " << file1.c_str() << " cannot be opened" << endl;
		exit(1);
	}
	fp2 = gzopen(file2.c_str(), "r");
	if (fp2 == Z_NULL) {
		cerr << "file " << file2.c_str() << " cannot be opened" << endl;
		exit(1);
	}
	kseq_t *kseq1 = kseq_init(fp1);
	kseq_t *kseq2 = kseq_init(fp2);
	FaRec rec1;
	FaRec rec2;

#pragma omp parallel private(rec1, rec2)
	for (int l1, l2;;) {
#pragma omp critical(kseq)
		{
			l1 = kseq_read(kseq1);
			if (l1 >= 0) {
				rec1.seq = string(kseq1->seq.s, l1);
				rec1.header = string(kseq1->name.s, kseq1->name.l);
				rec1.qual = string(kseq1->qual.s, kseq1->qual.l);
			}
			l2 = kseq_read(kseq2);
			if (l2 >= 0) {
				rec2.seq = string(kseq2->seq.s, l2);
				rec2.header = string(kseq2->name.s, kseq2->name.l);
				rec2.qual = string(kseq2->qual.s, kseq2->qual.l);
			}
		}
		if (l1 >= 0 && l2 >= 0) {
#pragma omp critical(totalReads)
			{
				++totalReads;
				if (totalReads % 10000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}
			}

			//hits results stored in hashmap of filter names and hits
			unordered_map<string, bool> hits1(m_filterNum);
			unordered_map<string, bool> hits2(m_filterNum);

			double score1 = 0;
			double score2 = 0;

			vector<double> scores1(m_filterNum, 0.0);
			vector<double> scores2(m_filterNum, 0.0);

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = m_hashSigs.begin();
					j != m_hashSigs.end(); ++j) {
				evaluateReadPair(rec1.seq, rec2.seq, *j, hits1, hits2, score1,
						score2, scores1, scores2);
			}

			//Evaluate hit data and record for summary

			const string &outputFileName = resSummary.updateSummaryData(hits1,
					hits2);
			printPair(rec1, rec2, score1, score2, outputFileName);
			printPairToFile(outputFileName, rec1, rec2, outputFiles, outputType,
					score1, score2, scores1, scores2);
		} else
			break;
	}
	kseq_destroy(kseq1);
	kseq_destroy(kseq2);

	//close sorting files
	for (unordered_map<string, boost::shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j) {
		j->second->close();
		cerr << "File written to: "
				<< m_prefix + "_" + j->first + "." + outputType + m_postfix
				<< endl;
	}

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPair(const string &file) {

	//results summary object
	ResultsManager resSummary(m_filterOrder, m_inclusive);

	unordered_map<string, FaRec> unPairedReads;

	size_t totalReads = 0;

	//print out header info and initialize variables for summary

	cerr << "Filtering Start" << "\n";

	gzFile fp = gzopen(file.c_str(), "r");
	if (fp == Z_NULL) {
		cerr << "file " << file.c_str() << " cannot be opened" << endl;
		exit(1);
	}

	kseq_t *kseq = kseq_init(fp);
	FaRec rec;

#pragma omp parallel private(rec)
	for (int l;;) {
#pragma omp critical(kseq_read)
		{
			l = kseq_read(kseq);
			if (l >= 0) {
				rec.header = string(kseq->name.s, kseq->name.l);
				rec.seq = string(kseq->seq.s, l);
				rec.qual = string(kseq->qual.s, kseq->qual.l);
			}
		}
		bool pairFound;
		if (l >= 0) {
			//TODO:prevent copy somehow?
			string readID = string(rec.header, rec.header.length() - 2);
#pragma omp critical(unPairedReads)
			{
				if (unPairedReads.find(readID) != unPairedReads.end()) {
					pairFound = true;
				} else {
					unPairedReads[readID] = rec;
				}
			}
			if (pairFound) {
				FaRec &rec1 =
						rec.header[rec.header.length() - 1] == '1' ?
								rec : unPairedReads[readID];
				FaRec &rec2 =
						rec.header[rec.header.length() - 1] == '1' ?
								unPairedReads[readID] : rec;
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				unordered_map<string, bool> hits1(m_filterNum);
				unordered_map<string, bool> hits2(m_filterNum);

				double score1 = 0;
				double score2 = 0;

				vector<double> scores1(m_filterNum, 0.0);
				vector<double> scores2(m_filterNum, 0.0);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = m_hashSigs.begin();
						j != m_hashSigs.end(); ++j) {
					evaluateReadPair(rec1.seq, rec2.seq, *j, hits1, hits2,
							score1, score2, scores1, scores2);
				}

				//Evaluate hit data and record for summary
				const string &outputFileName = resSummary.updateSummaryData(
						hits1, hits2);
				printPair(rec1, rec2, score1, score2, outputFileName);
				unPairedReads.erase(readID);
			}
		} else
			break;
	}
	kseq_destroy(kseq);

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * Prints reads into separate files
 */
void BioBloomClassifier::filterPairPrint(const string &file,
		const string &outputType) {

	//results summary object
	ResultsManager resSummary(m_filterOrder, m_inclusive);

	unordered_map<string, FaRec> unPairedReads;

	size_t totalReads = 0;

	unordered_map<string, boost::shared_ptr<Dynamicofstream> > outputFiles;
	boost::shared_ptr<Dynamicofstream> noMatch1(
			new Dynamicofstream(
					m_prefix + "_" + NO_MATCH + "_1." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> noMatch2(
			new Dynamicofstream(
					m_prefix + "_" + NO_MATCH + "_2." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> multiMatch1(
			new Dynamicofstream(
					m_prefix + "_" + MULTI_MATCH + "_1." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> multiMatch2(
			new Dynamicofstream(
					m_prefix + "_" + MULTI_MATCH + "_2." + outputType
							+ m_postfix));
	outputFiles[NO_MATCH + "_1"] = noMatch1;
	outputFiles[NO_MATCH + "_2"] = noMatch2;
	outputFiles[MULTI_MATCH + "_1"] = multiMatch1;
	outputFiles[MULTI_MATCH + "_2"] = multiMatch2;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = m_hashSigs.begin();
			j != m_hashSigs.end(); ++j) {
		const vector<string> idsInFilter = (*m_filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			boost::shared_ptr<Dynamicofstream> temp1(
					new Dynamicofstream(
							m_prefix + "_" + *i + "_1." + outputType
									+ m_postfix));
			boost::shared_ptr<Dynamicofstream> temp2(
					new Dynamicofstream(
							m_prefix + "_" + *i + "_2." + outputType
									+ m_postfix));
			outputFiles[*i + "_1"] = temp1;
			outputFiles[*i + "_2"] = temp2;
		}
	}

	//print out header info and initialize variables for summary
	cerr << "Filtering Start" << "\n";

	gzFile fp = gzopen(file.c_str(), "r");
	if (fp == Z_NULL) {
		cerr << "file " << file.c_str() << " cannot be opened" << endl;
		exit(1);
	}

	kseq_t *kseq = kseq_init(fp);
	FaRec rec;

#pragma omp parallel private(rec)
	for (int l;;) {
#pragma omp critical(kseq_read)
		{
			l = kseq_read(kseq);
			if (l >= 0) {
				rec.header = string(kseq->name.s, kseq->name.l);
				rec.seq = string(kseq->seq.s, l);
				rec.qual = string(kseq->qual.s, kseq->qual.l);
			}
		}
		bool pairFound;
		if (l >= 0) {
			//TODO:prevent copy somehow?
			string readID = string(rec.header, rec.header.length() - 2);
#pragma omp critical(unPairedReads)
			{
				if (unPairedReads.find(readID) != unPairedReads.end()) {
					pairFound = true;
				} else {
					unPairedReads[readID] = rec;
				}
			}
			if (pairFound) {
				FaRec &rec1 =
						rec.header[rec.header.length() - 1] == '1' ?
								rec : unPairedReads[readID];
				FaRec &rec2 =
						rec.header[rec.header.length() - 1] == '1' ?
								unPairedReads[readID] : rec;
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				unordered_map<string, bool> hits1(m_filterNum);
				unordered_map<string, bool> hits2(m_filterNum);

				double score1 = 0;
				double score2 = 0;

				vector<double> scores1(m_filterNum, 0.0);
				vector<double> scores2(m_filterNum, 0.0);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = m_hashSigs.begin();
						j != m_hashSigs.end(); ++j) {
					evaluateReadPair(rec1.seq, rec2.seq, *j, hits1, hits2,
							score1, score2, scores1, scores2);
				}

				//Evaluate hit data and record for summary
				const string &outputFileName = resSummary.updateSummaryData(
						hits1, hits2);
				printPairToFile(outputFileName, rec1, rec2, outputFiles,
						outputType, score1, score2, scores1, scores2);
				printPair(rec1, rec2, score1, score2, outputFileName);
				unPairedReads.erase(readID);
			}
		} else
			break;
	}
	kseq_destroy(kseq);
	//close sorting files
	for (unordered_map<string, boost::shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j) {
		j->second->close();
		cerr << "File written to: "
				<< m_prefix + "_" + j->first + "." + outputType + m_postfix
				<< endl;
	}

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

//helper methods

/*
 * Loads list of filters into memory
 * todo: Implement non-block I/O when loading multiple filters at once
 */
void BioBloomClassifier::loadFilters(const vector<string> &filterFilePaths) {
	cerr << "Starting to Load Filters." << endl;
	//load up files
	for (vector<string>::const_iterator it = filterFilePaths.begin();
			it != filterFilePaths.end(); ++it) {
		//check if files exist
		if (!fexists(*it)) {
			cerr << "Error: " + (*it) + " File cannot be opened" << endl;
			exit(1);
		}
		string infoFileName = (*it).substr(0, (*it).length() - 2) + "txt";
		if (!fexists(infoFileName)) {
			cerr
					<< "Error: " + (infoFileName)
							+ " File cannot be opened. A corresponding info file is needed."
					<< endl;
			exit(1);
		}

		//info file creation
		boost::shared_ptr<BloomFilterInfo> info(
				new BloomFilterInfo(infoFileName));
		//append kmer size to hash signature to insure correct kmer size is used
		stringstream hashSig;
		hashSig << info->getHashNum() << info->getKmerSize();

		//if hashSig exists add filter to list
		if (m_infoFiles.count(hashSig.str()) != 1) {
			m_hashSigs.push_back(hashSig.str());
			vector<boost::shared_ptr<BloomFilterInfo> > tempVect;
			boost::shared_ptr<MultiFilter> temp(
					new MultiFilter(info->getHashNum(), info->getKmerSize()));
			m_filters[hashSig.str()] = temp;
			m_infoFiles[hashSig.str()] = tempVect;
		}
		m_infoFiles[hashSig.str()].push_back(info);
		boost::shared_ptr<BloomFilter> filter(
				new BloomFilter(info->getCalcuatedFilterSize(),
						info->getHashNum(), info->getKmerSize(), *it));
		m_filters[hashSig.str()]->addFilter(info->getFilterID(), filter);
		m_filtersSingle[info->getFilterID()] = filter;
		m_filterOrder.push_back(info->getFilterID());
		cerr << "Loaded Filter: " + info->getFilterID() << endl;
	}
	if (m_scoreThreshold == 1 && m_hashSigs.size() > 1) {
		cerr
				<< "If -s = 1 (best hit mode) all filters must use the same k and same number of hash functions."
				<< endl;
		exit(1);
	}
	cerr << "Filter Loading Complete." << endl;
}

/*
 * checks if file exists
 */
bool BioBloomClassifier::fexists(const string &filename) const {
	ifstream ifile(filename.c_str());
	return ifile.good();
}

/*
 * Collaborative filtering method
 * Assume filters use the same k-mer size
 */
void BioBloomClassifier::evaluateReadCollab(const string &rec,
		const string &hashSig, unordered_map<string, bool> &hits) {
	//get filterIDs to iterate through has in a consistent order
	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();

	//todo: read proc possibly unneeded, see evalSingle
	ReadsProcessor proc(kmerSize);

	//create storage for hits per filter
	std::multimap<unsigned, string> firstPassHits;

	//base for each filter until one filter obtains hit threshold
	//TODO: staggered pattering
	for (vector<string>::reverse_iterator i = m_filterOrder.rbegin();
			i != m_filterOrder.rend(); ++i) {
		hits[*i] = false;
		unsigned screeningHits = 0;
		size_t screeningLoc = rec.length() % kmerSize / 2;
		//First pass filtering
		while (rec.length() >= screeningLoc + kmerSize) {
			const unsigned char* currentKmer = proc.prepSeq(rec, screeningLoc);
			if (currentKmer != NULL) {
				if (m_filtersSingle.at(*i)->contains(currentKmer)) {
					++screeningHits;
				}
			}
			screeningLoc += kmerSize;
		}
		firstPassHits.insert(pair<unsigned, string>(screeningHits, *i));
	}

	//evaluate promising group first
	for (multimap<unsigned, string>::reverse_iterator i =
			firstPassHits.rbegin(); i != firstPassHits.rend(); ++i) {
		string filterID = i->second;
		BloomFilter &tempFilter = *m_filtersSingle.at(filterID);
		if (SeqEval::evalRead(rec, kmerSize, tempFilter, m_scoreThreshold,
				1.0 - m_scoreThreshold, getEvalMode())) {
			hits[filterID] = true;
			break;
		}
	}
}

/*
 * Collaborative filtering method
 * Assume filters use the same k-mer size
 */
void BioBloomClassifier::evaluateReadCollabPair(const string &rec1,
		const string &rec2, const string &hashSig,
		unordered_map<string, bool> &hits1,
		unordered_map<string, bool> &hits2) {
	//get filterIDs to iterate through has in a consistent order
	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();

	//todo: read proc possibly unneeded, see evalSingle
	ReadsProcessor proc(kmerSize);

	//create storage for hits per filter
	std::multimap<unsigned, string> firstPassHits;

	//base for each filter until one filter obtains hit threshold
	//TODO: staggered pattering
	for (vector<string>::reverse_iterator i = m_filterOrder.rbegin();
			i != m_filterOrder.rend(); ++i) {
		hits1[*i] = false;
		hits2[*i] = false;
		unsigned screeningHits = 0;
		size_t screeningLoc = rec1.length() % kmerSize / 2;
		//First pass filtering
		while (rec1.length() >= screeningLoc + kmerSize) {
			const unsigned char* currentKmer1 = proc.prepSeq(rec1,
					screeningLoc);
			const unsigned char* currentKmer2 = proc.prepSeq(rec2,
					screeningLoc);
			if (currentKmer1 != NULL) {
				if (m_filtersSingle.at(*i)->contains(currentKmer1)) {
					++screeningHits;
				}
			}
			if (currentKmer2 != NULL) {
				if (m_filtersSingle.at(*i)->contains(currentKmer2)) {
					++screeningHits;
				}
			}
			screeningLoc += kmerSize;
		}
		firstPassHits.insert(pair<unsigned, string>(screeningHits, *i));
	}

	//evaluate promising group first
	for (multimap<unsigned, string>::reverse_iterator i =
			firstPassHits.rbegin(); i != firstPassHits.rend(); ++i) {
		BloomFilter &tempFilter = *m_filtersSingle.at(i->second);
		if (m_inclusive) {
			if (SeqEval::evalRead(rec1, kmerSize, tempFilter, m_scoreThreshold,
					1.0 - m_scoreThreshold, getEvalMode())
					|| SeqEval::evalRead(rec2, kmerSize, tempFilter,
							m_scoreThreshold, 1.0 - m_scoreThreshold,
							getEvalMode())) {
				hits1[i->second] = true;
				hits2[i->second] = true;
				break;
			}
		} else {
			if (SeqEval::evalRead(rec1, kmerSize, tempFilter, m_scoreThreshold,
					1.0 - m_scoreThreshold, getEvalMode())
					&& SeqEval::evalRead(rec2, kmerSize, tempFilter,
							m_scoreThreshold, 1.0 - m_scoreThreshold,
							getEvalMode())) {
				hits1[i->second] = true;
				hits2[i->second] = true;
				break;
			}
		}
	}
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Updates hits value to number of hits (hashSig is used to as key)
 * Faster variant that assume there a redundant tile of 0
 */
void BioBloomClassifier::evaluateReadMin(const string &rec,
		const string &hashSig, unordered_map<string, bool> &hits) {
	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*m_filters[hashSig]).getFilterIds();

	//get kmersize for set of info files
	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();

	unordered_map<string, unsigned> tempHits;

	//Establish tiling pattern
	unsigned startModifier1 = (rec.length() % kmerSize) / 2;
	size_t currentKmerNum = 0;

	for (vector<string>::const_iterator i = idsInFilter.begin();
			i != idsInFilter.end(); ++i) {
		tempHits[*i] = 0;
	}

	ReadsProcessor proc(kmerSize);
	//cut read into kmer size given
	while (rec.length() >= (currentKmerNum + 1) * kmerSize) {

		const unsigned char* currentKmer = proc.prepSeq(rec,
				currentKmerNum * kmerSize + startModifier1);

		//check to see if string is invalid
		if (currentKmer != NULL) {

			const unordered_map<string, bool> &results =
					m_filters[hashSig]->multiContains(currentKmer);

			//record hit number in order
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i) {
				if (results.at(*i)) {
					++tempHits[*i];
				}
			}
		}
		++currentKmerNum;
	}
	for (vector<string>::const_iterator i = idsInFilter.begin();
			i != idsInFilter.end(); ++i) {
		hits[*i] = tempHits.at(*i) >= m_minHit;
	}
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 */
void BioBloomClassifier::evaluateReadStd(const string &rec,
		const string &hashSig, unordered_map<string, bool> &hits) {

	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*m_filters[hashSig]).getFilterIds();

	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();

	//todo: read proc possibly unneeded, see evalSingle
	ReadsProcessor proc(kmerSize);

	for (vector<string>::const_iterator i = idsInFilter.begin();
			i != idsInFilter.end(); ++i) {
		bool pass = false;
		hits[*i] = false;
		if (m_minHit > 0) {
			unsigned screeningHits = 0;
			size_t screeningLoc = rec.length() % kmerSize / 2;
			//First pass filtering
			while (rec.length() >= screeningLoc + kmerSize) {
				const unsigned char* currentKmer = proc.prepSeq(rec,
						screeningLoc);
				if (currentKmer != NULL) {
					if (m_filtersSingle.at(*i)->contains(currentKmer)) {
						screeningHits++;
						if (screeningHits >= m_minHit) {
							pass = true;
							break;
						}
					}
				}
				screeningLoc += kmerSize;
			}
		} else {
			pass = true;
		}
		if (pass) {
			BloomFilter &tempFilter = *m_filtersSingle.at(*i);
			hits[*i] = SeqEval::evalRead(rec, kmerSize, tempFilter,
					m_scoreThreshold, 1.0 - m_scoreThreshold, getEvalMode());
		}
	}
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Reads are assigned to best hit
 */
double BioBloomClassifier::evaluateReadBestHit(const string &rec,
		const string &hashSig, unordered_map<string, bool> &hits,
		vector<double> &scores) {
	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*m_filters[hashSig]).getFilterIds();

	vector<string> bestFilters;
	double maxScore = 0;

	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();

	ReadsProcessor proc(kmerSize);

	for (unsigned i = 0; i < idsInFilter.size(); ++i) {
		bool pass = false;
		hits[idsInFilter[i]] = false;
		if (m_minHit > 0) {
			unsigned screeningHits = 0;
			size_t screeningLoc = rec.length() % kmerSize / 2;
			//First pass filtering
			while (rec.length() >= screeningLoc + kmerSize) {
				const unsigned char* currentKmer = proc.prepSeq(rec,
						screeningLoc);
				if (currentKmer != NULL) {
					if (m_filtersSingle.at(idsInFilter[i])->contains(
							currentKmer)) {
						screeningHits++;
						if (screeningHits >= m_minHit) {
							pass = true;
							break;
						}
					}
				}
				screeningLoc += kmerSize;
			}
		} else {
			pass = true;
		}
		if (pass) {
			BloomFilter &tempFilter = *m_filtersSingle.at(idsInFilter[i]);
			double score = SeqEval::evalSingleExhaust(rec, kmerSize,
					tempFilter);
			scores[i] = SeqEval::normalizeScore(score, kmerSize, rec.length());
			;
			if (maxScore < score) {
				maxScore = score;
				bestFilters.clear();
				bestFilters.push_back(idsInFilter[i]);
			} else if (maxScore == score) {
				bestFilters.push_back(idsInFilter[i]);
			}
		}
	}
	if (maxScore > 0) {
		for (unsigned i = 0; i < bestFilters.size(); ++i) {
			hits[bestFilters[i]] = true;
		}
	}
	return SeqEval::normalizeScore(maxScore, kmerSize, rec.length());
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Will return scores in vector
 * Will return partial score if threshold is not met
 */
void BioBloomClassifier::evaluateReadScore(const string &rec,
		const string &hashSig, unordered_map<string, bool> &hits,
		vector<double> &scores) {
	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*m_filters[hashSig]).getFilterIds();

	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();

	//todo: read proc possibly unneeded, see evalSingle
	ReadsProcessor proc(kmerSize);

	size_t numKmers = rec.length() - kmerSize + 1;
	unsigned hitCount = 0;

	vector<vector<size_t> > hashValues(numKmers);
	vector<bool> visited(numKmers);

	//position of sequences
	vector<unsigned> pos(idsInFilter.size(), 0);

	//first pass
	for (unsigned i = 0; i < idsInFilter.size(); ++i) {
		bool pass = false;
		hits[idsInFilter[i]] = false;
		if (m_minHit > 0) {
			unsigned screeningHits = 0;
			size_t screeningLoc = rec.length() % kmerSize / 2;
			//First pass filtering
			while (rec.length() >= screeningLoc + kmerSize) {
				const unsigned char* currentKmer = proc.prepSeq(rec,
						screeningLoc);
				if (currentKmer != NULL) {
					if (m_filtersSingle.at(idsInFilter[i])->contains(
							currentKmer)) {
						screeningHits++;
						if (screeningHits >= m_minHit) {
							pass = true;
							break;
						}
					}
				}
				screeningLoc += kmerSize;
			}
		} else {
			pass = true;
		}
		if (pass) {
			BloomFilter &tempFilter = *m_filtersSingle.at(idsInFilter[i]);

			//Evaluate sequences until threshold
			//record end location
			hits[idsInFilter[i]] = SeqEval::eval(rec, kmerSize, tempFilter,
					m_scoreThreshold, 1.0 - m_scoreThreshold, visited,
					hashValues, pos[i], scores[i], proc);
			hitCount += hits[idsInFilter[i]];
		}
	}

	//final pass if more than 2 reach threshold
	if (hitCount > 1) {
		for (unsigned i = 0; i < idsInFilter.size(); ++i) {
			BloomFilter &tempFilter = *m_filtersSingle.at(idsInFilter[i]);

			//Evaluate sequences until threshold
			//record end location
			SeqEval::eval(rec, kmerSize, tempFilter, 1, 0, visited, hashValues,
					pos[i], scores[i], proc);
		}
	}
}

void BioBloomClassifier::setMainFilter(const string &filtername) {
	if (m_filtersSingle.find(filtername) == m_filtersSingle.end()) {
		cerr << "Filter with this name \"" << filtername
				<< "\" does not exist\n";
		cerr << "Valid filter Names:\n";
		for (vector<string>::const_iterator itr = m_filterOrder.begin();
				itr != m_filterOrder.end(); ++itr) {
			cerr << *itr << endl;
		}
		exit(1);
	}
	m_mainFilter = filtername;
}

BioBloomClassifier::~BioBloomClassifier() {
}

