/*
 * BioBloomClassifier.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: cjustin
 */

#include <ResultsManager.hpp>
#include "BioBloomClassifier.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "Common/Options.h"
#include <boost/shared_ptr.hpp>
#include <unordered_map>
#if _OPENMP
# include <omp.h>
#endif

BioBloomClassifier::BioBloomClassifier(const vector<string> &filterFilePaths,
		double scoreThreshold, const string &prefix,
		const string &outputPostFix) :
		m_scoreThreshold(scoreThreshold), m_filterNum(filterFilePaths.size()), m_prefix(
				prefix), m_postfix(outputPostFix), m_stdout(false), m_inclusive(
				false) {
	loadFilters(filterFilePaths);
}

/*
 * Generic filtering function (single end, no fa or fq file outputs)
 */
void BioBloomClassifier::filter(const vector<string> &inputFiles) {

	//results summary object
	ResultsManager<unsigned> resSummary(m_filterOrder, m_inclusive);

	size_t totalReads = 0;

	//print out header info and initialize variables

	cerr << "Filtering Start" << endl;

	FaRec rec;
	vector<double> scores(m_filterNum, 0);
	vector<unsigned> hits;
	hits.reserve(m_filterNum);
	double score = 0;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it) {
		gzFile fp;
		fp = gzopen(it->c_str(), "r");
		if (fp == Z_NULL) {
			cerr << "file " << *it << " cannot be opened" << endl;
			exit(1);
		}
		kseq_t *kseq = kseq_init(fp);
#pragma omp parallel private(rec, scores, score, hits)
		for (int l;;) {
#pragma omp critical(kseq_read)
			{
				l = kseq_read(kseq);
				if (l >= 0) {
					rec.seq = string(kseq->seq.s, l);
					rec.header = string(kseq->name.s, kseq->name.l);
					rec.qual = string(kseq->qual.s, kseq->qual.l);
					rec.comment = string(kseq->comment.s, kseq->comment.l);
				}
			}
			if (l >= 0) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				hits.clear();
				score = 0;
				scores.clear();
				evaluateRead(rec.seq, hits, score, scores);
				//Evaluate hit data and record for summary and print if needed
				printSingle(rec, score, resSummary.updateSummaryData(hits));

			} else
				break;
		}
		kseq_destroy(kseq);
		gzclose(fp);
	}
	cerr << "Total Reads: " << totalReads << "\n";
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
	ResultsManager<unsigned> resSummary(m_filterOrder, m_inclusive);

	size_t totalReads = 0;

	vector<Dynamicofstream*> outputFiles(m_filterOrder.size() + 2, 0);
	//initialize variables
	unsigned index = 0;
	for (vector<string>::const_iterator i = m_filterOrder.begin();
			i != m_filterOrder.end(); ++i) {
		outputFiles[index++] = new Dynamicofstream(
				m_prefix + "_" + *i + "." + outputType + m_postfix);
	}
	outputFiles[index++] = new Dynamicofstream(
			m_prefix + "_" + NO_MATCH + "." + outputType + m_postfix);
	outputFiles[index] = new Dynamicofstream(
			m_prefix + "_" + MULTI_MATCH + "." + outputType + m_postfix);

	//print out header info and initialize variables

	cerr << "Filtering Start" << endl;

	FaRec rec;
	vector<double> scores(m_filterNum, 0);
	vector<unsigned> hits;
	hits.reserve(m_filterNum);
	double score = 0;
	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it) {
		gzFile fp;
		fp = gzopen(it->c_str(), "r");
		if (fp == Z_NULL) {
			cerr << "file " << *it << " cannot be opened" << endl;
			exit(1);
		}

		kseq_t *kseq = kseq_init(fp);
#pragma omp parallel private(rec, scores, score, hits)
		for (int l;;) {
#pragma omp critical(kseq_read)
			{
				l = kseq_read(kseq);
				if (l >= 0) {
					rec.seq = string(kseq->seq.s, l);
					rec.header = string(kseq->name.s, kseq->name.l);
					rec.qual = string(kseq->qual.s, kseq->qual.l);
					rec.comment = string(kseq->comment.s, kseq->comment.l);
				}
			}
			if (l >= 0) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				hits.clear();
				score = 0;
				scores.clear();
				evaluateRead(rec.seq, hits, score, scores);
				//Evaluate hit data and record for summary
				unsigned outputFileName = resSummary.updateSummaryData(hits);
				printSingle(rec, score, outputFileName);
				printSingleToFile(outputFileName, rec, outputFiles, outputType,
						score, scores, resSummary);

			} else
				break;
		}
		kseq_destroy(kseq);
		gzclose(fp);
	}

	//close sorting files
	for (unsigned i = 0; i < m_filterOrder.size(); ++i) {
		outputFiles[i]->close();
		delete (outputFiles[i]);
		cerr << "File written to: "
				<< m_prefix + "_" + m_filterOrder[i] + "." + outputType
						+ m_postfix << endl;
	}
	outputFiles[m_filterOrder.size()]->close();
	delete (outputFiles[m_filterOrder.size()]);
	cerr << "File written to: "
			<< m_prefix + "_" + NO_MATCH + "." + outputType + m_postfix << endl;
	outputFiles[m_filterOrder.size() + 1]->close();
	delete (outputFiles[m_filterOrder.size() + 1]);
	cerr << "File written to: "
			<< m_prefix + "_" + MULTI_MATCH + "." + outputType + m_postfix
			<< endl;
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
	ResultsManager<unsigned> resSummary(m_filterOrder, m_inclusive);

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
	vector<double> scores1(m_filterNum, 0.0);
	vector<double> scores2(m_filterNum, 0.0);
	vector<unsigned> hits1;
	hits1.reserve(m_filterNum);
	vector<unsigned> hits2;
	hits1.reserve(m_filterNum);
	double score1 = 0;
	double score2 = 0;

#pragma omp parallel private(rec1, rec2, scores1, score1, hits1, scores2, score2, hits2)
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
				rec2.comment = string(kseq2->comment.s, kseq2->comment.l);
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
			hits1.clear();
			hits2.clear();
			score1 = 0;
			score2 = 0;
			scores1.clear();
			scores2.clear();

			evaluateReadPair(rec1.seq, rec2.seq, hits1, hits2, score1, score2,
					scores1, scores2);

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
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * prints reads
 */
void BioBloomClassifier::filterPairPrint(const string &file1,
		const string &file2, const string &outputType) {

	//results summary object
	ResultsManager<unsigned> resSummary(m_filterOrder, m_inclusive);

	size_t totalReads = 0;

	vector<Dynamicofstream*> outputFiles1(m_filterOrder.size() + 2, 0);
	vector<Dynamicofstream*> outputFiles2(m_filterOrder.size() + 2, 0);
	//initialize variables
	unsigned index = 0;
	for (vector<string>::const_iterator i = m_filterOrder.begin();
			i != m_filterOrder.end(); ++i) {
		outputFiles1[index] = new Dynamicofstream(
				m_prefix + "_" + *i + "_1." + outputType + m_postfix);
		outputFiles2[index++] = new Dynamicofstream(
				m_prefix + "_" + *i + "_2." + outputType + m_postfix);
	}
	outputFiles1[index] = new Dynamicofstream(
			m_prefix + "_" + NO_MATCH + "_1." + outputType + m_postfix);
	outputFiles2[index++] = new Dynamicofstream(
			m_prefix + "_" + NO_MATCH + "_2." + outputType + m_postfix);
	outputFiles1[index] = new Dynamicofstream(
			m_prefix + "_" + MULTI_MATCH + "_1." + outputType + m_postfix);
	outputFiles2[index] = new Dynamicofstream(
			m_prefix + "_" + MULTI_MATCH + "_2." + outputType + m_postfix);

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
	vector<double> scores1(m_filterNum, 0.0);
	vector<double> scores2(m_filterNum, 0.0);
	vector<unsigned> hits1;
	hits1.reserve(m_filterNum);
	vector<unsigned> hits2;
	hits1.reserve(m_filterNum);
	double score1 = 0;
	double score2 = 0;

#pragma omp parallel private(rec1, rec2, scores1, score1, hits1, scores2, score2, hits2)
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
				rec2.comment = string(kseq2->comment.s, kseq2->comment.l);
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

			hits1.clear();
			hits2.clear();
			score1 = 0;
			score2 = 0;
			scores1.clear();
			scores2.clear();

			evaluateReadPair(rec1.seq, rec2.seq, hits1, hits2, score1, score2,
					scores1, scores2);

			unsigned outputFileIndex = resSummary.updateSummaryData(hits1,
					hits2);

			//Evaluate hit data and record for summary
			printPair(rec1, rec2, score1, score2, outputFileIndex);
			printPairToFile(outputFileIndex, rec1, rec2, outputFiles1,
					outputFiles2, outputType, score1, score2, scores1, scores2, resSummary);
		} else
			break;
	}
	kseq_destroy(kseq1);
	kseq_destroy(kseq2);

	//close sorting files
	for (unsigned i = 0; i < m_filterOrder.size(); ++i) {
		outputFiles1[i]->close();
		delete (outputFiles1[i]);
		cerr << "File written to: "
				<< m_prefix + "_" + m_filterOrder[i] + "_1." + outputType
						+ m_postfix << endl;
		outputFiles2[i]->close();
		delete (outputFiles2[i]);
		cerr << "File written to: "
				<< m_prefix + "_" + m_filterOrder[i] + "_2." + outputType
						+ m_postfix << endl;
	}
	outputFiles1[m_filterOrder.size()]->close();
	delete (outputFiles1[m_filterOrder.size()]);
	cerr << "File written to: "
			<< m_prefix + "_" + NO_MATCH + "_1." + outputType + m_postfix
			<< endl;
	outputFiles2[m_filterOrder.size()]->close();
	delete (outputFiles2[m_filterOrder.size()]);
	cerr << "File written to: "
			<< m_prefix + "_" + NO_MATCH + "_2." + outputType + m_postfix
			<< endl;
	outputFiles1[m_filterOrder.size() + 1]->close();
	delete (outputFiles1[m_filterOrder.size() + 1]);
	cerr << "File written to: "
			<< m_prefix + "_" + MULTI_MATCH + "_1." + outputType + m_postfix
			<< endl;
	outputFiles2[m_filterOrder.size() + 1]->close();
	delete (outputFiles2[m_filterOrder.size() + 1]);
	cerr << "File written to: "
			<< m_prefix + "_" + MULTI_MATCH + "_2." + outputType + m_postfix
			<< endl;
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
	ResultsManager<unsigned> resSummary(m_filterOrder, m_inclusive);

	unordered_map<string, boost::shared_ptr<FaRec> > unPairedReads;

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
	vector<double> scores1(m_filterNum, 0.0);
	vector<double> scores2(m_filterNum, 0.0);
	vector<unsigned> hits1;
	hits1.reserve(m_filterNum);
	vector<unsigned> hits2;
	hits1.reserve(m_filterNum);
	double score1 = 0;
	double score2 = 0;

#pragma omp parallel private(rec, scores1, scores2, hits1, hits2, score1, score2)
	for (int l;;) {
#pragma omp critical(kseq_read)
		{
			l = kseq_read(kseq);
			if (l >= 0) {
				rec.header = string(kseq->name.s,
						kseq->name.s[kseq->name.l - 1] == '1'
								|| kseq->name.s[kseq->name.l - 1] == '2' ?
								kseq->name.l - 2 : kseq->name.l);
				rec.seq = string(kseq->seq.s, l);
				rec.qual = string(kseq->qual.s, kseq->qual.l);
				rec.comment = string(kseq->comment.s, kseq->comment.l);
			}
		}
		if (l >= 0) {
			bool pairFound = false;
			FaRec rec1;
#pragma omp critical(unPairedReads)
			{
				unordered_map<string, boost::shared_ptr<FaRec> >::iterator itr =
						unPairedReads.find(rec.header);
				if (itr != unPairedReads.end()) {
					pairFound = true;
					rec1 = *itr->second;
				} else {
					unPairedReads[rec.header] = boost::shared_ptr<FaRec>(new FaRec(rec));
				}
			}
			if (pairFound) {
				FaRec &rec2 = rec;
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				hits1.clear();
				hits2.clear();
				score1 = 0;
				score2 = 0;
				scores1.clear();
				scores2.clear();

				evaluateReadPair(rec1.seq, rec2.seq, hits1, hits2, score1,
						score2, scores1, scores2);

				unsigned outputFileIndex = resSummary.updateSummaryData(hits1,
						hits2);

				//Evaluate hit data and record for summary
				printPair(rec1, rec2, score1, score2, outputFileIndex);
#pragma omp critical(unPairedReads)
				unPairedReads.erase(rec.header);
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
	ResultsManager<unsigned> resSummary(m_filterOrder, m_inclusive);

	unordered_map<string, boost::shared_ptr<FaRec> > unPairedReads;
	size_t totalReads = 0;

	vector<Dynamicofstream*> outputFiles1(m_filterOrder.size() + 2, 0);
	vector<Dynamicofstream*> outputFiles2(m_filterOrder.size() + 2, 0);
	//initialize variables
	unsigned index = 0;
	for (vector<string>::const_iterator i = m_filterOrder.begin();
			i != m_filterOrder.end(); ++i) {
		outputFiles1[index] = new Dynamicofstream(
				m_prefix + "_" + *i + "_1." + outputType + m_postfix);
		outputFiles2[index++] = new Dynamicofstream(
				m_prefix + "_" + *i + "_2." + outputType + m_postfix);
	}
	outputFiles1[index] = new Dynamicofstream(
			m_prefix + "_" + NO_MATCH + "_1." + outputType + m_postfix);
	outputFiles2[index++] = new Dynamicofstream(
			m_prefix + "_" + NO_MATCH + "_2." + outputType + m_postfix);
	outputFiles1[index] = new Dynamicofstream(
			m_prefix + "_" + MULTI_MATCH + "_1." + outputType + m_postfix);
	outputFiles2[index] = new Dynamicofstream(
			m_prefix + "_" + MULTI_MATCH + "_2." + outputType + m_postfix);

	//print out header info and initialize variables for summary
	cerr << "Filtering Start" << "\n";

	gzFile fp = gzopen(file.c_str(), "r");
	if (fp == Z_NULL) {
		cerr << "file " << file.c_str() << " cannot be opened" << endl;
		exit(1);
	}

	kseq_t *kseq = kseq_init(fp);
	FaRec rec;
	vector<double> scores1(m_filterNum, 0.0);
	vector<double> scores2(m_filterNum, 0.0);
	vector<unsigned> hits1;
	hits1.reserve(m_filterNum);
	vector<unsigned> hits2;
	hits1.reserve(m_filterNum);
	double score1 = 0;
	double score2 = 0;

#pragma omp parallel private(rec, scores1, scores2, hits1, hits2, score1, score2)
	for (int l;;) {
#pragma omp critical(kseq_read)
		{
			l = kseq_read(kseq);
			if (l >= 0) {
				rec.header = string(kseq->name.s,
						kseq->name.s[kseq->name.l - 1] == '1'
								|| kseq->name.s[kseq->name.l - 1] == '2' ?
								kseq->name.l - 2 : kseq->name.l);
				rec.seq = string(kseq->seq.s, l);
				rec.qual = string(kseq->qual.s, kseq->qual.l);
				rec.comment = string(kseq->comment.s, kseq->comment.l);
			}
		}
		if (l >= 0) {
			bool pairFound = false;
			FaRec rec1;
#pragma omp critical(unPairedReads)
			{
				unordered_map<string, boost::shared_ptr<FaRec> >::iterator itr =
						unPairedReads.find(rec.header);
				if (itr != unPairedReads.end()) {
					pairFound = true;
					rec1 = *itr->second;
				} else {
					unPairedReads[rec.header] = boost::shared_ptr<FaRec>(new FaRec(rec));
				}
			}
			if (pairFound) {
				FaRec &rec2 = rec;
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				hits1.clear();
				hits2.clear();
				score1 = 0;
				score2 = 0;
				scores1.clear();
				scores2.clear();

				evaluateReadPair(rec1.seq, rec2.seq, hits1, hits2, score1,
						score2, scores1, scores2);

				unsigned outputFileIndex = resSummary.updateSummaryData(hits1,
						hits2);

				//Evaluate hit data and record for summary
				printPair(rec1, rec2, score1, score2, outputFileIndex);
				printPairToFile(outputFileIndex, rec1, rec2, outputFiles1,
						outputFiles2, outputType, score1, score2, scores1,
						scores2, resSummary);
#pragma omp critical(unPairedReads)
				unPairedReads.erase(rec.header);
			}
		} else
			break;
	}
	kseq_destroy(kseq);
	//close sorting files
	for (unsigned i = 0; i < m_filterOrder.size(); ++i) {
		outputFiles1[i]->close();
		delete (outputFiles1[i]);
		cerr << "File written to: "
				<< m_prefix + "_" + m_filterOrder[i] + "_1." + outputType
						+ m_postfix << endl;
		outputFiles2[i]->close();
		delete (outputFiles2[i]);
		cerr << "File written to: "
				<< m_prefix + "_" + m_filterOrder[i] + "_2." + outputType
						+ m_postfix << endl;
	}
	outputFiles1[m_filterOrder.size()]->close();
	delete (outputFiles1[m_filterOrder.size()]);
	cerr << "File written to: "
			<< m_prefix + "_" + NO_MATCH + "_1." + outputType + m_postfix
			<< endl;
	outputFiles2[m_filterOrder.size()]->close();
	delete (outputFiles2[m_filterOrder.size()]);
	cerr << "File written to: "
			<< m_prefix + "_" + NO_MATCH + "_2." + outputType + m_postfix
			<< endl;
	outputFiles1[m_filterOrder.size() + 1]->close();
	delete (outputFiles1[m_filterOrder.size() + 1]);
	cerr << "File written to: "
			<< m_prefix + "_" + MULTI_MATCH + "_1." + outputType + m_postfix
			<< endl;
	outputFiles2[m_filterOrder.size() + 1]->close();
	delete (outputFiles2[m_filterOrder.size() + 1]);
	cerr << "File written to: "
			<< m_prefix + "_" + MULTI_MATCH + "_2." + outputType + m_postfix
			<< endl;
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
	ResultsManager<unsigned> resSummary(m_filterOrder, m_inclusive);

	size_t totalReads = 0;

	cerr << "Filtering Start" << "\n";

	vector<double> scores1(m_filterNum, 0.0);
	vector<double> scores2(m_filterNum, 0.0);
	vector<unsigned> hits1;
	hits1.reserve(m_filterNum);
	vector<unsigned> hits2;
	hits1.reserve(m_filterNum);
	double score1 = 0;
	double score2 = 0;

#pragma omp parallel for private(scores1, score1, hits1, scores2, score2, hits2)
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
					if (totalReads % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				hits1.clear();
				hits2.clear();
				score1 = 0;
				score2 = 0;
				scores1.clear();
				scores2.clear();
				evaluateReadPair(kseq1->seq.s, kseq2->seq.s, hits1, hits2,
						score1, score2, scores1, scores2);

				//Evaluate hit data and record for summary
				printPair(kseq1, kseq2, score1, score2,
						resSummary.updateSummaryData(hits1, hits2));
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
 * Filtering using sets of files
 */
void BioBloomClassifier::filterPairPrint(const vector<string> &inputFiles1,
		const vector<string> &inputFiles2, const string &outputType) {

	//results summary object
	ResultsManager<unsigned> resSummary(m_filterOrder, m_inclusive);

	size_t totalReads = 0;

	cerr << "Filtering Start" << "\n";

	vector<Dynamicofstream*> outputFiles1;
	vector<Dynamicofstream*> outputFiles2;
	//initialize variables
	unsigned index = 0;
	for (vector<string>::const_iterator i = m_filterOrder.begin();
			i != m_filterOrder.end(); ++i) {
		outputFiles1[index] = new Dynamicofstream(
				m_prefix + "_" + *i + "_1." + outputType + m_postfix);
		outputFiles2[index++] = new Dynamicofstream(
				m_prefix + "_" + *i + "_2." + outputType + m_postfix);
	}
	outputFiles1[index] = new Dynamicofstream(
			m_prefix + "_" + NO_MATCH + "_1." + outputType + m_postfix);
	outputFiles2[index++] = new Dynamicofstream(
			m_prefix + "_" + NO_MATCH + "_2." + outputType + m_postfix);
	outputFiles1[index] = new Dynamicofstream(
			m_prefix + "_" + MULTI_MATCH + "_1." + outputType + m_postfix);
	outputFiles2[index] = new Dynamicofstream(
			m_prefix + "_" + MULTI_MATCH + "_2." + outputType + m_postfix);

	vector<double> scores1(m_filterNum, 0.0);
	vector<double> scores2(m_filterNum, 0.0);
	vector<unsigned> hits1;
	hits1.reserve(m_filterNum);
	vector<unsigned> hits2;
	hits1.reserve(m_filterNum);
	double score1 = 0;
	double score2 = 0;

#pragma omp parallel for private(scores1, score1, hits1, scores2, score2, hits2)
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
					if (totalReads % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				hits1.clear();
				hits2.clear();
				score1 = 0;
				score2 = 0;
				scores1.clear();
				scores2.clear();
				evaluateReadPair(kseq1->seq.s, kseq2->seq.s, hits1, hits2,
						score1, score2, scores1, scores2);

				//Evaluate hit data and record for summary
				unsigned outputFileIndex = resSummary.updateSummaryData(hits1,
						hits2);
				printPair(kseq1, kseq2, score1, score2, outputFileIndex);
				printPairToFile(outputFileIndex, kseq1, kseq2, outputFiles1,
						outputFiles2, outputType, score1, score2, scores1,
						scores2, resSummary);
			} else
				break;
		}
		kseq_destroy(kseq1);
		kseq_destroy(kseq2);
	}
	for (unsigned i = 0; i < m_filterOrder.size(); ++i) {
		outputFiles1[i]->close();
		delete (outputFiles1[i]);
		cerr << "File written to: "
				<< m_prefix + "_" + m_filterOrder[i] + "_1." + outputType
						+ m_postfix << endl;
		outputFiles2[i]->close();
		delete (outputFiles2[i]);
		cerr << "File written to: "
				<< m_prefix + "_" + m_filterOrder[i] + "_2." + outputType
						+ m_postfix << endl;
	}
	outputFiles1[m_filterOrder.size()]->close();
	delete (outputFiles1[m_filterOrder.size()]);
	cerr << "File written to: "
			<< m_prefix + "_" + NO_MATCH + "_1." + outputType + m_postfix
			<< endl;
	outputFiles2[m_filterOrder.size()]->close();
	delete (outputFiles2[m_filterOrder.size()]);
	cerr << "File written to: "
			<< m_prefix + "_" + NO_MATCH + "_2." + outputType + m_postfix
			<< endl;
	outputFiles1[m_filterOrder.size() + 1]->close();
	delete (outputFiles1[m_filterOrder.size() + 1]);
	cerr << "File written to: "
			<< m_prefix + "_" + MULTI_MATCH + "_1." + outputType + m_postfix
			<< endl;
	outputFiles2[m_filterOrder.size() + 1]->close();
	delete (outputFiles2[m_filterOrder.size() + 1]);
	cerr << "File written to: "
			<< m_prefix + "_" + MULTI_MATCH + "_2." + outputType + m_postfix
			<< endl;
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
	m_infoFiles.reserve(filterFilePaths.size());
	m_filters.reserve(filterFilePaths.size());
	cerr << "Starting to Load Filters." << endl;
	//load up files
	for (vector<string>::const_iterator it = filterFilePaths.begin();
			it != filterFilePaths.end(); ++it) {
		string infoFileName = (*it).substr(0, (*it).length() - 2) + "txt";

		BloomFilterInfo *temp = new BloomFilterInfo(infoFileName);

		m_infoFiles.push_back(temp);
		m_filters.push_back(new BloomFilter(*it));
		m_filterOrder.push_back(temp->getFilterID());
		cerr << "Loaded Filter: " + temp->getFilterID();
		if(opt::scoringMethod == opt::BINOMIAL){
			cerr << " FPR: " << m_filters.back()->getFPR();
		}
		cerr << endl;
	}
	cerr << "Filter Loading Complete." << endl;
}

///*
// * Collaborative filtering method
// * Assume filters use the same k-mer size
// */
//void BioBloomClassifier::evaluateReadCollab(const string &rec,
//		vector<unsigned> &hits) {
//	//get filterIDs to iterate through has in a consistent order
//	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();
//
//	//create storage for hits per filter
//	std::multimap<unsigned, string> firstPassHits;
//
//	//base for each filter until one filter obtains hit threshold
//	//TODO: staggered pattering
//	for (vector<string>::reverse_iterator i = m_filterOrder.rbegin();
//			i != m_filterOrder.rend(); ++i) {
//		hits[*i] = false;
//		unsigned screeningHits = 0;
//		size_t screeningLoc = rec.length() % kmerSize / 2;
//		//First pass filtering
//		while (rec.length() >= screeningLoc + kmerSize) {
//			const unsigned char* currentKmer = proc.prepSeq(rec, screeningLoc);
//			if (currentKmer != NULL) {
//				if (m_filtersSingle.at(*i)->contains(currentKmer)) {
//					++screeningHits;
//				}
//			}
//			screeningLoc += kmerSize;
//		}
//		firstPassHits.insert(pair<unsigned, string>(screeningHits, *i));
//	}
//
//	//evaluate promising group first
//	for (multimap<unsigned, string>::reverse_iterator i =
//			firstPassHits.rbegin(); i != firstPassHits.rend(); ++i) {
//		string filterID = i->second;
//		BloomFilter &tempFilter = *m_filtersSingle.at(filterID);
//		if (SeqEval::evalRead1(rec, kmerSize, tempFilter, m_scoreThreshold,
//				1.0 - m_scoreThreshold)) {
//			hits[filterID] = true;
//			break;
//		}
//	}
//}

///*
// * Collaborative filtering method
// * Assume filters use the same k-mer size
// */
//void BioBloomClassifier::evaluateReadCollabPair(const string &rec1,
//		const string &rec2, vector<unsigned> &hits1, vector<unsigned> &hits2) {
//	//get filterIDs to iterate through has in a consistent order
//	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();
//
//	//todo: read proc possibly unneeded, see evalSingle
//	ReadsProcessor proc(kmerSize);
//
//	//create storage for hits per filter
//	std::multimap<unsigned, string> firstPassHits;
//
//	//base for each filter until one filter obtains hit threshold
//	//TODO: staggered pattering
//	for (vector<string>::reverse_iterator i = m_filterOrder.rbegin();
//			i != m_filterOrder.rend(); ++i) {
//		hits1[*i] = false;
//		hits2[*i] = false;
//		unsigned screeningHits = 0;
//		size_t screeningLoc = rec1.length() % kmerSize / 2;
//		//First pass filtering
//		while (rec1.length() >= screeningLoc + kmerSize) {
//			const unsigned char* currentKmer1 = proc.prepSeq(rec1,
//					screeningLoc);
//			const unsigned char* currentKmer2 = proc.prepSeq(rec2,
//					screeningLoc);
//			if (currentKmer1 != NULL) {
//				if (m_filtersSingle.at(*i)->contains(currentKmer1)) {
//					++screeningHits;
//				}
//			}
//			if (currentKmer2 != NULL) {
//				if (m_filtersSingle.at(*i)->contains(currentKmer2)) {
//					++screeningHits;
//				}
//			}
//			screeningLoc += kmerSize;
//		}
//		firstPassHits.insert(pair<unsigned, string>(screeningHits, *i));
//	}
//
//	//evaluate promising group first
//	for (multimap<unsigned, string>::reverse_iterator i =
//			firstPassHits.rbegin(); i != firstPassHits.rend(); ++i) {
//		BloomFilter &tempFilter = *m_filtersSingle.at(i->second);
//		if (m_inclusive) {
//			if (SeqEval::evalRead1(rec1, kmerSize, tempFilter, m_scoreThreshold,
//					1.0 - m_scoreThreshold)
//					|| SeqEval::evalRead1(rec2, kmerSize, tempFilter,
//							m_scoreThreshold, 1.0 - m_scoreThreshold)) {
//				hits1[i->second] = true;
//				hits2[i->second] = true;
//				break;
//			}
//		} else {
//			if (SeqEval::evalRead1(rec1, kmerSize, tempFilter, m_scoreThreshold,
//					1.0 - m_scoreThreshold)
//					&& SeqEval::evalRead1(rec2, kmerSize, tempFilter,
//							m_scoreThreshold, 1.0 - m_scoreThreshold)) {
//				hits1[i->second] = true;
//				hits2[i->second] = true;
//				break;
//			}
//		}
//	}
//}

///*
// * For a single read evaluate hits for a single hash signature
// * Sections with ambiguity bases are treated as misses
// * Updates hits value to number of hits (hashSig is used to as key)
// * Faster variant that assume there a redundant tile of 0
// */
//void BioBloomClassifier::evaluateReadMin(const string &rec,
//		vector<unsigned> &hits) {
//	//get filterIDs to iterate through has in a consistent order
//	const vector<string> &idsInFilter = (*m_filters[hashSig]).getFilterIds();
//
//	//get kmersize for set of info files
//	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();
//
//	unordered_map<string, unsigned> tempHits;
//
//	//Establish tiling pattern
//	unsigned startModifier1 = (rec.length() % kmerSize) / 2;
//	size_t currentKmerNum = 0;
//
//	for (vector<string>::const_iterator i = idsInFilter.begin();
//			i != idsInFilter.end(); ++i) {
//		tempHits[*i] = 0;
//	}
//
//	ReadsProcessor proc(kmerSize);
//	//cut read into kmer size given
//	while (rec.length() >= (currentKmerNum + 1) * kmerSize) {
//
//		const unsigned char* currentKmer = proc.prepSeq(rec,
//				currentKmerNum * kmerSize + startModifier1);
//
//		//check to see if string is invalid
//		if (currentKmer != NULL) {
//
//			const vector<unsigned> &results = m_filters[hashSig]->multiContains(
//					currentKmer);
//
//			//record hit number in order
//			for (vector<string>::const_iterator i = idsInFilter.begin();
//					i != idsInFilter.end(); ++i) {
//				if (results.at(*i)) {
//					++tempHits[*i];
//				}
//			}
//		}
//		++currentKmerNum;
//	}
//	for (vector<string>::const_iterator i = idsInFilter.begin();
//			i != idsInFilter.end(); ++i) {
//		hits[*i] = tempHits.at(*i) >= m_minHit;
//	}
//}

/*
 * Ordered filtering method
 */
void BioBloomClassifier::evaluateReadOrdered(const string &rec,
		vector<unsigned> &hits) {
	for (unsigned i = 0; i != m_filters.size(); ++i) {
		if (SeqEval::evalRead(rec, *m_filters[i], m_scoreThreshold)) {
			hits.push_back(i);
			break;
		}
	}
}

/*
 * Collaborative filtering method
 * Assume filters use the same k-mer size
 */
void BioBloomClassifier::evaluateReadOrderedPair(const string &rec1,
		const string &rec2, vector<unsigned> &hits1, vector<unsigned> &hits2) {
	for (unsigned i = 0; i != m_filters.size(); ++i) {
		if (m_inclusive) {
			if (SeqEval::evalRead(rec1, *m_filters[i], m_scoreThreshold)
					|| SeqEval::evalRead(rec2, *m_filters[i],
							m_scoreThreshold)) {
				hits1.push_back(i);
				hits2.push_back(i);
				break;
			}
		} else {
			if (SeqEval::evalRead(rec1, *m_filters[i], m_scoreThreshold)
					&& SeqEval::evalRead(rec2, *m_filters[i],
							m_scoreThreshold)) {
				hits1.push_back(i);
				hits2.push_back(i);
				break;
			}
		}
	}
}

void BioBloomClassifier::evaluateReadStd(const string &rec,
		vector<unsigned> &hits) {
	for (unsigned i = 0; i != m_filters.size(); ++i) {
		if (SeqEval::evalRead(rec, *m_filters[i], m_scoreThreshold)) {
			hits.push_back(i);
		}
	}
}

/*
 * Reads are assigned to best hit
 */
double BioBloomClassifier::evaluateReadBestHit(const string &rec,
		vector<unsigned> &hits, vector<double> &scores) {
	assert(scores.size() == 0);
	vector<unsigned> bestFilters;
	double maxScore = 0;

	for (unsigned i = 0; i < m_filters.size(); ++i) {
		double score = log10(SeqEval::evalScore(rec, *m_filters[i]))*-10;
		if (maxScore < score) {
			maxScore = score;
			bestFilters.clear();
			bestFilters.push_back(i);
		} else if (maxScore == score) {
			bestFilters.push_back(i);
		}
	}
	if (maxScore > 0) {
		for (vector<unsigned>::iterator i = bestFilters.begin();
				i != bestFilters.end(); ++i) {
			hits.push_back(*i);
			scores.push_back(maxScore);
		}
	}
	return maxScore;
}

/*
 * Will return scores in vector
 * TODO: OPTIMIZE WITH HASH VALUE STORING/starting where I left off
 */
void BioBloomClassifier::evaluateReadScore(const string &rec,
		vector<unsigned> &hits, vector<double> &scores) {
	for (unsigned i = 0; i < m_filters.size(); ++i) {
		bool hit = SeqEval::evalRead(rec, *m_filters[i], m_scoreThreshold);
		if (hit) {
			hits.emplace_back(i);
		}
	}
	//compute score for multimatches
	for (unsigned i = 0; i < m_filters.size(); ++i) {
		double score = SeqEval::evalScore(rec, *m_filters[i], 0);
		scores.emplace_back(log10(score)*-10.0);
	}
}

BioBloomClassifier::~BioBloomClassifier() {
}

