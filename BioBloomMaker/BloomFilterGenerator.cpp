/*
 * BloomFilterGenerator.cpp
 *
 * Using an input file (currently fasta), generates a bloom filter
 * The filter is returned. Intended to create a filter information
 * object as well (not implemented).
 *
 *  Created on: Jun 20, 2012
 *      Author: cjustin
 */
#include "BloomFilterGenerator.h"
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include "Common/BloomFilterInfo.h"
#include <cassert>
#include <cmath>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

/*
 * Constructor:
 * User must specify kmer size used in sliding window and give it a list of
 * filenames with corresponding headers to make filter from.
 */
BloomFilterGenerator::BloomFilterGenerator(vector<string> const &filenames,
		unsigned kmerSize, unsigned hashNum) :
		m_fileNames(filenames), m_kmerSize(kmerSize), m_hashNum(hashNum), m_expectedEntries(
				0), m_filterSize(0), m_totalEntries(0) {
	m_expectedEntries = calcExpectedEntries();
}

/*
 * Constructor:
 * User must specify kmer size used in sliding window and give it a list of
 * filenames with corresponding headers to make filter from.
 * Variant allows users to set a specific filter size before hand
 */
BloomFilterGenerator::BloomFilterGenerator(vector<string> const &filenames,
		unsigned kmerSize, unsigned hashNum, size_t numElements) :
		m_fileNames(filenames), m_kmerSize(kmerSize), m_hashNum(hashNum), m_expectedEntries(
				numElements), m_filterSize(0), m_totalEntries(0) {
}

/*
 * Generates a bloom filter outputting it to a filename
 * Returns the m_redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files.
 *
 * Outputs to fileName path
 */
size_t BloomFilterGenerator::generate(const string &filename) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	size_t redundancy = 0;
	redundancy += loadFilter(filter, m_totalEntries);
	cerr
			<< "Approximated (due to false positives) total unique k-mers in reference files "
			<< m_totalEntries << endl;

	filter.storeFilter(filename);
	return redundancy;
}

/*
 * Generates a bloom filter outputting it to a filename
 * Returns the m_redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files for initial seeding
 *
 * Uses fastq files to stream in additional sequence
 * Stops after m_expectedEntries of entries
 *
 * Outputs to fileName path
 */
size_t BloomFilterGenerator::generateProgressive(const string &filename,
		double score, const string &file1, const string &file2, createMode mode,
		bool printReads, const string &subtractFilter) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	BloomFilter* filterSub = NULL;

	if (subtractFilter != "") {
		filterSub = new BloomFilter(subtractFilter);
		checkFilters(filter, *filterSub);
	}
	//for each file loop over all headers and obtain seq
	//load input file + make filter
	size_t redundancy = 0;

	if (opt::noRep && filterSub != NULL) {
		redundancy += loadFilterSubtract(filter, *filterSub, m_totalEntries);
	} else {
		redundancy += loadFilter(filter, m_totalEntries);
	}
	cerr
			<< "Approximated (due to false positives) total unique k-mers in reference files "
			<< m_totalEntries << endl;

	size_t totalReads = 0;
	size_t taggedReads = 0;
	for (unsigned i = 0; i < opt::progItrns; ++i) {
		cerr << "Iteration " << i + 1 << endl;

		gzFile fp1;
		gzFile fp2;

		fp1 = gzopen(file1.c_str(), "r");
		if (fp1 == Z_NULL) {
#pragma omp critical(cerr)
			cerr << "file " << file1.c_str() << " cannot be opened" << endl;
			exit(1);
		} else {
#pragma omp critical(cerr)
			cerr << "Reading file " << file1.c_str() << endl;
		}
		fp2 = gzopen(file2.c_str(), "r");
		if (fp2 == Z_NULL) {
#pragma omp critical(cerr)
			cerr << "file " << file2.c_str() << " cannot be opened" << endl;
			exit(1);
		} else {
#pragma omp critical(cerr)
			cerr << "Reading file " << file2.c_str() << endl;
		}
		kseq_t *seq1 = kseq_init(fp1);
		kseq_t *seq2 = kseq_init(fp2);
		int l1;
		FqRec rec1;
		int l2;
		FqRec rec2;
#pragma omp parallel private(l1, l2, rec1, rec2)
		for (;;) {
#pragma omp critical(kseq_read)
			{
				l1 = kseq_read(seq1);
				if (l1 >= 0) {
					rec1.seq = string(seq1->seq.s, l1);
					rec1.header = string(seq1->name.s, seq1->name.l);
					rec1.qual = string(seq1->qual.s, seq1->qual.l);
				}
				l2 = kseq_read(seq2);
				if (l2 >= 0) {
					rec2.seq = string(seq2->seq.s, l2);
					rec2.header = string(seq2->name.s, seq2->name.l);
					rec2.qual = string(seq2->qual.s, seq2->qual.l);
				}
				if (l1 >= 0 && l2 >= 0) {
					++totalReads;
					if (totalReads % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< "\tUnique k-mers Added: " << m_totalEntries
								<< "\tReads Used in Tagging: " << taggedReads
								<< endl;
					}
				}
			}

			if (l1 >= 0 && l2 >= 0 && m_totalEntries < m_expectedEntries) {
				size_t numKmers1 =
						rec1.seq.length() > m_kmerSize ?
								l1 - m_kmerSize + 1 : 0;
				size_t numKmers2 =
						rec2.seq.length() > m_kmerSize ?
								l2 - m_kmerSize + 1 : 0;
				switch (mode) {
				case PROG_INC: {
					if (numKmers1 > score
							&& (SeqEval::evalRead(rec1.seq, filter, score,
									filterSub))) {
#pragma omp atomic
						++taggedReads;
						if (printReads) {
							unsigned taggedKmers1 = loadFilter(filter,
									rec1.seq);
							unsigned taggedKmers2 = loadFilter(filter,
									rec2.seq);
							unsigned repCount1 = checkFilter(filterSub,
									rec1.seq);
							unsigned repCount2 = checkFilter(filterSub,
									rec2.seq);
#pragma omp critical(debugPrint)
							{
								printDebug(rec1, taggedKmers1, repCount1,
										taggedReads, totalReads);
								printDebug(rec2, taggedKmers2, repCount2,
										taggedReads, totalReads);
							}
						} else {
							loadFilter(filter, rec1.seq);
							loadFilter(filter, rec2.seq);
						}
					} else if (numKmers2 > score
							&& (SeqEval::evalRead(rec2.seq, filter, score,
									filterSub))) {
#pragma omp atomic
						++taggedReads;
						if (printReads) {
							unsigned taggedKmers1 = loadFilter(filter,
									rec1.seq);
							unsigned taggedKmers2 = loadFilter(filter,
									rec2.seq);
							unsigned repCount1 = checkFilter(filterSub,
									rec1.seq);
							unsigned repCount2 = checkFilter(filterSub,
									rec2.seq);
#pragma omp critical(debugPrint)
							{
								printDebug(rec1, taggedKmers1, repCount1,
										taggedReads, totalReads);
								printDebug(rec2, taggedKmers2, repCount2,
										taggedReads, totalReads);
							}
						} else {
							loadFilter(filter, rec1.seq);
							loadFilter(filter, rec2.seq);
						}
					}
					break;
				}
				case PROG_STD: {
					if (SeqEval::evalRead(rec1.seq, filter, score, filterSub)
							&& SeqEval::evalRead(rec2.seq, filter, score,
									filterSub)) {
#pragma omp atomic
						++taggedReads;
						if (printReads) {
							unsigned taggedKmers1 = loadFilter(filter,
									rec1.seq);
							unsigned taggedKmers2 = loadFilter(filter,
									rec2.seq);
							unsigned repCount1 = checkFilter(filterSub,
									rec1.seq);
							unsigned repCount2 = checkFilter(filterSub,
									rec2.seq);
#pragma omp critical(debugPrint)
							{
								printDebug(rec1, taggedKmers1, repCount1,
										taggedReads, totalReads);
								printDebug(rec2, taggedKmers2, repCount2,
										taggedReads, totalReads);
							}
						} else {
							loadFilter(filter, rec1.seq);
							loadFilter(filter, rec2.seq);
						}
					}
					break;
				}
				}
			} else
				break;
		}
		kseq_destroy(seq1);
		kseq_destroy(seq2);
		gzclose(fp1);
		gzclose(fp2);
	}
	if (m_totalEntries >= m_expectedEntries) {
		cerr << "K-mer threshold reached at read " << totalReads << endl;
	} else {
		cerr << "K-mer threshold not reached, number of k-mers: "
				<< m_totalEntries << endl;
	}

	filter.storeFilter(filename);
	if (filterSub != NULL) {
		delete (filterSub);
	}
	return redundancy;
}

/*
 * Generates a bloom filter outputting it to a filename
 * Returns the m_redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files for initial seeding
 *
 * Uses fastq files to stream in additional sequence
 * Stops after m_expectedEntries of entries
 *
 * Outputs to fileName path
 */
size_t BloomFilterGenerator::generateProgressiveBait(const string &filename,
		double score, const string &file1, const string &file2, createMode mode,
		bool printReads, const string &subtractFilter) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	BloomFilter* filterSub = NULL;

	if (subtractFilter != "") {
		filterSub = new BloomFilter(subtractFilter);
		checkFilters(filter, *filterSub);
	}

	size_t baitFilterElements = calcExpectedEntries();

	//secondary bait filter
	BloomFilter baitFilter(
			BloomFilterInfo::calcOptimalSize(baitFilterElements, opt::fpr,
					m_hashNum), m_hashNum, m_kmerSize);

	size_t totalEntriesBait = 0;
	loadFilter(baitFilter, totalEntriesBait);

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	size_t redundancy = 0;
	if (opt::noRep && filterSub != NULL) {
		redundancy += loadFilterSubtract(filter, *filterSub, m_totalEntries);
	} else {
		redundancy += loadFilter(filter, m_totalEntries);
	}
	cerr
			<< "Approximated (due to false positives) total unique k-mers in reference files "
			<< m_totalEntries << endl;

	size_t totalReads = 0;
	size_t taggedReads = 0;
	for (unsigned i = 0; i < opt::progItrns; ++i) {
		cerr << "Iteration " << i + 1 << endl;

		gzFile fp1;
		gzFile fp2;

		fp1 = gzopen(file1.c_str(), "r");
		if (fp1 == Z_NULL) {
#pragma omp critical(cerr)
			cerr << "file " << file1.c_str() << " cannot be opened" << endl;
			exit(1);
		} else {
#pragma omp critical(cerr)
			cerr << "Reading file " << file1.c_str() << endl;
		}
		fp2 = gzopen(file2.c_str(), "r");
		if (fp2 == Z_NULL) {
#pragma omp critical(cerr)
			cerr << "file " << file2.c_str() << " cannot be opened" << endl;
			exit(1);
		} else {
#pragma omp critical(cerr)
			cerr << "Reading file " << file2.c_str() << endl;
		}
		kseq_t *seq1 = kseq_init(fp1);
		kseq_t *seq2 = kseq_init(fp2);
		int l1;
		FqRec rec1;
		int l2;
		FqRec rec2;
#pragma omp parallel private(l1, l2, rec1, rec2)
		for (;;) {
#pragma omp critical(kseq_read)
			{
				l1 = kseq_read(seq1);
				if (l1 >= 0) {
					rec1.seq = string(seq1->seq.s, l1);
					rec1.header = string(seq1->name.s, seq1->name.l);
					rec1.qual = string(seq1->qual.s, seq1->qual.l);
				}
				l2 = kseq_read(seq2);
				if (l2 >= 0) {
					rec2.seq = string(seq2->seq.s, l2);
					rec2.header = string(seq2->name.s, seq2->name.l);
					rec2.qual = string(seq2->qual.s, seq2->qual.l);
				}
				if (l1 >= 0 && l2 >= 0) {
					++totalReads;
					if (totalReads % opt::fileInterval == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< "\tUnique k-mers Added: " << m_totalEntries
								<< "\tReads Used in Tagging: " << taggedReads
								<< endl;
					}
				}
			}

			if (l1 >= 0 && l2 >= 0 && m_totalEntries < m_expectedEntries) {
				size_t numKmers1 =
						rec1.seq.length() > m_kmerSize ?
								l1 - m_kmerSize + 1 : 0;
				size_t numKmers2 =
						rec2.seq.length() > m_kmerSize ?
								l2 - m_kmerSize + 1 : 0;
				switch (mode) {
				case PROG_INC: {
					if (numKmers1 > score
							&& (SeqEval::evalRead(rec1.seq, filter, score,
									filterSub)
									|| SeqEval::evalRead(rec1.seq, baitFilter,
											opt::baitThreshold, filterSub))) {
#pragma omp atomic
						++taggedReads;
						if (printReads) {
							unsigned taggedKmers1 = loadFilter(filter,
									rec1.seq);
							unsigned taggedKmers2 = loadFilter(filter,
									rec2.seq);
							unsigned repCount1 = checkFilter(filterSub,
									rec1.seq);
							unsigned repCount2 = checkFilter(filterSub,
									rec2.seq);
#pragma omp critical(debugPrint)
							{
								printDebug(rec1, taggedKmers1, repCount1,
										taggedReads, totalReads);
								printDebug(rec2, taggedKmers2, repCount2,
										taggedReads, totalReads);
							}
						} else {
							loadFilter(filter, rec1.seq);
							loadFilter(filter, rec2.seq);
						}
					} else if (numKmers2 > score
							&& (SeqEval::evalRead(rec2.seq, filter, score,
									filterSub)
									|| SeqEval::evalRead(rec2.seq, baitFilter,
											opt::baitThreshold, filterSub))) {
#pragma omp atomic
						++taggedReads;
						if (printReads) {
							unsigned taggedKmers1 = loadFilter(filter,
									rec1.seq);
							unsigned taggedKmers2 = loadFilter(filter,
									rec2.seq);
							unsigned repCount1 = checkFilter(filterSub,
									rec1.seq);
							unsigned repCount2 = checkFilter(filterSub,
									rec2.seq);
#pragma omp critical(debugPrint)
							{
								printDebug(rec1, taggedKmers1, repCount1,
										taggedReads, totalReads);
								printDebug(rec2, taggedKmers2, repCount2,
										taggedReads, totalReads);
							}
						} else {
							loadFilter(filter, rec1.seq);
							loadFilter(filter, rec2.seq);
						}
					}
					break;
				}
				case PROG_STD: {
					if ((SeqEval::evalRead(rec1.seq, filter, score, filterSub)
							|| SeqEval::evalRead(rec1.seq, baitFilter,
									opt::baitThreshold, filterSub))
							&& (SeqEval::evalRead(rec2.seq, filter, score,
									filterSub)
									|| SeqEval::evalRead(rec2.seq, baitFilter,
											opt::baitThreshold, filterSub))) {
#pragma omp atomic
						++taggedReads;
						if (printReads) {
							unsigned taggedKmers1 = loadFilter(filter,
									rec1.seq);
							unsigned taggedKmers2 = loadFilter(filter,
									rec2.seq);
							unsigned repCount1 = checkFilter(filterSub,
									rec1.seq);
							unsigned repCount2 = checkFilter(filterSub,
									rec2.seq);
#pragma omp critical(debugPrint)
							{
								printDebug(rec1, taggedKmers1, repCount1,
										taggedReads, totalReads);
								printDebug(rec2, taggedKmers2, repCount2,
										taggedReads, totalReads);
							}
						} else {
							loadFilter(filter, rec1.seq);
							loadFilter(filter, rec2.seq);
						}
					}
					break;
				}
				}
			} else
				break;
		}
		kseq_destroy(seq1);
		kseq_destroy(seq2);
		gzclose(fp1);
		gzclose(fp2);
	}
	if (m_totalEntries >= m_expectedEntries) {
		cerr << "K-mer threshold reached at read " << totalReads << endl;
	} else {
		cerr << "K-mer threshold not reached, number of k-mers: "
				<< m_totalEntries << endl;
	}

	filter.storeFilter(filename);
	if (filterSub != NULL) {
		delete (filterSub);
	}
	return redundancy;
}

size_t BloomFilterGenerator::generateProgressiveBait(const string &filename,
		double score, const vector<string> &files1,
		const vector<string> &files2, createMode mode, bool printReads,
		const string &subtractFilter) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	BloomFilter* filterSub = NULL;

	if (subtractFilter != "") {
		filterSub = new BloomFilter(subtractFilter);
		checkFilters(filter, *filterSub);
	}

	size_t baitFilterElements = calcExpectedEntries();

	//secondary bait filter
	BloomFilter baitFilter(
			BloomFilterInfo::calcOptimalSize(baitFilterElements, opt::fpr,
					m_hashNum), m_hashNum, m_kmerSize);

	size_t totalEntriesBait = 0;
	loadFilter(baitFilter, totalEntriesBait);

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	size_t redundancy = 0;

	if (opt::noRep && filterSub != NULL) {
	} else {
		redundancy += loadFilter(filter, m_totalEntries);
	}
	cerr
			<< "Approximated (due to false positives) total unique k-mers in reference files "
			<< m_totalEntries << endl;

	unsigned iter = 0;

	size_t totalReads = 0;
	size_t taggedReads = 0;
	if (m_totalEntries < m_expectedEntries) {
		for (; iter < opt::progItrns; ++iter) {
			cerr << "Iteration " << iter + 1 << endl;

#pragma omp parallel for
			for (unsigned i = 0; i < files1.size(); ++i) {
				gzFile fp1;
				gzFile fp2;

				fp1 = gzopen(files1[i].c_str(), "r");
				if (fp1 == Z_NULL) {
#pragma omp critical(cerr)
					cerr << "file " << files1[i].c_str() << " cannot be opened"
							<< endl;
					exit(1);
				} else {
#pragma omp critical(cerr)
					cerr << "Reading File: " << files1[i].c_str() << endl;
				}
				fp2 = gzopen(files2[i].c_str(), "r");
				if (fp2 == Z_NULL) {
#pragma omp critical(cerr)
					cerr << "file " << files2[i].c_str() << " cannot be opened"
							<< endl;
					exit(1);
				} else {
#pragma omp critical(cerr)
					cerr << "Reading File: " << files2[i].c_str() << endl;
				}
				kseq_t *seq1 = kseq_init(fp1);
				kseq_t *seq2 = kseq_init(fp2);
				int l1;
				int l2;
				for (;;) {
					{
						l1 = kseq_read(seq1);
						l2 = kseq_read(seq2);
						if (l1 >= 0 && l2 >= 0) {
#pragma omp critical(err)
							if (++totalReads % opt::fileInterval == 0) {
								cerr << "Currently Reading Read Number: "
										<< totalReads
										<< "\tUnique k-mers Added: "
										<< m_totalEntries
										<< "\tReads Used in Tagging: "
										<< taggedReads << endl;
							}
						}
					}

					if (l1 >= 0 && l2 >= 0
							&& m_totalEntries < m_expectedEntries) {
						size_t numKmers1 =
								seq1->seq.l > m_kmerSize ?
										seq1->seq.l - m_kmerSize + 1 : 0;
						size_t numKmers2 =
								seq2->seq.l > m_kmerSize ?
										seq2->seq.l - m_kmerSize + 1 : 0;
						switch (mode) {
						case PROG_INC: {
							if (numKmers1 > score
									&& (SeqEval::evalRead(seq1->seq.s, filter,
											score, filterSub)
											|| SeqEval::evalRead(seq1->seq.s,
													baitFilter,
													opt::baitThreshold,
													filterSub))) {
#pragma omp atomic
								++taggedReads;
								if (printReads) {
									unsigned taggedKmers1 = loadFilter(filter,
											seq1->seq.s);
									unsigned taggedKmers2 = loadFilter(filter,
											seq2->seq.s);
									unsigned repCount1 = checkFilter(filterSub,
											seq1->seq.s);
									unsigned repCount2 = checkFilter(filterSub,
											seq2->seq.s);
#pragma omp critical(debugPrint)
									{
										FqRec rec1, rec2;
										rec1.seq = string(seq1->seq.s, l1);
										rec1.header = string(seq1->name.s,
												seq1->name.l);
										rec1.qual = string(seq1->qual.s,
												seq1->qual.l);
										rec2.seq = string(seq2->seq.s, l2);
										rec2.header = string(seq2->name.s,
												seq2->name.l);
										rec2.qual = string(seq2->qual.s,
												seq2->qual.l);
										printDebug(rec1, taggedKmers1,
												repCount1, taggedReads,
												totalReads);
										printDebug(rec2, taggedKmers2,
												repCount2, taggedReads,
												totalReads);
									}
								} else {
									loadFilter(filter, seq1->seq.s);
									loadFilter(filter, seq2->seq.s);
								}
							} else if (numKmers2 > score
									&& (SeqEval::evalRead(seq2->seq.s, filter,
											score, filterSub)
											|| SeqEval::evalRead(seq2->seq.s,
													baitFilter,
													opt::baitThreshold,
													filterSub))) {
#pragma omp atomic
								++taggedReads;
								if (printReads) {
									unsigned taggedKmers1 = loadFilter(filter,
											seq1->seq.s);
									unsigned taggedKmers2 = loadFilter(filter,
											seq2->seq.s);
									unsigned repCount1 = checkFilter(filterSub,
											seq1->seq.s);
									unsigned repCount2 = checkFilter(filterSub,
											seq2->seq.s);
#pragma omp critical(debugPrint)
									{
										FqRec rec1, rec2;
										rec1.seq = string(seq1->seq.s, l1);
										rec1.header = string(seq1->name.s,
												seq1->name.l);
										rec1.qual = string(seq1->qual.s,
												seq1->qual.l);
										rec2.seq = string(seq2->seq.s, l2);
										rec2.header = string(seq2->name.s,
												seq2->name.l);
										rec2.qual = string(seq2->qual.s,
												seq2->qual.l);
										printDebug(rec1, taggedKmers1,
												repCount1, taggedReads,
												totalReads);
										printDebug(rec2, taggedKmers2,
												repCount2, taggedReads,
												totalReads);
									}
								} else {
									loadFilter(filter, seq1->seq.s);
									loadFilter(filter, seq2->seq.s);
								}
							}
							break;
						}
						case PROG_STD: {
							if ((SeqEval::evalRead(seq1->seq.s, filter, score,
									filterSub)
									|| SeqEval::evalRead(seq1->seq.s,
											baitFilter, opt::baitThreshold,
											filterSub))
									&& (SeqEval::evalRead(seq2->seq.s, filter,
											score, filterSub)
											|| SeqEval::evalRead(seq2->seq.s,
													baitFilter,
													opt::baitThreshold,
													filterSub))) {
#pragma omp atomic
								++taggedReads;
								if (printReads) {
									unsigned taggedKmers1 = loadFilter(filter,
											seq1->seq.s);
									unsigned taggedKmers2 = loadFilter(filter,
											seq2->seq.s);
									unsigned repCount1 = checkFilter(filterSub,
											seq1->seq.s);
									unsigned repCount2 = checkFilter(filterSub,
											seq2->seq.s);
#pragma omp critical(debugPrint)
									{
										FqRec rec1, rec2;
										rec1.seq = string(seq1->seq.s, l1);
										rec1.header = string(seq1->name.s,
												seq1->name.l);
										rec1.qual = string(seq1->qual.s,
												seq1->qual.l);
										rec2.seq = string(seq2->seq.s, l2);
										rec2.header = string(seq2->name.s,
												seq2->name.l);
										rec2.qual = string(seq2->qual.s,
												seq2->qual.l);
										printDebug(rec1, taggedKmers1,
												repCount1, taggedReads,
												totalReads);
										printDebug(rec2, taggedKmers2,
												repCount2, taggedReads,
												totalReads);
									}
								} else {
									loadFilter(filter, seq1->seq.s);
									loadFilter(filter, seq2->seq.s);
								}
							}
							break;
						}
						}
					} else
						break;
				}
				kseq_destroy(seq1);
				kseq_destroy(seq2);
				gzclose(fp1);
				gzclose(fp2);
			}
		}
	}
	if (m_totalEntries >= m_expectedEntries) {
		cerr << "K-mer threshold reached at read " << totalReads << endl;
	} else {
		cerr << "K-mer threshold not reached, number of k-mers: "
				<< m_totalEntries << endl;
	}

	filter.storeFilter(filename);
	if (filterSub != NULL) {
		delete (filterSub);
	}
	return redundancy;
}

size_t BloomFilterGenerator::generateProgressive(const string &filename,
		double score, const vector<string> &files1,
		const vector<string> &files2, createMode mode, bool printReads,
		const string &subtractFilter) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	BloomFilter* filterSub = NULL;

	if (subtractFilter != "") {
		filterSub = new BloomFilter(subtractFilter);
		checkFilters(filter, *filterSub);
	}

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	size_t redundancy = 0;
	if (opt::noRep) {
		redundancy += loadFilterSubtract(filter, *filterSub, m_totalEntries);
	} else {
		redundancy += loadFilter(filter, m_totalEntries);
	}
	cerr
			<< "Approximated (due to false positives) total unique k-mers in reference files "
			<< m_totalEntries << endl;

	unsigned iter = 0;

	size_t totalReads = 0;
	size_t taggedReads = 0;
	if (m_totalEntries < m_expectedEntries) {
		for (; iter < opt::progItrns; ++iter) {
			cerr << "Iteration " << iter + 1 << endl;

#pragma omp parallel for
			for (unsigned i = 0; i < files1.size(); ++i) {
				gzFile fp1;
				gzFile fp2;

				fp1 = gzopen(files1[i].c_str(), "r");
				if (fp1 == Z_NULL) {
#pragma omp critical(cerr)
					cerr << "file " << files1[i].c_str() << " cannot be opened"
							<< endl;
					exit(1);
				} else {
#pragma omp critical(cerr)
					cerr << "Reading File: " << files1[i].c_str() << endl;
				}
				fp2 = gzopen(files2[i].c_str(), "r");
				if (fp2 == Z_NULL) {
#pragma omp critical(cerr)
					cerr << "file " << files2[i].c_str() << " cannot be opened"
							<< endl;
					exit(1);
				} else {
#pragma omp critical(cerr)
					cerr << "Reading File: " << files2[i].c_str() << endl;
				}
				kseq_t *seq1 = kseq_init(fp1);
				kseq_t *seq2 = kseq_init(fp2);
				int l1;
				int l2;
				for (;;) {
					{
						l1 = kseq_read(seq1);
						l2 = kseq_read(seq2);
						if (l1 >= 0 && l2 >= 0) {
#pragma omp critical(err)
							if (++totalReads % opt::fileInterval == 0) {
								cerr << "Currently Reading Read Number: "
										<< totalReads
										<< "\tUnique k-mers Added: "
										<< m_totalEntries
										<< "\tReads Used in Tagging: "
										<< taggedReads << endl;
							}
						}
					}

					if (l1 >= 0 && l2 >= 0
							&& m_totalEntries < m_expectedEntries) {
						size_t numKmers1 =
								seq1->seq.l > m_kmerSize ?
										seq1->seq.l - m_kmerSize + 1 : 0;
						size_t numKmers2 =
								seq2->seq.l > m_kmerSize ?
										seq2->seq.l - m_kmerSize + 1 : 0;
						switch (mode) {
						case PROG_INC: {
							if (numKmers1 > score
									&& (SeqEval::evalRead(seq1->seq.s, filter,
											score, filterSub))) {
#pragma omp atomic
								++taggedReads;
								if (printReads) {
									unsigned taggedKmers1 = loadFilter(filter,
											seq1->seq.s);
									unsigned taggedKmers2 = loadFilter(filter,
											seq2->seq.s);
									unsigned repCount1 = checkFilter(filterSub,
											seq1->seq.s);
									unsigned repCount2 = checkFilter(filterSub,
											seq2->seq.s);
#pragma omp critical(debugPrint)
									{
										FqRec rec1, rec2;
										rec1.seq = string(seq1->seq.s, l1);
										rec1.header = string(seq1->name.s,
												seq1->name.l);
										rec1.qual = string(seq1->qual.s,
												seq1->qual.l);
										rec2.seq = string(seq2->seq.s, l2);
										rec2.header = string(seq2->name.s,
												seq2->name.l);
										rec2.qual = string(seq2->qual.s,
												seq2->qual.l);
										printDebug(rec1, taggedKmers1,
												repCount1, taggedReads,
												totalReads);
										printDebug(rec2, taggedKmers2,
												repCount2, taggedReads,
												totalReads);
									}
								} else {
									loadFilter(filter, seq1->seq.s);
									loadFilter(filter, seq2->seq.s);
								}
							} else if (numKmers2 > score
									&& (SeqEval::evalRead(seq2->seq.s, filter,
											score, filterSub))) {
#pragma omp atomic
								++taggedReads;
								if (printReads) {
									unsigned taggedKmers1 = loadFilter(filter,
											seq1->seq.s);
									unsigned taggedKmers2 = loadFilter(filter,
											seq2->seq.s);
									unsigned repCount1 = checkFilter(filterSub,
											seq1->seq.s);
									unsigned repCount2 = checkFilter(filterSub,
											seq2->seq.s);
#pragma omp critical(debugPrint)
									{
										FqRec rec1, rec2;
										rec1.seq = string(seq1->seq.s, l1);
										rec1.header = string(seq1->name.s,
												seq1->name.l);
										rec1.qual = string(seq1->qual.s,
												seq1->qual.l);
										rec2.seq = string(seq2->seq.s, l2);
										rec2.header = string(seq2->name.s,
												seq2->name.l);
										rec2.qual = string(seq2->qual.s,
												seq2->qual.l);
										printDebug(rec1, taggedKmers1,
												repCount1, taggedReads,
												totalReads);
										printDebug(rec2, taggedKmers2,
												repCount2, taggedReads,
												totalReads);
									}
								} else {
									loadFilter(filter, seq1->seq.s);
									loadFilter(filter, seq2->seq.s);
								}
							}
							break;
						}
						case PROG_STD: {
							if (SeqEval::evalRead(seq1->seq.s, filter, score,
									filterSub)
									&& SeqEval::evalRead(seq2->seq.s, filter,
											score, filterSub)) {
#pragma omp atomic
								++taggedReads;
								if (printReads) {
									unsigned taggedKmers1 = loadFilter(filter,
											seq1->seq.s);
									unsigned taggedKmers2 = loadFilter(filter,
											seq2->seq.s);
									unsigned repCount1 = checkFilter(filterSub,
											seq1->seq.s);
									unsigned repCount2 = checkFilter(filterSub,
											seq2->seq.s);
#pragma omp critical(debugPrint)
									{
										FqRec rec1, rec2;
										rec1.seq = string(seq1->seq.s, l1);
										rec1.header = string(seq1->name.s,
												seq1->name.l);
										rec1.qual = string(seq1->qual.s,
												seq1->qual.l);
										rec2.seq = string(seq2->seq.s, l2);
										rec2.header = string(seq2->name.s,
												seq2->name.l);
										rec2.qual = string(seq2->qual.s,
												seq2->qual.l);
										printDebug(rec1, taggedKmers1,
												repCount1, taggedReads,
												totalReads);
										printDebug(rec2, taggedKmers2,
												repCount2, taggedReads,
												totalReads);
									}
								} else {
									loadFilter(filter, seq1->seq.s);
									loadFilter(filter, seq2->seq.s);
								}
							}
							break;
						}
						}
					} else
						break;
				}
				kseq_destroy(seq1);
				kseq_destroy(seq2);
				gzclose(fp1);
				gzclose(fp2);
			}
		}
	}
	if (m_totalEntries >= m_expectedEntries) {
		cerr << "K-mer threshold reached at read " << totalReads << endl;
	} else {
		cerr << "K-mer threshold not reached, number of k-mers: "
				<< m_totalEntries << endl;
	}

	filter.storeFilter(filename);
	if (filterSub != NULL) {
		delete (filterSub);
	}
	return redundancy;
}

size_t BloomFilterGenerator::generateProgressive(const string &filename,
		double score, const vector<string> &files, bool printReads,
		const string &subtractFilter) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	BloomFilter* filterSub = NULL;

	if (subtractFilter != "") {
		filterSub = new BloomFilter(subtractFilter);
		checkFilters(filter, *filterSub);
	}

	size_t baitFilterElements = calcExpectedEntries();

	//secondary bait filter
	BloomFilter baitFilter(
			BloomFilterInfo::calcOptimalSize(baitFilterElements, opt::fpr,
					m_hashNum), m_hashNum, m_kmerSize);

	size_t totalEntriesBait = 0;
	loadFilter(baitFilter, totalEntriesBait);

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	size_t redundancy = 0;

	if (opt::noRep && filterSub != NULL) {
	} else {
		redundancy += loadFilter(filter, m_totalEntries);
	}
	cerr
			<< "Approximated (due to false positives) total unique k-mers in reference files "
			<< m_totalEntries << endl;

	unsigned iter = 0;

	size_t totalReads = 0;
	size_t taggedReads = 0;
	if (m_totalEntries < m_expectedEntries) {
		for (; iter < opt::progItrns; ++iter) {
			cerr << "Iteration " << iter + 1 << endl;

#pragma omp parallel for
			for (unsigned i = 0; i < files.size(); ++i) {
				gzFile fp = gzopen(files[i].c_str(), "r");
				if (fp == Z_NULL) {
#pragma omp critical(cerr)
					cerr << "file " << files[i].c_str() << " cannot be opened"
							<< endl;
					exit(1);
				} else {
#pragma omp critical(cerr)
					cerr << "Reading File: " << files[i].c_str() << endl;
				}

				kseq_t *seq = kseq_init(fp);
				int l;
				for (;;) {
					l = kseq_read(seq);

					if (l >= 0 && m_totalEntries < m_expectedEntries) {
#pragma omp atomic
						++totalReads;
#pragma omp critical(err)
						if (totalReads % opt::fileInterval == 0) {
							cerr << "Currently Reading Read Number: "
									<< totalReads << "\tUnique k-mers Added: "
									<< m_totalEntries
									<< "\tReads Used in Tagging: "
									<< taggedReads << endl;
						}
						size_t numKmers =
								seq->seq.l > m_kmerSize ?
										seq->seq.l - m_kmerSize + 1 : 0;
						if (numKmers > score
								&& (SeqEval::evalRead(seq->seq.s, filter, score,
										filterSub))) {
#pragma omp atomic
							++taggedReads;
							if (printReads) {
								unsigned taggedKmers = loadFilter(filter,
										seq->seq.s);
								unsigned repCount = checkFilter(filterSub,
										seq->seq.s);
#pragma omp critical(debugPrint)
								{
									FqRec rec;
									rec.seq = string(seq->seq.s, l);
									rec.header = string(seq->name.s,
											seq->name.l);
									rec.qual = string(seq->qual.s, seq->qual.l);
									printDebug(rec, taggedKmers, repCount,
											taggedReads, totalReads);
								}
							} else {
								loadFilter(filter, seq->seq.s);
							}
						}
					} else
						break;
				}
				kseq_destroy(seq);
				gzclose(fp);
			}
		}
	}
	if (m_totalEntries >= m_expectedEntries) {
		cerr << "K-mer threshold reached at read " << totalReads << endl;
	} else {
		cerr << "K-mer threshold not reached, number of k-mers: "
				<< m_totalEntries << endl;
	}

	filter.storeFilter(filename);
	if (filterSub != NULL) {
		delete (filterSub);
	}
	return redundancy;
}

/*
 * Generates a bloom filter outputting it to a filename
 * Input a filename to use as a subtractive filter
 * Returns the m_redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files.
 *
 * Outputs to fileName path
 */
size_t BloomFilterGenerator::generate(const string &filename,
		const string &subtractFilter) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	//load other bloom filter info
	string infoFileName = (subtractFilter).substr(0,
			(subtractFilter).length() - 2) + "txt";
	BloomFilterInfo subInfo(infoFileName);

	//load other bloomfilter
	BloomFilter filterSub(subtractFilter);

	size_t redundancy = loadFilterSubtract(filter, filterSub, m_totalEntries);

	filter.storeFilter(filename);
	return redundancy;
}

//setters
void BloomFilterGenerator::setFilterSize(size_t bits) {
	m_filterSize = bits;
}

//getters

/*
 * Returns the total number of inserted filter entries
 */
size_t BloomFilterGenerator::getTotalEntries() const {
	return m_totalEntries;
}

/*
 * Returns the maximum possible number of expected filter entries based on inputs
 */
size_t BloomFilterGenerator::getExpectedEntries() const {
	return m_expectedEntries;
}

//destructor
BloomFilterGenerator::~BloomFilterGenerator() {
}
