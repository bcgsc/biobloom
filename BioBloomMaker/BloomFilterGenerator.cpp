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

//Todo Refactor to remove repetitive and potentially error prone parts of code

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

	if (opt::fastIO) {
		redundancy = loadFilterFast(filter);
	} else {
		redundancy = loadFilterLowMem(filter);
	}

	filter.storeFilter(filename);
	return redundancy;
}

void BloomFilterGenerator::printReadPair(const string &rec1,
		const string &header1, const string &rec2, const string &header2) {
#pragma omp critical(cout)
	{
		cout << ">" << header1 << "\n" << rec1 << "\n";
		cout << ">" << header2 << "\n" << rec2 << "\n";
	}
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
		const SeqEval::EvalMode evalMode, bool printReads,
		const string &subtractFilter) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	size_t baitFilterElements = calcExpectedEntries();

	//secondary bait filter
	BloomFilter baitFilter(
			BloomFilterInfo::calcOptimalSize(baitFilterElements, 0.0001,
					m_hashNum), m_hashNum, m_kmerSize);

	//load other bloom filter info
	string infoFileName = (subtractFilter).substr(0,
			(subtractFilter).length() - 2) + "txt";
	BloomFilterInfo subInfo(infoFileName);

	//load other bloomfilter
	BloomFilter filterSub(subInfo.getCalcuatedFilterSize(),
			subInfo.getHashNum(), subInfo.getKmerSize(), subtractFilter);

	if (subInfo.getKmerSize() != m_kmerSize) {
		cerr
				<< "Error: Subtraction filter's different from current filter's k-mer size."
				<< endl;
		exit(1);
	}

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	size_t redundancy = 0;

	if (opt::fastIO) {
		redundancy += loadFilterFast(filter);
	} else {
		redundancy += loadFilterFast(filter);
	}

	vector<boost::shared_ptr<ReadsProcessor> > procs(opt::threads);
	//each thread gets its own thread processor
	for (unsigned i = 0; i < opt::threads; ++i) {
		procs[i] = boost::shared_ptr<ReadsProcessor>(
				new ReadsProcessor(m_kmerSize));
	}

	size_t totalReads = 0;
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
		string rec1;
		string header1;
		int l2;
		string rec2;
		string header2;
#pragma omp parallel private(l1, l2, rec1, header1, rec2, header2)
		for (;;) {
#pragma omp critical(kseq_read)
			{
				l1 = kseq_read(seq1);
				if (l1 >= 0) {
					rec1 = string(seq1->seq.s, l1);
					header1 = string(seq1->name.s, seq1->name.l);
				}
				l2 = kseq_read(seq2);
				if (l2 >= 0) {
					rec2 = string(seq2->seq.s, l2);
					header2 = string(seq2->name.s, seq2->name.l);
				}
				if (l1 >= 0 && l2 >= 0) {
					++totalReads;
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
			}

			if (l1 >= 0 && l2 >= 0 && m_totalEntries < m_expectedEntries) {
				size_t numKmers1 =
						rec1.length() > m_kmerSize ? l1 - m_kmerSize + 1 : 0;
				size_t numKmers2 =
						rec2.length() > m_kmerSize ? l2 - m_kmerSize + 1 : 0;
				vector<vector<size_t> > hashValues1(numKmers1);
				vector<vector<size_t> > hashValues2(numKmers2);
				switch (mode) {
				case PROG_INC: {
					if (numKmers1 > score
							&& (SeqEval::evalRead(rec1, m_kmerSize, filter,
									score, 1.0 - score, m_hashNum, hashValues1,
									filterSub, evalMode)
									|| SeqEval::evalRead(rec1, m_kmerSize,
											baitFilter, opt::baitThreshold,
											1.0 - opt::baitThreshold, m_hashNum,
											hashValues1, filterSub, evalMode))) {
						//load remaining sequences
						for (unsigned i = 0; i < numKmers1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq =
										procs[omp_get_thread_num()]->prepSeq(
												rec1, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < numKmers2; ++i) {
							const unsigned char* currentSeq =
									procs[omp_get_thread_num()]->prepSeq(rec2,
											i);
							insertKmer(currentSeq, filter);
						}
						if (printReads)
							printReadPair(header1, rec1, header2, rec2);
					} else if (numKmers2 > score
							&& (SeqEval::evalRead(rec2, m_kmerSize, filter,
									score, 1.0 - score, m_hashNum, hashValues2,
									filterSub, evalMode)
									|| SeqEval::evalRead(rec2, m_kmerSize,
											baitFilter, opt::baitThreshold,
											1.0 - opt::baitThreshold, m_hashNum,
											hashValues2, filterSub, evalMode))) {
						//load remaining sequences
						for (unsigned i = 0; i < numKmers1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq =
										procs[omp_get_thread_num()]->prepSeq(
												rec1, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < numKmers2; ++i) {
							if (hashValues2[i].empty()) {
								const unsigned char* currentSeq =
										procs[omp_get_thread_num()]->prepSeq(
												rec2, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues2[i], filter);
							}
						}
						if (printReads)
							printReadPair(rec1, header1, rec2, header2);
					}
					break;
				}
				case PROG_STD: {
					if (SeqEval::evalRead(rec1, m_kmerSize, filter, score,
							1.0 - score, m_hashNum, hashValues1, filterSub,
							evalMode)
							&& SeqEval::evalRead(rec2, m_kmerSize, filter,
									score, 1.0 - score, m_hashNum, hashValues2,
									filterSub, evalMode)) {
						//load remaining sequences
						for (unsigned i = 0; i < numKmers1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq =
										procs[omp_get_thread_num()]->prepSeq(
												rec1, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < numKmers2; ++i) {
							if (hashValues2[i].empty()) {
								const unsigned char* currentSeq =
										procs[omp_get_thread_num()]->prepSeq(
												rec2, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues2[i], filter);
							}
						}
						if (printReads)
							printReadPair(header1, rec1, header2, rec2);
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
	return redundancy;
}

//TODO REFACTOR to minimize redundant code!
size_t BloomFilterGenerator::generateProgressive(const string &filename,
		double score, const vector<string> &files1,
		const vector<string> &files2, createMode mode,
		const SeqEval::EvalMode evalMode, bool printReads,
		const string &subtractFilter) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	size_t baitFilterElements = calcExpectedEntries();

	//secondary bait filter
	BloomFilter baitFilter(
			BloomFilterInfo::calcOptimalSize(baitFilterElements, 0.0001,
					m_hashNum), m_hashNum, m_kmerSize);

	//load other bloom filter info
	string infoFileName = (subtractFilter).substr(0,
			(subtractFilter).length() - 2) + "txt";
	BloomFilterInfo subInfo(infoFileName);

	//load other bloomfilter
	BloomFilter filterSub(subInfo.getCalcuatedFilterSize(),
			subInfo.getHashNum(), subInfo.getKmerSize(), subtractFilter);

	if (subInfo.getKmerSize() != m_kmerSize) {
		cerr
				<< "Error: Subtraction filter's different from current filter's k-mer size."
				<< endl;
		exit(1);
	}

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	size_t redundancy = 0;

	if (opt::fastIO) {
		redundancy += loadFilterFast(filter);
	} else {
		redundancy += loadFilterFast(filter);
	}

	vector<boost::shared_ptr<ReadsProcessor> > procs(opt::threads);
	//each thread gets its own thread processor
	for (unsigned i = 0; i < opt::threads; ++i) {
		procs[i] = boost::shared_ptr<ReadsProcessor>(
				new ReadsProcessor(m_kmerSize));
	}
	unsigned iter = 0;

	size_t totalReads = 0;
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
					l1 = kseq_read(seq1);
					l2 = kseq_read(seq2);
					if (l1 >= 0 && l2 >= 0
							&& m_totalEntries < m_expectedEntries) {
#pragma omp atomic
						++totalReads;
#pragma omp critical(totalReads)
						if (totalReads % 10000000 == 0) {
							cerr << "Currently Reading Read Number: "
									<< totalReads << endl;
						}
						size_t numKmers1 =
								seq1->seq.l > m_kmerSize ?
										l1 - m_kmerSize + 1 : 0;
						size_t numKmers2 =
								seq2->seq.l > m_kmerSize ?
										l2 - m_kmerSize + 1 : 0;
						vector<vector<size_t> > hashValues1(numKmers1);
						vector<vector<size_t> > hashValues2(numKmers2);
						switch (mode) {
						case PROG_INC: {
							if (numKmers1 > score
									&& (SeqEval::evalRead(seq1->seq.s,
											m_kmerSize, filter, score,
											1.0 - score, m_hashNum, hashValues1,
											filterSub, evalMode)
											|| SeqEval::evalRead(seq1->seq.s,
													m_kmerSize, baitFilter,
													opt::baitThreshold,
													1.0 - opt::baitThreshold,
													m_hashNum, hashValues1,
													filterSub, evalMode))) {
								//load remaining sequences
								for (unsigned i = 0; i < numKmers1; ++i) {
									if (hashValues1[i].empty()) {
										const unsigned char* currentSeq =
												procs[omp_get_thread_num()]->prepSeq(
														seq1->seq.s, i);
										insertKmer(currentSeq, filter);
									} else {
										insertKmer(hashValues1[i], filter);
									}
								}
								//load store second read
								for (unsigned i = 0; i < numKmers2; ++i) {
									const unsigned char* currentSeq =
											procs[omp_get_thread_num()]->prepSeq(
													seq2->seq.s, i);
									insertKmer(currentSeq, filter);
								}
								if (printReads)
									printReadPair(seq1->name.s, seq1->seq.s,
											seq2->name.s, seq2->seq.s);
							} else if (numKmers2 > score
									&& (SeqEval::evalRead(seq2->seq.s,
											m_kmerSize, filter, score,
											1.0 - score, m_hashNum, hashValues2,
											filterSub, evalMode)
											|| SeqEval::evalRead(seq2->seq.s,
													m_kmerSize, baitFilter,
													opt::baitThreshold,
													1.0 - opt::baitThreshold,
													m_hashNum, hashValues2,
													filterSub, evalMode))) {
								//load remaining sequences
								for (unsigned i = 0; i < numKmers1; ++i) {
									if (hashValues1[i].empty()) {
										const unsigned char* currentSeq =
												procs[omp_get_thread_num()]->prepSeq(
														seq1->seq.s, i);
										insertKmer(currentSeq, filter);
									} else {
										insertKmer(hashValues1[i], filter);
									}
								}
								//load store second read
								for (unsigned i = 0; i < numKmers2; ++i) {
									if (hashValues2[i].empty()) {
										const unsigned char* currentSeq =
												procs[omp_get_thread_num()]->prepSeq(
														seq2->seq.s, i);
										insertKmer(currentSeq, filter);
									} else {
										insertKmer(hashValues2[i], filter);
									}
								}
								if (printReads)
									printReadPair(seq1->name.s, seq1->seq.s,
											seq2->name.s, seq2->seq.s);
							}
							break;
						}
						case PROG_STD: {
							if (SeqEval::evalRead(seq1->seq.s, m_kmerSize,
									filter, score, 1.0 - score, m_hashNum,
									hashValues1, filterSub, evalMode)
									&& SeqEval::evalRead(seq2->seq.s,
											m_kmerSize, filter, score,
											1.0 - score, m_hashNum, hashValues2,
											filterSub, evalMode)) {
								//load remaining sequences
								for (unsigned i = 0; i < numKmers1; ++i) {
									if (hashValues1[i].empty()) {
										const unsigned char* currentSeq =
												procs[omp_get_thread_num()]->prepSeq(
														seq1->seq.s, i);
										insertKmer(currentSeq, filter);
									} else {
										insertKmer(hashValues1[i], filter);
									}
								}
								//load store second read
								for (unsigned i = 0; i < numKmers2; ++i) {
									if (hashValues2[i].empty()) {
										const unsigned char* currentSeq =
												procs[omp_get_thread_num()]->prepSeq(
														seq2->seq.s, i);
										insertKmer(currentSeq, filter);
									} else {
										insertKmer(hashValues2[i], filter);
									}
								}
								if (printReads)
									printReadPair(seq1->name.s, seq1->seq.s,
											seq2->name.s, seq2->seq.s);
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
		const SeqEval::EvalMode evalMode, bool printReads) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	size_t redundancy = 0;

	if (opt::fastIO) {
		redundancy += loadFilterFast(filter);
	} else {
		redundancy += loadFilterLowMem(filter);
	}

	size_t baitFilterElements = calcExpectedEntries();

	//secondary bait filter
	BloomFilter baitFilter(
			BloomFilterInfo::calcOptimalSize(baitFilterElements, 0.0001,
					m_hashNum), m_hashNum, m_kmerSize);

	vector<boost::shared_ptr<ReadsProcessor> > procs(opt::threads);
	//each thread gets its own thread processor
	for (unsigned i = 0; i < opt::threads; ++i) {
		procs[i] = boost::shared_ptr<ReadsProcessor>(
				new ReadsProcessor(m_kmerSize));
	}

	size_t totalReads = 0;
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
		string rec1;
		string header1;
		int l2;
		string rec2;
		string header2;

		if (m_totalEntries < m_expectedEntries) {
#pragma omp parallel private(l1, l2, rec1, header1, rec2, header2)
			for (;;) {
#pragma omp critical(kseq_read)
				{
					l1 = kseq_read(seq1);
					if (l1 >= 0) {
						rec1 = string(seq1->seq.s, l1);
						header1 = string(seq1->name.s, seq1->name.l);
					}
					l2 = kseq_read(seq2);
					if (l2 >= 0) {
						rec2 = string(seq2->seq.s, l2);
						header2 = string(seq2->name.s, seq2->name.l);
					}
					if (l1 >= 0 && l2 >= 0) {
						++totalReads;
						if (totalReads % 10000000 == 0) {
							cerr << "Currently Reading Read Number: "
									<< totalReads << endl;
						}
					}
				}

				if (l1 >= 0 && l2 >= 0 && m_totalEntries < m_expectedEntries) {
					size_t numKmers1 =
							rec1.length() > m_kmerSize ?
									rec1.length() - m_kmerSize + 1 : 0;
					size_t numKmers2 =
							rec2.length() > m_kmerSize ?
									rec2.length() - m_kmerSize + 1 : 0;
					vector<vector<size_t> > hashValues1(numKmers1);
					vector<vector<size_t> > hashValues2(numKmers2);
					switch (mode) {
					case PROG_INC: {
						if (numKmers1 > score
								&& (SeqEval::evalRead(rec1, m_kmerSize, filter,
										score, 1.0 - score, m_hashNum,
										hashValues1, evalMode)
										|| SeqEval::evalRead(rec1, m_kmerSize,
												baitFilter, opt::baitThreshold,
												1.0 - opt::baitThreshold,
												m_hashNum, hashValues1,
												evalMode))) {
							//load remaining sequences
							for (unsigned i = 0; i < numKmers1; ++i) {
								if (hashValues1[i].empty()) {
									const unsigned char* currentSeq =
											procs[omp_get_thread_num()]->prepSeq(
													rec1, i);
									insertKmer(currentSeq, filter);
								} else {
									insertKmer(hashValues1[i], filter);
								}
							}
							//load store second read
							for (unsigned i = 0; i < numKmers2; ++i) {
								const unsigned char* currentSeq =
										procs[omp_get_thread_num()]->prepSeq(
												rec2, i);
								insertKmer(currentSeq, filter);
							}
							if (printReads)
								printReadPair(rec1, header1, rec2, header2);
						} else if (numKmers2 > score
								&& (SeqEval::evalRead(rec1, m_kmerSize, filter,
										score, 1.0 - score, m_hashNum,
										hashValues1, evalMode)
										|| SeqEval::evalRead(rec1, m_kmerSize,
												baitFilter, opt::baitThreshold,
												1.0 - opt::baitThreshold,
												m_hashNum, hashValues1,
												evalMode))) {
							//load remaining sequences
							for (unsigned i = 0; i < numKmers1; ++i) {
								if (hashValues1[i].empty()) {
									const unsigned char* currentSeq =
											procs[omp_get_thread_num()]->prepSeq(
													rec1, i);
									insertKmer(currentSeq, filter);
								} else {
									insertKmer(hashValues1[i], filter);
								}
							}
							//load store second read
							for (unsigned i = 0; i < numKmers2; ++i) {
								if (hashValues2[i].empty()) {
									const unsigned char* currentSeq =
											procs[omp_get_thread_num()]->prepSeq(
													rec2, i);
									insertKmer(currentSeq, filter);
								} else {
									insertKmer(hashValues2[i], filter);
								}
							}
							if (printReads)
								printReadPair(rec1, header1, rec2, header2);
						}
						break;
					}
					case PROG_STD: {
						if (numKmers1 > score
								&& (SeqEval::evalRead(rec1, m_kmerSize, filter,
										score, 1.0 - score, m_hashNum,
										hashValues1, evalMode)
										|| SeqEval::evalRead(rec1, m_kmerSize,
												baitFilter, opt::baitThreshold,
												1.0 - opt::baitThreshold,
												m_hashNum, hashValues1,
												evalMode)) && numKmers2 > score
								&& (SeqEval::evalRead(rec1, m_kmerSize, filter,
										score, 1.0 - score, m_hashNum,
										hashValues1, evalMode)
										|| SeqEval::evalRead(rec1, m_kmerSize,
												baitFilter, opt::baitThreshold,
												1.0 - opt::baitThreshold,
												m_hashNum, hashValues1,
												evalMode))) {
							//load remaining sequences
							for (unsigned i = 0; i < numKmers1; ++i) {
								if (hashValues1[i].empty()) {
									const unsigned char* currentSeq =
											procs[omp_get_thread_num()]->prepSeq(
													rec1, i);
									insertKmer(currentSeq, filter);
								} else {
									insertKmer(hashValues1[i], filter);
								}
							}
							//load store second read
							for (unsigned i = 0; i < numKmers2; ++i) {
								if (hashValues2[i].empty()) {
									const unsigned char* currentSeq =
											procs[omp_get_thread_num()]->prepSeq(
													rec2, i);
									insertKmer(currentSeq, filter);
								} else {
									insertKmer(hashValues2[i], filter);
								}
							}
							if (printReads)
								printReadPair(rec1, header1, rec2, header2);
						}
						break;
					}
					}
				} else
					break;
			}
		}
		if (m_totalEntries >= m_expectedEntries) {
			cerr << "K-mer threshold reached at read " << totalReads << endl;
		} else {
			cerr << "K-mer threshold not reached, number of k-mers: "
					<< m_totalEntries << endl;
		}
	}

	filter.storeFilter(filename);
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
	BloomFilter filterSub(subInfo.getCalcuatedFilterSize(),
			subInfo.getHashNum(), subInfo.getKmerSize(), subtractFilter);

	if (subInfo.getKmerSize() != m_kmerSize) {
		cerr
				<< "Error: Subtraction filter's k-mer size is a different size from than output filter's k-mer size."
				<< endl;
		exit(1);
	}

	size_t redundancy = loadFilterFastSubtract(filter, filterSub);

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
