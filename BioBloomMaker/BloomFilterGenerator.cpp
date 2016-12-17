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
		unsigned kmerSize, unsigned hashNum):
		m_fileNames(filenames), m_kmerSize(kmerSize), m_hashNum(hashNum), m_expectedEntries(0), m_filterSize(0), m_totalEntries(
				0){

	for (unsigned i = 0; i < m_fileNames.size(); ++i) {
		gzFile fp;
		fp = gzopen(m_fileNames[i].c_str(), "r");
		kseq_t *seq = kseq_init(fp);
		int l;
		size_t length;
#pragma omp parallel private(l, length)
		for (;;) {
#pragma omp critical(kseq_read)
			{
				l = kseq_read(seq);
				length = seq->seq.l;
			}
			if (l >= 0) {
#pragma omp atomic
				m_expectedEntries += length - m_kmerSize + 1;
			} else {
				break;
			}
		}
		kseq_destroy(seq);
		gzclose(fp);
	}
}

/*
 * Constructor:
 * User must specify kmer size used in sliding window and give it a list of
 * filenames with corresponding headers to make filter from.
 * Variant allows users to set a specific filter size before hand
 */
BloomFilterGenerator::BloomFilterGenerator(vector<string> const &filenames,
		unsigned kmerSize, unsigned hashNum, size_t numElements) :
		m_fileNames(filenames), m_kmerSize(kmerSize), m_hashNum(hashNum),  m_expectedEntries(numElements), m_filterSize(
				0), m_totalEntries(0) {
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

	if(opt::fastIO){
		redundancy = loadFilterFast(filter);
	}
	else{
		redundancy = loadFilterLowMem(filter);
	}

	filter.storeFilter(filename);
	return redundancy;
}

void BloomFilterGenerator::printReadPair(const FastqRecord& rec1,
		const FastqRecord& rec2)
{
#pragma omp critical(cout)
	{
		cout << rec1;
		cout << rec2;
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
		const string &subtractFilter )
{

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
				<< "Error: Subtraction filter's different from current filter's k-mer size."
				<< endl;
		exit(1);
	}

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	size_t redundancy = 0;

	if(opt::fastIO){
		redundancy += loadFilterFast(filter);
	}
	else{
		redundancy += loadFilterLowMem(filter);
	}

	size_t totalReads = 0;

	vector< boost::shared_ptr<ReadsProcessor> > procs(opt::threads);
	//each thread gets its own thread processor
	for(unsigned i = 0; i < opt::threads; ++i){
		procs[i] = boost::shared_ptr<ReadsProcessor>(new ReadsProcessor(m_kmerSize));
	}

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
	for (FastqRecord rec1;;) {
		FastqRecord rec2;
		bool good1;
		bool good2;

		if (m_totalEntries >= m_expectedEntries) {
			//so threshold message only printed once
			if (sequence1.eof() || sequence2.eof()) {
				break;
			}
#pragma omp critical(breakClose)
			{
				sequence1.breakClose();
				sequence2.breakClose();
				cerr << "K-mer threshold reached at read " << totalReads << endl;
			}
		}

#pragma omp critical(sequence)
		{
			good1 = sequence1 >> rec1;
			good2 = sequence2 >> rec2;
		}

		if (good1 && good2) {
#pragma omp critical(totalReads)
			{
				if (totalReads++ % 10000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}
			}
			string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
			string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
			if (tempStr1 == tempStr2) {
				size_t numKmers1 = rec1.seq.length() > m_kmerSize ? rec1.seq.length() - m_kmerSize + 1 : 0;
				size_t numKmers2 = rec2.seq.length() > m_kmerSize ? rec2.seq.length() - m_kmerSize + 1 : 0;
				vector< vector<size_t> > hashValues1(numKmers1);
				vector< vector<size_t> > hashValues2(numKmers2);
				switch (mode) {
				case PROG_INC: {
					if (numKmers1 > score && SeqEval::evalRead(rec1, m_kmerSize, filter,
							score, 1.0 - score, m_hashNum,
							hashValues1, filterSub, evalMode)) {
						//load remaining sequences
						for (unsigned i = 0; i < numKmers1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
										rec1.seq, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < numKmers2; ++i) {
							const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
									rec2.seq, i);
							insertKmer(currentSeq, filter);
						}
						if (printReads)
							printReadPair(rec1, rec2);
					} else if (numKmers2 > score && SeqEval::evalRead(rec2, m_kmerSize, filter,
							score, 1.0 - score, m_hashNum,
							hashValues2, filterSub, evalMode)) {
						//load remaining sequences
						for (unsigned i = 0; i < numKmers1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
										rec1.seq, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < numKmers2; ++i) {
							if (hashValues2[i].empty()) {
								const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
										rec2.seq, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues2[i], filter);
							}
						}
						if (printReads)
							printReadPair(rec1, rec2);
					}
					break;
				}
				case PROG_STD: {
					if (SeqEval::evalRead(rec1, m_kmerSize, filter,
							score, 1.0 - score, m_hashNum,
							hashValues1, filterSub, evalMode)
							&& SeqEval::evalRead(rec2, m_kmerSize, filter,
								score , 1.0 - score, m_hashNum, hashValues2,
								filterSub, evalMode)) {
						//load remaining sequences
						for (unsigned i = 0; i < numKmers1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
										rec1.seq, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < numKmers2; ++i) {
							if (hashValues2[i].empty()) {
								const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
										rec2.seq, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues2[i], filter);
							}
						}
						if (printReads)
							printReadPair(rec1, rec2);
					}
					break;
				}
				}
			} else {
				cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
						<< tempStr2 << endl;
				exit(1);
			}
		} else
			break;
	}
	if (!sequence1.eof() || !sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
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
		const SeqEval::EvalMode evalMode, bool printReads)
{

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	size_t redundancy = 0;

	if(opt::fastIO){
		redundancy += loadFilterFast(filter);
	}
	else{
		redundancy += loadFilterLowMem(filter);
	}

	size_t totalReads = 0;

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);

	vector< boost::shared_ptr<ReadsProcessor> > procs(opt::threads);
	//each thread gets its own thread processor
	for(unsigned i = 0; i < opt::threads; ++i){
		procs[i] = boost::shared_ptr<ReadsProcessor>(new ReadsProcessor(m_kmerSize));
	}

#pragma omp parallel
	for (FastqRecord rec1;;) {
		FastqRecord rec2;
		bool good1;
		bool good2;

		if (m_totalEntries >= m_expectedEntries) {
			//so threshold message only printed once
			if (sequence1.eof() || sequence2.eof()) {
				break;
			}
#pragma omp critical(breakClose)
			{
				sequence1.breakClose();
				sequence2.breakClose();
				cerr << "K-mer threshold reached at read " << totalReads << endl;
			}
		}

#pragma omp critical(sequence)
		{
			good1 = sequence1 >> rec1;
			good2 = sequence2 >> rec2;
		}

		if (good1 && good2) {
#pragma omp critical(totalReads)
			{
				++totalReads;
				if (totalReads % 10000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}
			}
			string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
			string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
			if (tempStr1 == tempStr2) {
				size_t numKmers1 = rec1.seq.length() > m_kmerSize ? rec1.seq.length() - m_kmerSize + 1 : 0;
				size_t numKmers2 = rec2.seq.length() > m_kmerSize ? rec2.seq.length() - m_kmerSize + 1 : 0;
//				if(numKmers2 == 0){
//					cerr << rec1.seq.length() << endl;;
//					cerr << rec2.seq.length() << endl;;
//					exit(1);
//				}
				vector<vector<size_t> > hashValues1(numKmers1);
				vector<vector<size_t> > hashValues2(numKmers2);
				switch (mode) {
				case PROG_INC: {
					if (numKmers1 > score && SeqEval::evalRead(rec1, m_kmerSize, filter,
							score, 1.0 - score, m_hashNum,
							hashValues1, evalMode)) {
						//load remaining sequences
						for (unsigned i = 0; i < numKmers1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
										rec1.seq, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < numKmers2; ++i) {
							const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
									rec2.seq, i);
							insertKmer(currentSeq, filter);
						}
						if (printReads)
							printReadPair(rec1, rec2);
					} else if (numKmers2 > score && SeqEval::evalRead(rec2, m_kmerSize, filter,
							score, 1.0 - score, m_hashNum, hashValues2, evalMode)) {
						//load remaining sequences
						for (unsigned i = 0; i < numKmers1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
										rec1.seq, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < numKmers2; ++i) {
							if (hashValues2[i].empty()) {
								const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
										rec2.seq, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues2[i], filter);
							}
						}
						if (printReads)
							printReadPair(rec1, rec2);
					}
					break;
				}
				case PROG_STD: {
					try{
					if (numKmers1 > score && SeqEval::evalRead(rec1, m_kmerSize, filter,
							score, 1.0 - score, m_hashNum, hashValues1, evalMode)
							&& numKmers2 > score && SeqEval::evalRead(rec2, m_kmerSize, filter,
								score, 1.0 - score, m_hashNum, hashValues2, evalMode)) {
						//load remaining sequences
						for (unsigned i = 0; i < numKmers1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
										rec1.seq, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < numKmers2; ++i) {
							if (hashValues2[i].empty()) {
								const unsigned char* currentSeq = procs[omp_get_thread_num()]->prepSeq(
										rec2.seq, i);
								insertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues2[i], filter);
							}
						}
						if (printReads)
							printReadPair(rec1, rec2);
					}
					} catch (const std::bad_alloc&) {
						cerr << "2 " << totalReads << endl;
						exit(6);
					}
					break;
				}
				}
			} else {
				cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
						<< tempStr2 << endl;
				exit(1);
			}
		} else
			break;
	}
	if (!sequence1.eof() || !sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
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
