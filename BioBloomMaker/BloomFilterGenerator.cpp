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
#include "WindowedFileParser.h"
#include "Common/BloomFilterInfo.h"
#include <cassert>
#include <cmath>
#include <FastaReader.h>
#include "Common/SeqEval.h"

/*
 * Constructor:
 * User must specify kmer size used in sliding window and give it a list of
 * filenames with corresponding headers to make filter from.
 */
BloomFilterGenerator::BloomFilterGenerator(vector<string> const &filenames,
		unsigned kmerSize, unsigned hashNum):
		m_kmerSize(kmerSize), m_hashNum(hashNum), m_expectedEntries(0), m_filterSize(0), m_totalEntries(
				0), m_redundancy(0){

	//for each file loop over all headers and obtain max number of elements
	for (vector<string>::const_iterator i = filenames.begin();
			i != filenames.end(); ++i) {
		WindowedFileParser parser(*i, kmerSize);
		m_fileNamesAndHeaders[*i] = parser.getHeaders();
		for (vector<string>::iterator j = m_fileNamesAndHeaders[*i].begin();
				j != m_fileNamesAndHeaders[*i].end(); ++j) {
			//subtract kmer size for max number of possible kmers
			m_expectedEntries += parser.getSequenceSize(*j) - kmerSize + 1;
		}
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
		m_kmerSize(kmerSize), m_hashNum(hashNum),  m_expectedEntries(numElements), m_filterSize(
				0), m_totalEntries(0), m_redundancy(0) {
	//for each file loop over all headers and obtain max number of elements
	for (vector<string>::const_iterator i = filenames.begin();
			i != filenames.end(); ++i) {
		WindowedFileParser parser(*i, kmerSize);
		m_fileNamesAndHeaders[*i] = parser.getHeaders();
	}
}

/*
 * Generates a bloom filter outputting it to a filename
 * Returns the m_redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files.
 *
 * Outputs to fileName path
 * TODO: AGAIN TRY TO MAKE THIS LESS BOILERPLATE-LIKE
 */
size_t BloomFilterGenerator::generate(const string &filename) {

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	for (boost::unordered_map<string, vector<string> >::iterator i =
			m_fileNamesAndHeaders.begin(); i != m_fileNamesAndHeaders.end(); ++i) {
		//let user know that files are being read
		cerr << "Processing File: " << i->first << endl;
		WindowedFileParser parser(i->first, m_kmerSize);
		for (vector<string>::iterator j = i->second.begin();
				j != i->second.end(); ++j) {
			parser.setLocationByHeader(*j);
			//object to process reads
			//insert elements into filter
			//read fasta file line by line and split using sliding window
			while (parser.notEndOfSeqeunce()) {
				const unsigned char* currentSeq = parser.getNextSeq();
				checkAndInsertKmer(currentSeq, filter);
			}
		}
	}
	filter.storeFilter(filename);
	return m_redundancy;
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
		double score, const string &file1, const string &file2, createMode mode)
{

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	for (boost::unordered_map<string, vector<string> >::iterator i =
			m_fileNamesAndHeaders.begin(); i != m_fileNamesAndHeaders.end(); ++i) {
		//let user know that files are being read
		cerr << "Processing File: " << i->first << endl;
		WindowedFileParser parser(i->first, m_kmerSize);
		for (vector<string>::iterator j = i->second.begin();
				j != i->second.end(); ++j) {
			parser.setLocationByHeader(*j);
			//object to process reads
			//insert elements into filter
			//read fasta file line by line and split using sliding window
			while (parser.notEndOfSeqeunce()) {
				const unsigned char* currentSeq = parser.getNextSeq();
				checkAndInsertKmer(currentSeq, filter);
			}
		}
	}

	size_t totalReads = 0;

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
	for (FastqRecord rec1;;) {
		FastqRecord rec2;
		bool good1;
		bool good2;

		if (m_totalEntries >= m_expectedEntries) {
			//so threshold message only printed once
			if (sequence1.eof() || sequence1.eof()) {
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
			ReadsProcessor proc(m_kmerSize);
			string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
			string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
			if (tempStr1 == tempStr2) {
				unsigned size1 = rec1.seq.length() - m_kmerSize + 1;
				unsigned size2 = rec2.seq.length() - m_kmerSize + 1;
				vector<vector<size_t> > hashValues1(size1);
				vector<vector<size_t> > hashValues2(size2);
				switch (mode) {
				case PROG_INC: {
					if (SeqEval::evalSingle(rec1, m_kmerSize, filter,
									score * double(size1),
									(1.0 - score) * double(size1), m_hashNum,
									hashValues1)) {
						//load remaining sequences
						for (unsigned i = 0; i < size1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq = proc.prepSeq(
										rec1.seq, i);
								checkAndInsertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < size2; ++i) {
							const unsigned char* currentSeq = proc.prepSeq(
									rec2.seq, i);
							checkAndInsertKmer(currentSeq, filter);
						}
					} else if (SeqEval::evalSingle(rec2, m_kmerSize, filter,
									score * size2, (1.0 - score) * size2, m_hashNum,
									hashValues2)) {
						//load remaining sequences
						for (unsigned i = 0; i < size1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq = proc.prepSeq(
										rec1.seq, i);
								checkAndInsertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < size2; ++i) {
							if (hashValues2[i].empty()) {
								const unsigned char* currentSeq = proc.prepSeq(
										rec2.seq, i);
								checkAndInsertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues2[i], filter);
							}
						}
					}
					break;
				}
				case PROG_STD: {
					if (SeqEval::evalSingle(rec1, m_kmerSize, filter,
							score * double(size1),
							(1.0 - score) * double(size1), m_hashNum,
							hashValues1)
							&& SeqEval::evalSingle(rec2, m_kmerSize, filter,
									score * size2, (1.0 - score) * size2,
									m_hashNum, hashValues2)) {
						//load remaining sequences
						for (unsigned i = 0; i < size1; ++i) {
							if (hashValues1[i].empty()) {
								const unsigned char* currentSeq = proc.prepSeq(
										rec1.seq, i);
								checkAndInsertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues1[i], filter);
							}
						}
						//load store second read
						for (unsigned i = 0; i < size2; ++i) {
							if (hashValues2[i].empty()) {
								const unsigned char* currentSeq = proc.prepSeq(
										rec2.seq, i);
								checkAndInsertKmer(currentSeq, filter);
							} else {
								insertKmer(hashValues2[i], filter);
							}
						}
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
	return m_redundancy;
}


/*
 * Generates a bloom filter outputting it to a filename
 * Input a filename to use as a subtractive filter
 * Returns the m_redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files.
 *
 * Outputs to fileName path
 */
//TODO refactor to remove boilerplate-ness to method above
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

	if (subInfo.getKmerSize() > m_kmerSize) {
		cerr
				<< "Error: Subtraction filter's k-mer size is larger than output filter's k-mer size."
				<< endl;
		exit(1);
	}

	//ReadProcessor for subtraction filter
	ReadsProcessor subProc(subInfo.getKmerSize());

	size_t kmerRemoved = 0;

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	for (boost::unordered_map<string, vector<string> >::iterator i =
			m_fileNamesAndHeaders.begin(); i != m_fileNamesAndHeaders.end(); ++i) {
		//let user know that files are being read
		cerr << "Processing File: " << i->first << endl;
		WindowedFileParser parser(i->first, m_kmerSize);
		for (vector<string>::iterator j = i->second.begin();
				j != i->second.end(); ++j) {
			parser.setLocationByHeader(*j);
			//object to process reads
			//insert elements into filter
			//read fasta file line by line and split using sliding window
			while (parser.notEndOfSeqeunce()) {
				const unsigned char* currentSeq = parser.getNextSeq();
				if (currentSeq != NULL) {
					//allow kmer into filter?
					bool allowKmer = false;

					//Check if kmer or subkmers are located in filter
					if (subInfo.getKmerSize() == m_kmerSize) {
						//if kmer does not exists set allowance to true
						allowKmer = !filterSub.contains(currentSeq);
					} else {
						//TODO make compatable with smaller kmer sizes
						cerr
								<< "ERROR: Must use identical size k-mers in subtractive filter"
								<< endl;
//						uint16_t subSections = kmerSize - kmerSize;
//						for (uint16_t i = 0; i <= subSections; ++i) {
//							if (!filterSub.contains(subProc.prepSeq(currentSeq, i)))
//							{
//								//if any sub kmer does not exists set allowance to true
//								allowKmer = true;
//								break;
//							}
//						}
					}

					if (allowKmer) {
						const vector<size_t> &tempHash = multiHash(currentSeq,
								m_hashNum, m_kmerSize);
						if (filter.contains(tempHash)) {
							m_redundancy++;
						} else {
							filter.insert(tempHash);
							m_totalEntries++;
						}
					} else {
						++kmerRemoved;
					}
				}
			}
		}
	}

	cerr << "Total Number of K-mers not added: " << kmerRemoved << endl;

	filter.storeFilter(filename);
	return m_redundancy;
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
