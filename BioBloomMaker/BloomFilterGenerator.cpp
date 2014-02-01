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
#include "Common/BloomFilter.h"
#include "Common/BloomFilterInfo.h"
#include <cassert>
#include <cmath>

/*
 * Constructor:
 * User must specify kmer size used in sliding window and give it a list of
 * filenames with corresponding headers to make filter from.
 */
BloomFilterGenerator::BloomFilterGenerator(vector<string> const &filenames,
		uint16_t kmerSize, uint8_t hashNum):
		kmerSize(kmerSize), hashNum(hashNum), expectedEntries(0), filterSize(0), totalEntries(
				0) {

	//for each file loop over all headers and obtain max number of elements
	for (vector<string>::const_iterator i = filenames.begin();
			i != filenames.end(); ++i) {
		WindowedFileParser parser(*i, kmerSize);
		fileNamesAndHeaders[*i] = parser.getHeaders();
		for (vector<string>::iterator j = fileNamesAndHeaders[*i].begin();
				j != fileNamesAndHeaders[*i].end(); ++j) {
			//subtract kmer size for max number of possible kmers
			expectedEntries += parser.getSequenceSize(*j) - kmerSize;
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
		uint16_t kmerSize, uint8_t hashNum, size_t numElements) :
		kmerSize(kmerSize), hashNum(hashNum),  expectedEntries(numElements), filterSize(
				0), totalEntries(0) {
	//for each file loop over all headers and obtain max number of elements
	for (vector<string>::const_iterator i = filenames.begin();
			i != filenames.end(); ++i) {
		WindowedFileParser parser(*i, kmerSize);
		fileNamesAndHeaders[*i] = parser.getHeaders();
	}
}

/*
 * Generates a bloom filter outputting it to a filename
 * Returns the redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files.
 *
 * Outputs to fileName path
 */
size_t BloomFilterGenerator::generate(const string &filename) {

	//need the number of hash functions used to be greater than 0
	assert(hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(filterSize > expectedEntries);

	//setup bloom filter
	BloomFilter filter(filterSize, hashNum, kmerSize);

	//redundancy metric value
	size_t redundancy = 0;

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	for (boost::unordered_map<string, vector<string> >::iterator i =
			fileNamesAndHeaders.begin(); i != fileNamesAndHeaders.end(); ++i) {
		//let user know that files are being read
		cerr << "Processing File: " << i->first << endl;
		WindowedFileParser parser(i->first, kmerSize);
		for (vector<string>::iterator j = i->second.begin();
				j != i->second.end(); ++j) {
			parser.setLocationByHeader(*j);
			//object to process reads
			//insert elements into filter
			//read fasta file line by line and split using sliding window
			while (parser.notEndOfSeqeunce()) {
				const unsigned char* currentSeq = parser.getNextSeq();
				if (currentSeq != NULL) {
					const vector<size_t> &tempHash = multiHash(currentSeq,
							hashNum, kmerSize);
					if (filter.contains(tempHash)) {
						redundancy++;
					} else {
						filter.insert(tempHash);
						totalEntries++;
					}
				}
			}
		}
	}
	filter.storeFilter(filename);
	return redundancy;
}

/*
 * Generates a bloom filter outputting it to a filename
 * Input a filename to use as a subtractive filter
 * Returns the redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files.
 *
 * Outputs to fileName path
 */
//TODO refactor to remove boilerplate-ness to method above
size_t BloomFilterGenerator::generate(const string &filename,
		const string &subtractFilter) {

	//need the number of hash functions used to be greater than 0
	assert(hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(filterSize > expectedEntries);

	//setup bloom filter
	BloomFilter filter(filterSize, hashNum, kmerSize);

	//load other bloom filter info
	string infoFileName = (subtractFilter).substr(0,
			(subtractFilter).length() - 2) + "txt";
	BloomFilterInfo subInfo(infoFileName);

	//load other bloomfilter
	BloomFilter filterSub(subInfo.getCalcuatedFilterSize(),
			subInfo.getHashNum(), subInfo.getKmerSize(), subtractFilter);

	if (subInfo.getKmerSize() > kmerSize) {
		cerr
				<< "Error: Subtraction filter's k-mer size is larger than output filter's k-mer size."
				<< endl;
		exit(1);
	}

	//ReadProcessor for subtraction filter
	ReadsProcessor subProc(subInfo.getKmerSize());

	//redundancy metric value
	size_t redundancy = 0;

	size_t kmerRemoved = 0;

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	for (boost::unordered_map<string, vector<string> >::iterator i =
			fileNamesAndHeaders.begin(); i != fileNamesAndHeaders.end(); ++i) {
		//let user know that files are being read
		cerr << "Processing File: " << i->first << endl;
		WindowedFileParser parser(i->first, kmerSize);
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
					if (subInfo.getKmerSize() == kmerSize) {
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
								hashNum, kmerSize);
						if (filter.contains(tempHash)) {
							redundancy++;
						} else {
							filter.insert(tempHash);
							totalEntries++;
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
	return redundancy;
}

//setters
void BloomFilterGenerator::setFilterSize(size_t bits) {
	filterSize = bits;
}

//getters

/*
 * Returns the total number of inserted filter entries
 */
const size_t BloomFilterGenerator::getTotalEntries() const {
	return totalEntries;
}

/*
 * Returns the maximum possible number of expected filter entries based on inputs
 */
const size_t BloomFilterGenerator::getExpectedEntries() const {
	return expectedEntries;
}

//destructor
BloomFilterGenerator::~BloomFilterGenerator() {
}
