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
#include <deque>
#include <cassert>

/*
 * Constructor:
 * User must specify kmer size used in sliding window and give it a list of
 * filenames with corresponding headers to make filter from.
 */
BloomFilterGenerator::BloomFilterGenerator(vector<string> const &filenames,
		uint16_t kmer) :
		kmerSize(kmer), expectedEntries(0), filterSize(0)
{

	//for each file loop over all headers and obtain max number of elements
	for (vector<string>::const_iterator i = filenames.begin();
			i != filenames.end(); ++i)
	{
		WindowedFileParser parser(*i, kmerSize);
		fileNamesAndHeaders[*i] = parser.getHeaders();
		for (vector<string>::iterator j = fileNamesAndHeaders[*i].begin();
				j != fileNamesAndHeaders[*i].end(); ++j)
		{
			//subtract kmer size for max number of possible kmers
			expectedEntries += parser.getSequenceSize(*j) - kmerSize;
		}
	}
}

/*
 * Returns a Bloom Filter generated from a file.
 * Currently only supports fasta files.
 *
 * Outputs to fileName path
 */
size_t BloomFilterGenerator::generate(string fileName)
{

	//need the number of hash functions used to be greater than 0
	assert(!hashFunctionNames.empty());

	//need the filter to be greater than the size of the number of expected entries
	assert(filterSize > expectedEntries);

	//setup bloom filter
	BloomFilter filter(filterSize, multiHash);

	//redundancy metric value
	size_t redundancy = 0;

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	for (boost::unordered_map<string, vector<string> >::iterator i =
			fileNamesAndHeaders.begin(); i != fileNamesAndHeaders.end(); ++i)
	{
		//let user know that files are being read
		cerr << "Reading File: " << i->first << endl;
		WindowedFileParser parser(i->first, kmerSize);
		for (vector<string>::iterator j = i->second.begin();
				j != i->second.end(); ++j)
		{
			parser.setLocationByHeader(*j);
			//object to process reads
			//insert elements into filter
			//read fasta file line by line and split using sliding window
			while (parser.notEndOfSeqeunce()) {
				const string &currentSeq = parser.getNextSeq();
				if (!currentSeq.empty()) {
					const vector<size_t> &tempHash = multiHash.multiHash(currentSeq);
					if (filter.contains(tempHash)) {
						redundancy++;
					} else {
						filter.insert(tempHash);
					}
				}
			}
		}
	}
	filter.storeFilter(fileName);
	return redundancy;
}

//setters
void BloomFilterGenerator::setFilterSize(size_t bits)
{
	filterSize = bits;
}

/*
 * Adds numFunc number of seed values to be used in filter
 * Returns the seed values used
 */
const vector<size_t> BloomFilterGenerator::addHashFuncs(uint16_t numFunc)
{
	vector<size_t> seedVals;
	for (uint16_t i = 0; i < numFunc; ++i) {
		multiHash.addHashFunction("CityHash64", i);
		hashFunctionNames.push_back("CityHash64");
		seedVals.push_back(i);
	}
	return seedVals;
}

//getters

/*
 * Returns the maximum possible number of expected filter entries based on inputs
 */
const size_t BloomFilterGenerator::getExpectedEntries() const
{
	return expectedEntries;
}

const vector<string> BloomFilterGenerator::getHashFuncNames() const
{
	return hashFunctionNames;
}

// todo: Tweak calculations as they are approximations and may not be 100% optimal
// see http://en.wikipedia.org/wiki/Bloom_filter

/*
 * Calculation will return optimal number of hash functions
 * to achieve lowest FPR for a given ratio of filter size and entries
 */
//Note: the floor is take because in practice you want to calculate as few hash values as possible
//NOTE: Not currently used.
const uint16_t BloomFilterGenerator::calcOptiHashNum(size_t size,
		size_t entries) const
{
	return uint16_t(((double) size / (double) entries) * log(2));
}

/*
 * Calculation assumes optimal ratio of bytes per entry given a fpr
 */
//Note: Rounded down because in practice you want to calculate as few hash values as possible
const uint16_t BloomFilterGenerator::calcOptiHashNum(float fpr) const
{
	return uint16_t(-log(fpr) / log(2));
}

//destructor
BloomFilterGenerator::~BloomFilterGenerator()
{
}
