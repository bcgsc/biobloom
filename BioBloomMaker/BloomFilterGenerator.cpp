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
#include <deque>
#include <cassert>

/*
 * Constructor:
 * User must specify kmer size used in sliding window and give it a list of
 * filenames with corresponding headers to make filter from.
 */
BloomFilterGenerator::BloomFilterGenerator(vector<string> const &filenames,
		uint16_t kmer) :
		kmerSize(kmer), expectedEntries(0), filterSize(0), totalEntries(0)
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
			expectedEntries += parser.getSequenceSize(*j) - (kmerSize - 1);
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
		uint16_t kmer, size_t numElements) :
		kmerSize(kmer), expectedEntries(numElements), filterSize(0), totalEntries(0)
{
	//for each file loop over all headers and obtain max number of elements
	for (vector<string>::const_iterator i = filenames.begin();
			i != filenames.end(); ++i)
	{
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
size_t BloomFilterGenerator::generate(const string &filename)
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
		cerr << "Processing File: " << i->first << endl;
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
					const vector<size_t> &tempHash = multiHash.multiHash(
							currentSeq);
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
 * Returns the redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files.
 *
 * Outputs to fileName path
 */
size_t BloomFilterGenerator::generate(const string &filename,
		const string &subtractFilter)
{

	//need the number of hash functions used to be greater than 0
	assert(!hashFunctionNames.empty());

	//need the filter to be greater than the size of the number of expected entries
	assert(filterSize > expectedEntries);

	//setup bloom filter
	BloomFilter filter(filterSize, multiHash);

	//load other bloom filter info
	string infoFileName = (subtractFilter).substr(0,
			(subtractFilter).length() - 2) + "txt";
	BloomFilterInfo subInfo(infoFileName);

	//createHashManager for subtraction filter
	HashManager subMan;

	vector<size_t>::const_iterator seedsItr = subInfo.getSeedValues().begin();

	for (vector<string>::const_iterator hashFnItr =
			subInfo.getHashFunctionNames().begin();
			hashFnItr != subInfo.getHashFunctionNames().end(); ++hashFnItr)
	{
		subMan.addHashFunction(*hashFnItr, *seedsItr);
		++seedsItr;
	}

	//load other bloomfilter
	BloomFilter filterSub(subInfo.getCalcuatedFilterSize(), subtractFilter,
			subMan);

	if (subInfo.getKmerSize() > kmerSize) {
		cerr
				<< "Error: Subtraction filter's kmer size is larger than output filter's kmer size."
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
			fileNamesAndHeaders.begin(); i != fileNamesAndHeaders.end(); ++i)
	{
		//let user know that files are being read
		cerr << "Processing File: " << i->first << endl;
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
					//allow kmer into filter?
					bool allowKmer = false;

					//Check if kmer or subkmers are located in filter
					if (subInfo.getKmerSize() == kmerSize) {
						//if kmer does not exists set allowance to true
						allowKmer = !filterSub.contains(currentSeq);
					} else {
						uint16_t subSections = kmerSize - subInfo.getKmerSize();

						for (uint16_t i = 0; i <= subSections; ++i) {
							if (!filterSub.contains(
									subProc.prepSeq(currentSeq, i)))
							{
								//if any sub kmer does not exists set allowance to true
								allowKmer = true;
								break;
							}
						}
					}

					if (allowKmer) {
						const vector<size_t> &tempHash = multiHash.multiHash(
								currentSeq);
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
 * Returns the total number of inserted filter entries
 */
const size_t BloomFilterGenerator::getTotalEntries() const
{
	return totalEntries;
}

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
