/*
 * BloomFilter.cpp
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#include "BloomFilter.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <stdio.h>
#include <cstring>

/* De novo filter constructor.
 * precondition: filterSize must be a multiple of 64
 */
BloomFilter::BloomFilter(size_t filterSize, HashManager const &hashFns) :
		size(filterSize), multiHasher(hashFns)
{
	initSize(size);
	memset(filter, 0, sizeInBytes);
}

/*
 * Loads the filter (file is a .bf file) from path specified
 */
BloomFilter::BloomFilter(size_t filterSize, string const &filterFilePath,
		HashManager const &hashFns) :
		size(filterSize), multiHasher(hashFns)
{
	initSize(size);

	FILE *file = fopen(filterFilePath.c_str(), "rb");
	if (file == NULL) {
		cerr << "file \"" << filterFilePath << "\" could not be read." << endl;
		exit(1);
	}

	long int lCurPos = ftell(file);
	fseek(file, 0, 2);
	size_t fileSize = ftell(file);
	fseek(file, lCurPos, 0);
	if (fileSize != sizeInBytes) {
		cerr << "Error: " << filterFilePath
				<< " does not match size given by its information file. Size: "
				<< fileSize << " vs " << sizeInBytes << " bytes." << endl;
		exit(1);
	}

	fread(filter, fileSize, 1, file);
	if (fclose(file) != 0) {
		cerr << "file \"" << filterFilePath << "\" could not be read." << endl;
		exit(1);
	}
}

/*
 * Checks filter size and initializes filter
 */
void BloomFilter::initSize(size_t size)
{
	if (size % 8 != 0) {
		cerr << "ERROR: Filter Size \"" << size << "\" is not a multiple of 8."
				<< endl;
		exit(1);
	}
	sizeInBytes = size / bitsPerChar;
	filter = new char[sizeInBytes];
}

/*
 * Accepts a list of precomputed hash values. Faster than rehashing each time.
 */
void BloomFilter::insert(vector<size_t> const &precomputed)
{
	//iterates through hashed values adding it to the filter

	for (vector<size_t>::const_iterator it = precomputed.begin();
			it != precomputed.end(); ++it)
	{
		size_t normalizedValue = *it % size;
		filter[normalizedValue / bitsPerChar] |= bitMask[normalizedValue
				% bitsPerChar];
	}
}

/*
 * Stores the filter as a binary file to the path specified
 * Stores uncompressed because the random data tend to
 * compress poorly anyway
 */
void BloomFilter::storeFilter(string const &filterFilePath) const
{
	ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

	assert(myFile.good());
	//write out each block
	myFile.write(filter, sizeInBytes);
	myFile.close();
	assert(myFile);
}

/*
 * Accepts a list of precomputed hash values. Faster than rehashing each time.
 */
const bool BloomFilter::contains(vector<size_t> const &values)
{
	for (vector<size_t>::const_iterator it = values.begin(); it != values.end();
			++it)
	{
		size_t normalizedValue = *it % size;
		char bit = bitMask[normalizedValue % bitsPerChar];
		if ((filter[normalizedValue / bitsPerChar] & bit) != bit) {
			return false;
		}
	}
	return true;
}

BloomFilter::~BloomFilter()
{
	delete[] filter;
}
