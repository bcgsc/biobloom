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

//todo: Performance optimization potential: make size of bitset a power of 2
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
//todo: OptPot: add piped decompression/compression support for loading speed
BloomFilter::BloomFilter(size_t filterSize, string const &filterFilePath,
		HashManager const &hashFns) :
		size(filterSize), multiHasher(hashFns)
{
	initSize(size);

	//Check file size is correct size
	struct stat sb;
	stat(filterFilePath.c_str(), &sb);
	if (sb.st_size != filterSize / 8) {
		cerr << "Error: " << filterFilePath
				<< " does not match size given by its information file. Size: "
				<< sb.st_size << "/" << filterSize / 8 << " bytes." << endl;
		exit(1);
	}

	//load in blocks to a vector
	ifstream binaryFile(filterFilePath.c_str(), ios::binary);
	if (binaryFile.is_open()) {
		binaryFile.read(filter, size);
	} else {
		cerr << "file \"" << filterFilePath << "\" could not be read." << endl;
		binaryFile.close();
		exit(1);
	}
	binaryFile.close();
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

void BloomFilter::insert(string const &kmer)
{
	vector<size_t> values = multiHasher.multiHash(kmer);
	//iterates through hashed values adding it to the filter
	for (vector<size_t>::iterator it = values.begin(); it != values.end(); ++it)
	{
		size_t normalizedValue = *it % size;
		filter[normalizedValue / bitsPerChar] |= bitMask[normalizedValue
				% bitsPerChar];
	}
}

/*
 * Stores the filter as a binary file to the path specified
 */
void BloomFilter::storeFilter(string const &filterFilePath) const
{
	ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

	//todo: OptPot: send data to a compression stream to decrease bytes needed to be written

	//write out each block
	for (int i = 0; i < sizeInBytes; i++) {
		myFile << filter[i];
	}
	myFile.close();
}

const bool BloomFilter::contains(string const &kmer)
{
	return contains(multiHasher.multiHash(kmer));
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

/*
 * need to cleanup dynamically allocated filter
 */
BloomFilter::~BloomFilter()
{
	delete (filter);
}
