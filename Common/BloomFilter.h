/*
 *
 * BloomFilter.h
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_
#include <string>
#include <vector>
#include <stdint.h>
#include <omp.h>
#include "MurmurHash3.h"
#include <math.h>

using namespace std;

static const uint8_t bitsPerChar = 0x08;
static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20,
		0x40, 0x80 };

//TODO Work out better way to deal with kmerSize since conversion to kmerSize in bytes is needed

/*
 * For precomputing hash values. kmerSize is the number of bytes of the string used.
 */
static inline vector<size_t> multiHash(const unsigned char* kmer, size_t num,
		uint16_t kmerSize) {
	vector<size_t> tempHashValues(num);
	//use raw kmer number as first hash value
	size_t kmerSizeInBytes = (kmerSize + 4 - 1) / 4;

	for (size_t i = 0; i < num; ++i) {
		uint64_t hashVals[2];
		MurmurHash3_x64_128(kmer, kmerSizeInBytes, i, hashVals);
		tempHashValues[i] = hashVals[0];
		if(++i < num)
			tempHashValues[i] = hashVals[1];
	}
	return tempHashValues;
}

class BloomFilter {
public:
	//for generating a new filter
	explicit BloomFilter(size_t filterSize, uint8_t hashNum, uint16_t kmerSize);
	void insert(vector<size_t> const &precomputed);
	void insert(const unsigned char* kmer);
	const bool contains(vector<size_t> const &precomputed) const;
	const bool contains(const unsigned char* kmer) const;

	uint8_t getHashNum() const;
	uint8_t getKmerSize() const;

	//for storing/restoring the filter
	void storeFilter(string const &filterFilePath) const;
	explicit BloomFilter(size_t filterSize, uint8_t hashNum, uint16_t kmerSize,
			string const &filterFilePath);

	virtual ~BloomFilter();
private:
	void initSize(size_t size);
	unsigned char* filter;
	size_t size;
	size_t sizeInBytes;
	uint8_t hashNum;
	uint16_t kmerSize;
	uint16_t kmerSizeInBytes;
};

#endif /* BLOOMFILTER_H_ */
