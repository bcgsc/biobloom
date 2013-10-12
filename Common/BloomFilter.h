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
#include "city.h"

using namespace std;

static const uint8_t bitsPerChar = 0x08;
static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20,
		0x40, 0x80 };

/*
 * For precomputing hash values. kmerSize is the number of bytes of the string used.
 */
static inline vector<size_t> multiHash(const char* kmer, size_t num,
		size_t kmerSize) {
	vector<size_t> tempHashValues(num);
	//use raw kmer number as first hash value
	tempHashValues[0] = kmer[0] | (kmer[1] << 8) | (kmer[2] << 16)
			| (kmer[3] << 24);
//	omp_set_num_threads(num);
//	int i;
//	#pragma omp parallel for shared(tempHashValues) private(i) schedule(static,1)
	for (int i = 1; i < num; ++i) {
		tempHashValues[i] = CityHash64WithSeed(kmer, kmerSize, i);
	}
	return tempHashValues;
}

class BloomFilter {
public:
	//for generating a new filter
	explicit BloomFilter(size_t filterSize, size_t hashNum, size_t kmerSize);
	void insert(vector<size_t> const &precomputed);
	void insert(const char* kmer);
	const bool contains(vector<size_t> const &precomputed);
	const bool contains(const char* kmer);

	uint8_t getHashNum() const;
	uint8_t getKmerSize() const;

	//for storing/restoring the filter
	void storeFilter(string const &filterFilePath) const;
	explicit BloomFilter(size_t filterSize, size_t hashNum, size_t kmerSize,
			string const &filterFilePath);

	virtual ~BloomFilter();
private:
	void initSize(size_t size);
	unsigned char* filter;
	size_t size;
	size_t sizeInBytes;
	uint8_t hashNum;
	uint8_t kmerSize;
};

#endif /* BLOOMFILTER_H_ */
