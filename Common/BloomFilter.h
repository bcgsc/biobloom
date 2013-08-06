/*
 * BloomFilter.h
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_
#include <string>
#include <vector>

using namespace std;

static const size_t bitsPerChar = 0x08;
static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20,
		0x40, 0x80 };

class BloomFilter {
public:
	//for generating a new filter
	explicit BloomFilter(size_t filterSize, size_t hashNum);
	void insert(vector<size_t> const &precomputed);

	//for storing/restoring the filter
	void storeFilter(string const &filterFilePath) const;
	BloomFilter(size_t filterSize, size_t hashNum, string const &filterFilePath);

	const bool contains(vector<size_t> const &precomputed);
	virtual ~BloomFilter();
private:
	void initSize(size_t size);
	char* filter;
	size_t size;
	size_t sizeInBytes;
	size_t hashNum;
};

#endif /* BLOOMFILTER_H_ */
