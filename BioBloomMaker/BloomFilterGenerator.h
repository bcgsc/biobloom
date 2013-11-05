/*
 * BloomFilterGenerator.h
 *
 *  Created on: Jun 20, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTERGENERATOR_H_
#define BLOOMFILTERGENERATOR_H_
#include <boost/unordered/unordered_map.hpp>
#include <vector>
#include "Common/HashManager.h"
using namespace std;

class BloomFilterGenerator {
public:
	explicit BloomFilterGenerator(vector<string> const &filenames,
			uint16_t kmer);
	explicit BloomFilterGenerator(vector<string> const &filenames,
			uint16_t kmer, size_t numElements);
	size_t generate(const string &filename);
	size_t generate(const string &filename, const string &subtractFilter);
	void setFilterSize(size_t bits);
	const vector<size_t> addHashFuncs(uint16_t numFunc);
	const size_t getTotalEntries() const;
	const size_t getExpectedEntries() const;
	const vector<string> getHashFuncNames() const;
	const uint16_t calcOptiHashNum(size_t size, size_t entries) const;
	const uint16_t calcOptiHashNum(float fpr) const;

	virtual ~BloomFilterGenerator();
private:
	size_t filterSize;
	size_t expectedEntries;
	size_t totalEntries;
	int16_t kmerSize;
	HashManager multiHash;
	vector<string> hashFunctionNames;
	boost::unordered_map<string, vector<string> > fileNamesAndHeaders;
};

#endif /* BLOOMFILTERGENERATOR_H_ */
