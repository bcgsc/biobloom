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
	explicit BloomFilterGenerator(
			vector<string> const &filenames,
			int16_t kmer);
	void generate(string filename);
	void setFilterSize(size_t bits);
	const vector<size_t> addHashFuncs(int16_t numFunc);
	const size_t getExpectedEntries() const;
	const vector<string> getHashFuncNames() const;
	const int16_t calcOptiHashNum(size_t size, size_t entries) const;
	const int16_t calcOptiHashNum(float fpr) const;

	virtual ~BloomFilterGenerator();
private:
	size_t filterSize;
	size_t expectedEntries;
	int16_t kmerSize;
	HashManager multiHash;
	vector<string> hashFunctionNames;
	boost::unordered_map<string, vector<string> > fileNamesAndHeaders;
};

#endif /* BLOOMFILTERGENERATOR_H_ */
