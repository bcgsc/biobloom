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
using namespace std;

class BloomFilterGenerator {
public:
	explicit BloomFilterGenerator(vector<string> const &filenames,uint8_t kmer, uint8_t numFunc, size_t numElements);
	explicit BloomFilterGenerator(vector<string> const &filenames,uint8_t kmer, uint8_t numFunc);
	size_t generate(const string &filename);
	size_t generate(const string &filename, const string &subtractFilter);
	void setFilterSize(size_t bits);
	const vector<size_t> addHashFuncs(uint16_t numFunc);
	const size_t getTotalEntries() const;
	const size_t getExpectedEntries() const;

	virtual ~BloomFilterGenerator();
private:
	size_t filterSize;
	size_t expectedEntries;
	size_t totalEntries;
	int16_t kmerSize;
	boost::unordered_map<string, vector<string> > fileNamesAndHeaders;
};

#endif /* BLOOMFILTERGENERATOR_H_ */
