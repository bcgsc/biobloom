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
	explicit BloomFilterGenerator(vector<string> const &filenames,
			uint16_t kmerSize, uint8_t hashNum, size_t numElements);

	explicit BloomFilterGenerator(vector<string> const &filenames,
			uint16_t kmerSize, uint8_t hashNum);

	size_t generate(const string &filename);
	size_t generate(const string &filename, const string &subtractFilter);
	void setFilterSize(size_t bits);

	void setHashFuncs(uint8_t numFunc);
	const size_t getTotalEntries() const;
	const size_t getExpectedEntries() const;

	virtual ~BloomFilterGenerator();
private:
	size_t filterSize;
	size_t expectedEntries;
	size_t totalEntries;
	int16_t kmerSize;
	int8_t hashNum;
	boost::unordered_map<string, vector<string> > fileNamesAndHeaders;
};

#endif /* BLOOMFILTERGENERATOR_H_ */
