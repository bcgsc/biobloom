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
			uint8_t kmer, uint8_t numFunc);
	size_t generate(string filename);
	void setFilterSize(size_t bits);
	const size_t getExpectedEntries() const;

	virtual ~BloomFilterGenerator();
private:
	size_t filterSize;
	size_t expectedEntries;
	int16_t kmerSize;
	boost::unordered_map<string, vector<string> > fileNamesAndHeaders;
};

#endif /* BLOOMFILTERGENERATOR_H_ */
