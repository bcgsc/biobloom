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
#include "Common/BloomFilter.h"
using namespace std;

enum createMode{PROG_STD, PROG_INC};

class BloomFilterGenerator {
public:
	explicit BloomFilterGenerator(vector<string> const &filenames,
			unsigned kmerSize, unsigned hashNum, size_t numElements);

	explicit BloomFilterGenerator(vector<string> const &filenames,
			unsigned kmerSize, unsigned hashNum);


	//TODO: THREAD ME!
	size_t generate(const string &filename);
	size_t generate(const string &filename, const string &subtractFilter);
	size_t generateProgressive(const string &filename, double score, const string &file1,
			const string &file2, createMode mode);
	void setFilterSize(size_t bits);

	void setHashFuncs(unsigned numFunc);
	size_t getTotalEntries() const;
	size_t getExpectedEntries() const;

	virtual ~BloomFilterGenerator();
private:
	unsigned m_kmerSize;
	unsigned m_hashNum;
	size_t m_expectedEntries;
	size_t m_filterSize;
	size_t m_totalEntries;
	size_t m_redundancy;

	boost::unordered_map<string, vector<string> > m_fileNamesAndHeaders;

	inline void checkAndInsertKmer(const unsigned char* currentSeq,
			BloomFilter &filter)
	{
		if (currentSeq != NULL) {
			const vector<size_t> &tempHash = multiHash(currentSeq, m_hashNum,
					m_kmerSize);
			insertKmer(tempHash, filter);
		}
	}

	inline void insertKmer(const vector<size_t> &hashVals,
			BloomFilter &filter)
	{
		if (filter.contains(hashVals)) {
#pragma omp atomic
			m_redundancy++;
		} else {
			filter.insert(hashVals);
#pragma omp atomic
			m_totalEntries++;
		}
	}
};

#endif /* BLOOMFILTERGENERATOR_H_ */
