/*
 * BloomFilter.h
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_
#include "Common/HashManager.h"
#include <string>
#include <boost/dynamic_bitset.hpp>

using namespace std;

class BloomFilter {
public:
	//for generating a new filter
	BloomFilter(size_t const &filterSize, HashManager const &hashFns);
	void insert(string const &kmer);

	//for storing/restoring the filter
	void storeFilter(string const &filterFilePath) const;
	BloomFilter(string const &filterFilePath, HashManager const &hashFns);

	const bool contains(string const &kmer);
	const bool contains(vector<size_t> const &precomputed);

	//bool contains(vector<size_t>); //If values have been pre-hashed

	virtual ~BloomFilter();
private:
	boost::dynamic_bitset<> filter;
	HashManager multiHasher;
};

#endif /* BLOOMFILTER_H_ */
