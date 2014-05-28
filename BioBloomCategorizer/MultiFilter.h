/*
 * MultiFilter.h
 *
 *  Created on: Sep 7, 2012
 *      Author: cjustin
 */

#ifndef MULTIFILTER_H_
#define MULTIFILTER_H_
#include <string>
#include <vector>
#include "Common/BloomFilter.h"
#include "boost/unordered/unordered_map.hpp"
#include "boost/shared_ptr.hpp"
using namespace std;

class MultiFilter {
public:
	MultiFilter(uint16_t hashNum, uint16_t kmerSize);
	void addFilter(string const &filterID, boost::shared_ptr<BloomFilter> filter);
	const boost::unordered_map<string, bool> multiContains(const unsigned char* kmer);
	const boost::unordered_map<string, bool> multiContains(const unsigned char* kmer,
			vector<string> const &tempFilters);
	const BloomFilter &getFilter(const string &filterID);
	const vector<string> &getFilterIds() const;
	virtual ~MultiFilter();
private:
	boost::unordered_map<string, boost::shared_ptr<BloomFilter> > filters;
	uint16_t hashNum;
	uint16_t kmerSize;
	vector<string> filterIDs;
};

#endif /* MULTIFILTER_H_ */
