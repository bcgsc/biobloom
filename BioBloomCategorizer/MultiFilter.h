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
#include "Common/HashManager.h"
#include "Common/BloomFilter.h"
#include "boost/unordered/unordered_map.hpp"
#include "boost/shared_ptr.hpp"
using namespace std;

class MultiFilter {
public:
	MultiFilter(HashManager const &hashManager);
	void addFilter(size_t filterSize, string const &filterID,
			string const &filePath);
	const boost::unordered_map<string, bool> &multiContains(string const &kmer);
	const vector<string> &getFilterIds() const;
	virtual ~MultiFilter();
private:
	//so don't have to reallocated memory multiple times
	boost::unordered_map<string, bool> tempResults;

	boost::unordered_map<string, boost::shared_ptr<BloomFilter> > filters;
	HashManager hashMan;
	vector<string> filterIDs;
};

#endif /* MULTIFILTER_H_ */
