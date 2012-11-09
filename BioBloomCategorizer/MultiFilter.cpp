/*
 * MultiFilter.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: cjustin
 */

#include "MultiFilter.h"

MultiFilter::MultiFilter(HashManager const &hashManager) :
		hashMan(hashManager) {
}

void MultiFilter::addFilter(string const filterID, string const &filePath) {
	boost::shared_ptr<BloomFilter> filter(new BloomFilter(filePath, hashMan));
	filters[filterID] = filter;
	filterIDs.push_back(filterID);
}

//todo: implement partial hash function hashing (ie. Only half the number of hashing  values for one filter)
/*
 * checks filters for kmer, hashing only single time
 */
const boost::unordered_map<string, bool> &MultiFilter::multiContains(
		string const &kmer) {
	const vector<size_t> &hashResults = hashMan.multiHash(kmer);

	for (boost::unordered_map<string, boost::shared_ptr<BloomFilter> >::iterator it =
			filters.begin(); it != filters.end(); ++it) {
		tempResults[(*it).first] = ((*it).second)->contains(hashResults);
	}
	return tempResults;
}

const vector<string> &MultiFilter::getFilterIds() const {
	return filterIDs;
}

MultiFilter::~MultiFilter() {
	// TODO Auto-generated destructor stub
}

