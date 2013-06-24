/*
 * MultiFilter.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: cjustin
 */

#include "MultiFilter.h"

MultiFilter::MultiFilter(vector<string> const &hashFn,
		vector<size_t> const &seeds) :
		hashMan(HashManager())
{
	vector<size_t>::const_iterator seedsItr = seeds.begin();
	//Create HashManager for MultiFilter
	for (vector<string>::const_iterator hashFnItr = hashFn.begin();
			hashFnItr != hashFn.end(); ++hashFnItr)
	{
		hashMan.addHashFunction(*hashFnItr, *seedsItr);
		++seedsItr;
	}

}

void MultiFilter::addFilter(size_t filterSize, string const &filterID,
		string const &filePath)
{
	boost::shared_ptr<BloomFilter> filter(
			new BloomFilter(filterSize, filePath, hashMan));
	filters[filterID] = filter;
	filterIDs.push_back(filterID);
}

//todo: implement partial hash function hashing (ie. Only half the number of hashing values for one filter)
/*
 * checks filters for kmer, hashing only single time
 */
const boost::unordered_map<string, bool> &MultiFilter::multiContains(
		string const &kmer)
{
	const vector<size_t> &hashResults = hashMan.multiHash(kmer);

	for (boost::unordered_map<string, boost::shared_ptr<BloomFilter> >::iterator it =
			filters.begin(); it != filters.end(); ++it)
	{
		tempResults[(*it).first] = ((*it).second)->contains(hashResults);
	}
	return tempResults;
}

/*
 * checks filters for kmer, given a list of filterIDs, hashing only single time
 */
const boost::unordered_map<string, bool> &MultiFilter::multiContains(
		string const &kmer, vector<string> const &tempFilters)
{
	const vector<size_t> &hashResults = hashMan.multiHash(kmer);

	for (vector<string>::const_iterator it = tempFilters.begin();
			it != tempFilters.end(); ++it)
	{
		tempResults[*it] = (filters[*it])->contains(hashResults);
	}

	return tempResults;
}

const vector<string> &MultiFilter::getFilterIds() const
{
	return filterIDs;
}

MultiFilter::~MultiFilter()
{
}

