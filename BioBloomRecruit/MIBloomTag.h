/*
 * MIBloomTag.h
 *
 *  Created on: Jan 15, 2017
 *      Author: cjustin
 */

#ifndef MIBLOOMTAG_H_
#define MIBLOOMTAG_H_

#include "bloomfilter/MIBloomFilter.hpp"
#include <vector>
#include <string>
#include <google/dense_hash_map>

class MIBloomTag {
public:
	typedef typename google::dense_hash_map<ID, ID> IDMap;
	explicit MIBloomTag(vector<string> const &filenames, size_t numElements, IDMap idMap);

	MIBloomFilter generate(const string &filePrefix, double fpr);

private:
	size_t m_totalEntries;
	size_t m_expectedEntries;
	vector<string> m_filenames;

};

#endif /* MIBLOOMTAG_H_ */
