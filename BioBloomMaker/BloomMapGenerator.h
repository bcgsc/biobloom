/*
 * BloomMapGenerator.h
 *
 *  Created on: Mar 17, 2016
 *      Author: cjustin
 */

#ifndef BLOOMMAPGENERATOR_H_
#define BLOOMMAPGENERATOR_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>
#include <google/dense_hash_map>
#include "bloomfilter/BloomMapSSBitVec.hpp"
#include "bloomfilter/RollingHashIterator.h"
#include "Common/Options.h"
#include <sdsl/int_vector.hpp>
#include <boost/shared_ptr.hpp>
#include "bloomfilter/BloomFilter.hpp"

using namespace std;

class BloomMapGenerator {
public:
	explicit BloomMapGenerator(vector<string> const &filenames,
			unsigned kmerSize, size_t numElements);

	void generate(const string &filePrefix, double fpr);

	virtual ~BloomMapGenerator();
private:

	unsigned m_kmerSize;
//	unsigned m_hashNum;
	size_t m_expectedEntries;
	size_t m_totalEntries;
	vector<string> m_fileNames;
	//TODO: replace with vectors?
	google::dense_hash_map<string, ID> m_headerIDs;
	vector<boost::shared_ptr<google::dense_hash_map<ID, ID> > > m_colliIDs;

//
//	//TODO MAKE INTO OPTION
	double m_colliThresh = 0.2;

	inline BloomMapSSBitVec<ID> generateBV(double fpr,
			const vector<vector<unsigned> > &ssVal);

	inline vector<boost::shared_ptr<google::dense_hash_map<ID, ID> > > generateGroups(
			std::ofstream &file);

	//helper methods
	inline void loadSeq(BloomMapSSBitVec<ID> &bloomMap, const string& seq,
			ID value) {
		/* init rolling hash state and compute hash values for first k-mer */
		//TODO FIX -> NEED TO VALIDATE THREAD SAFETY!
		if (opt::colliIDs) {
			for (RollingHashIterator itr(seq, m_kmerSize,
					bloomMap.getSeedValues()); itr != itr.end(); ++itr) {
				bloomMap.insert(*itr, value, m_colliIDs);
			}
		} else {
			for (RollingHashIterator itr(seq, m_kmerSize,
					bloomMap.getSeedValues()); itr != itr.end(); ++itr) {
				bloomMap.insert(*itr, value);
			}
		}
	}

	/*
	 * Returns count of collisions (counts unique k-mers)
	 */
	inline size_t loadSeq(sdsl::bit_vector &bv, const string& seq,
			const vector<vector<unsigned>> &seedVal) {
		size_t count = 0;
		if (seq.size() < m_kmerSize)
			return count;
		/* init rolling hash state and compute hash values for first k-mer */
		for (RollingHashIterator itr(seq, m_kmerSize, seedVal);
				itr != itr.end(); ++itr) {
			unsigned colliCount = 0;
			for (size_t i = 0; i < seedVal.size(); ++i) {
				size_t pos = itr->at(i) % bv.size();
				uint64_t *dataIndex = bv.data() + (pos >> 6);
				uint64_t bitMaskValue = (uint64_t) 1 << (pos & 0x3F);
				colliCount += __sync_fetch_and_or(dataIndex, bitMaskValue)
						>> (pos & 0x3F) & 1;
			}
			if (colliCount == seedVal.size()) {
				++count;
			}
		}
		return count;
	}

	inline void writeIDs(std::ofstream &file, google::dense_hash_map<string,ID> headerIDs) {
		assert(file);
		for (google::dense_hash_map<string, ID>::iterator itr =
				headerIDs.begin(); itr != headerIDs.end(); ++itr) {
			file << itr->second << "\t" << itr->first << "\n";
			assert(file);
		}
	}
};

#endif /* BLOOMMAPGENERATOR_H_ */
