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
#include "bloomfilter/BloomMapSS.hpp"
#include "bloomfilter/BloomMapSSBitVec.hpp"
#include "bloomfilter/RollingHashIterator.h"
#include "Common/Options.h"
#include <sdsl/int_vector.hpp>

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
//	google::dense_hash_map<ID, PairID> m_collisionIn;
//	google::dense_hash_map<PairID, ID> m_collisionOut;
//	google::dense_hash_map<PairID, &std::vector<size_t>> collisionTrack;
//
//	//TODO MAKE INTO OPTION
//	unsigned m_collisionTheshold = 10000;
//	unsigned m_countTheshold = 1000;

//	size_t m_collisionCount;
//	ID m_currentID;

	//helper methods
	inline void loadSeq(BloomMapSSBitVec<ID> &bloomMap, const string& seq, ID value) {
		/* init rolling hash state and compute hash values for first k-mer */
		for (RollingHashIterator itr(seq, m_kmerSize, bloomMap.getSeedValues());
				itr != itr.end(); ++itr) {
			bloomMap.insert(*itr, value);
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
			for (size_t i = 0; i < itr->size(); ++i) {
				size_t pos = itr->at(i) % bv.size();
				uint64_t *dataIndex = bv.data() + (pos >> 6);
				uint64_t bitMaskValue = (uint64_t) 1 << (pos & 0x3F);
				count = __sync_fetch_and_or(dataIndex, bitMaskValue) >> (pos & 0x3F) &1;
			}
		}
		return count;
	}

	inline void writeIDs(const string& filename, google::dense_hash_map<string,ID> headerIDs) {
		std::ofstream file;
		file.open(filename.c_str());
		assert(file);
		for (google::dense_hash_map<string, ID>::iterator itr =
				headerIDs.begin(); itr != headerIDs.end(); ++itr) {
			file << itr->second << "\t" << itr->first << "\n";
			assert(file);
		}
		file.close();
	}
};

#endif /* BLOOMMAPGENERATOR_H_ */
