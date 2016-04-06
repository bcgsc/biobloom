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
//#include <google/dense_hash_set>
#include "bloomfilter/BloomMapSS.hpp"
#include "bloomfilter/RollingHashIterator.h"
#include "Common/Options.h"

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
	google::dense_hash_map<ID, string> m_headerIDs;
	google::dense_hash_map<ID, PairID> m_collisionIn;
	google::dense_hash_map<PairID, ID> m_collisionOut;

	size_t m_collisionCount;
	ID m_currentID;

	//helper methods
	//TODO: collision detection
	inline size_t loadSeq(BloomMapSS<ID> &bloomMap, const string& seq, ID value) {
		if (seq.size() < m_kmerSize)
			return 0;
		size_t count = 0;
		//TODO: make into option
		unsigned colliThresh = 0;

		/* init rolling hash state and compute hash values for first k-mer */
		RollingHashIterator itr(seq, m_kmerSize, bloomMap.getSeedValues());
		while (itr != itr.end()) {
			vector<ID> temp = bloomMap.at(*itr);
			google::dense_hash_map<ID, unsigned> tempVect;
			tempVect.set_empty_key(0);
			unsigned maxCount = 0;
			unsigned maxID = 0;

			for (unsigned i = 0; i < temp.size(); ++i) {
				if (temp[i] != 0 && temp[i] != value) {
					google::dense_hash_map<ID, unsigned>::iterator tempItr =
							tempVect.find(temp[i]);
					if (tempItr == tempVect.end()) {
						tempVect[temp[i]] = 1;
					}
					else{
						++tempItr->second;
						if (maxCount < tempVect[temp[i]]) {
							maxCount = tempVect[temp[i]];
							maxID = temp[i];
						}
					}
				}
			}
			//full collision
			if (maxCount >= bloomMap.getSeedValues().size() - colliThresh) {
				//check for if collision ID already exists
				PairID pairID = maxID << ID_BITS | value;
				if(m_collisionOut.find(pairID) == m_collisionOut.end()){
					ID newID =  ++m_currentID;
					assert(newID < opt::COLLI);
					cerr << newID << endl;
					m_collisionIn[newID] = pairID;
					m_collisionOut[pairID] = newID;
					count += !bloomMap.insertAndCheck(*itr, newID, m_collisionCount);
				}
				else{
					count += !bloomMap.insertAndCheck(*itr, m_collisionOut[pairID], m_collisionCount);
				}
//				cout << value;
//				for (unsigned i = 0; i < temp.size(); ++i) {
//					if (temp[i] != 0) {
//						cout << "\t" << temp[i];
//					} else {
//						cout << "\t-";
//					}
//				}
//				cout << endl;
			}
			else{
				count += !bloomMap.insertAndCheck(*itr, value, m_collisionCount);
			}
			++itr;
		}
		return count;
	}

	inline void writeIDs(const string& filename, google::dense_hash_map<ID,string> headerIDs) {
		std::ofstream file;
		file.open(filename.c_str());
		assert(file);
		for (google::dense_hash_map<ID,string>::iterator itr = headerIDs.begin(); itr != headerIDs.end(); ++itr) {
			file << (*itr).first << "\t" << (*itr).second << "\n";
			assert(file);
		}
		file.close();
	}

		
};

#endif /* BLOOMMAPGENERATOR_H_ */
