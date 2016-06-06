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
#include <google/dense_hash_set>
#include "bloomfilter/BloomMapSSBitVec.hpp"
#include "bloomfilter/RollingHashIterator.h"
#include "Common/Options.h"
#include <sdsl/int_vector.hpp>
#include <boost/shared_ptr.hpp>
#include "bloomfilter/BloomFilter.hpp"
#include "tree.hh"
#include "newick_file.hh"
#include <queue>

using namespace std;

class BloomMapGenerator {
public:
	explicit BloomMapGenerator(vector<string> const &filenames,
			unsigned kmerSize, size_t numElements);

	void generate(const string &filePrefix, double fpr);

	virtual ~BloomMapGenerator();
private:

	typedef typename google::dense_hash_set<ID> IDSet;
	typedef typename google::dense_hash_map<ID, ID> IDMap;

	unsigned m_kmerSize;
	size_t m_expectedEntries;
	size_t m_totalEntries;
	vector<string> m_fileNames;
	google::dense_hash_map<string, ID> m_headerIDs;
	vector<boost::shared_ptr<IDMap > > m_colliIDs;

//	//TODO MAKE INTO OPTION
	double m_colliThresh = 0.15;

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

	/*
	 * Traverses nodes, assigning an ID to each node
	 * Higher nodes have a higher ID value
	 */
	//TODO: figure out if branch lengths matter
	void assignIDBFS(BiRC::treelib::Tree &tree, ID startID) {
		queue<int> travQueue;
		travQueue.push(tree.root());
		while (!travQueue.empty()) {
			int currentID = travQueue.front();
			travQueue.pop();
			if (tree.label(currentID).empty()) {
				std::stringstream ss;
				ss << startID--;
				tree.setLabel(currentID, std::string(ss.str()));
			}
			if (tree.left_child(currentID) != -1) {
				travQueue.push(tree.left_child(currentID));
			}
			if (tree.right_child(currentID) != -1) {
				travQueue.push(tree.right_child(currentID));
			}
		}
	}

	//TODO: optimize me!
	void setColliIds(const BiRC::treelib::Tree &tree, ID numNodes,
			std::ofstream &file) {
		google::dense_hash_map<ID, boost::shared_ptr<IDSet> > colliIDs;
		colliIDs.set_empty_key(opt::EMPTY);
		m_colliIDs.resize(numNodes);

		ID numLeaf = (numNodes + 1)/2;

		for (ID readID = 0; readID < numNodes; ++readID) {
			m_colliIDs[readID] = boost::shared_ptr<IDMap>(new IDMap());
			m_colliIDs[readID]->set_empty_key(opt::EMPTY);
		}

		ID root = stoi(tree.label(tree.root()));

		//for every possible combinination of IDs
		//do not bother assigning when root node is involved
		for (ID i = 1; i < numNodes - 1; ++i) {
			for (ID j = i + 1; j < numNodes - 1; ++j) {
				ID common = getLCAID(tree, i, j);
				if (common != root) {
					if (colliIDs.find(common) == colliIDs.end()) {
						colliIDs[common] = boost::shared_ptr<IDSet>(
								new IDSet());
						colliIDs[common]->set_empty_key(opt::EMPTY);
					}
					colliIDs[common]->insert(i);
					colliIDs[common]->insert(j);
					(*m_colliIDs[i])[j] = common;
					(*m_colliIDs[j])[i] = common;
				}
			}
		}

		assert(file);
		for (typename google::dense_hash_map<ID, boost::shared_ptr<IDSet> >::iterator itr =
				colliIDs.begin(); itr != colliIDs.end(); ++itr) {
			IDSet idSet = *(itr->second);
			file << itr->first;
			for (typename IDSet::iterator idItr = idSet.begin();
					idItr != idSet.end(); ++idItr) {
				if(*idItr <= numLeaf)
					file << "\t" << *idItr;
			}
			file << endl;
		}
	}

	/*
	 * Traverses tree from current positions to find LCA
	 */
	//TODO: optimize me!
	ID getLCAID(const BiRC::treelib::Tree &tree, ID id1, ID id2) {
		std::stringstream ss1;
		ss1 << id1;
		std::stringstream ss2;
		ss2 << id2;
		int nodeID1 = tree.node(std::string(ss1.str()));
		int nodeID2 = tree.node(std::string(ss2.str()));

		google::dense_hash_set<string> foundIDs;
		foundIDs.set_empty_key("");
		int currParentID = tree.parent(nodeID1);
		while(currParentID != -1){
			foundIDs.insert(tree.label(currParentID));
			currParentID = tree.parent(currParentID);
		}
		currParentID = tree.parent(nodeID2);
		while(currParentID != -1){
			if(foundIDs.find(tree.label(currParentID)) != foundIDs.end()){
				return stoi(tree.label(currParentID));
			}
			currParentID = tree.parent(currParentID);
		}
		return stoi(tree.label(tree.root()));
	}

	/*
	 * checks if file exists
	 */
	inline bool fexists(const string &filename)
	{
		ifstream ifile(filename.c_str());
		return ifile.good();
	}
};

#endif /* BLOOMMAPGENERATOR_H_ */
