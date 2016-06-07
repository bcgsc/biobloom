/*
 * SimMat.h
 *
 *  Created on: Apr 13, 2015
 *      Author: cjustin
 */

#ifndef SIMMAT_H_
#define SIMMAT_H_

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cmath>
#include <stdint.h>
#include <string>
#include <vector>
#include <unistd.h>
#include <assert.h>
#include <omp.h>
#include "SpacedSeedIndex.h"
#include <google/dense_hash_map>
#include <google/dense_hash_set>

using namespace std;

template<typename T, typename ID>
class SimMat {
public:
	typedef typename google::dense_hash_map<ID, ID> IDMap;
	typedef typename google::dense_hash_set<ID> IDSet;

	SimMat(const SpacedSeedIndex<ID> &index, unsigned numRead) :
			m_numRead(numRead), m_simMat(
					vector<T>(m_numRead * (m_numRead - 1) / 2)), m_indexTbl(
					index) {
		cerr << "Size of Similarity Matrix: "
				<< m_numRead * (m_numRead - 1) * sizeof(T) << " bytes"
				<< endl;
		getsimMat(m_simMat, m_indexTbl);
	}

//	inline void writeSimMat() {
//		outsimMat(m_file, m_fSimMat, m_rSimMat);
//	}

	/*
	 * Note: IDs in matrix are offset by 1 (+1 to get true ID)
	 */
	//TODO refactor to simplify code
	inline vector<boost::shared_ptr<google::dense_hash_map<ID, ID> > > getGroupings(
			double threshold, std::ofstream &file) {
		assert(threshold);

		//ID mapping to its associated collision ID
		google::dense_hash_map<ID, boost::shared_ptr<IDSet> > colliIDs;
		colliIDs.set_empty_key(opt::EMPTY);
		vector<boost::shared_ptr<IDMap> > groupIndex(m_numRead + 1);
		for (size_t readID = 0; readID < m_numRead + 1; ++readID) {
			groupIndex[readID] = boost::shared_ptr<IDMap>(new IDMap());
			groupIndex[readID]->set_empty_key(opt::EMPTY);
		}

		//first level
		ID lastID = m_numRead + 1;

//		while (lastID != opt::COLLI) {
//			vector<double> normSim(m_simMat.size());
//			//convert counts into matrix
//			for (ID readID1 = 0; readID1 < m_numRead; ++readID1) {
//				for (ID readID2 = 0; readID2 < m_numRead; ++readID2) {
//					size_t simInd = getSimIdx(std::max(readID1, readID2),
//							std::min(readID1, readID2));
//					double minSize = min(m_indexTbl.getUnique(readID1),
//							m_indexTbl.getUnique(readID2));
//					if (readID1 != readID2) {
//						normSim[simInd] = m_simMat[simInd] / minSize;
//					}
//				}
//			}
//			//generate ordered list of similar elements
//			vector<T> indexOrder = sort_indexes(normSim);
//
//			for (vector<T>::const_iterator itr = indexOrder.begin();
//					itr != indexOrder.end(); ++itr) {
//				//assign collision ID
//				ID readID1;
//				ID readID2;
//				getReadIdx(*itr, readID1, readID2);
//				//create collision IDs for most similar elements
//				//if already joined
//				if (normSim[*itr] > threshold) {
//					boost::shared_ptr<IDMap> &currMap1 = groupIndex[readID1 + 1];
//					boost::shared_ptr<IDMap> &currMap2 = groupIndex[readID2 + 1];
//					(*currMap1)[readID2 + 1] = lastID;
//					(*currMap2)[readID1 + 1] = lastID;
//					if (colliIDs.find(lastID) == colliIDs.end()) {
//						colliIDs[lastID] = boost::shared_ptr<IDSet>(
//								new IDSet());
//						colliIDs[lastID]->set_empty_key(opt::EMPTY);
//					}
//					colliIDs[lastID]->insert(readID1 + 1);
//					colliIDs[lastID]->insert(readID2 + 1);
//					++lastID;
//					if(lastID == opt::COLLI){
//						cerr << "ran out of IDs" << endl;
//						break;
//					}
//				} else {
//					break;
//				}
//				//end at threshold
//			}
//		}
//		//collapse collision IDs
//		//take average? of all significant collisions
//		//populate list of collapsed IDs


//		//populate normalized matrix
//		for (ID readID1 = 0; readID1 < m_numRead; ++readID1) {
//			for (ID readID2 = 0; readID2 < readID1; ++readID2) {
//				size_t simInd = getSimIdx(readID1, readID2);
//				double minSize = min(m_indexTbl.getReadLength(readID1),
//						m_indexTbl.getReadLength(readID2));
//				double simVal = m_simMat[simInd] / minSize;
//				tempSimMat[readID1][readID2] = simVal;
//			}
//		}
//
//		double bestSim = 0;
//		ID bestRead1 = 0;
//		ID bestRead2 = 0;
//
//		//Find most similar elements
//		for (ID readID1 = 0; readID1 < m_numRead; ++readID1) {
//			for (ID readID2 = 0; readID2 < readID1; ++readID2) {
//				if(bestSim < tempSimMat[readID1][readID2]){
//					bestSim = tempSimMat[readID1][readID2];
//					bestRead1 = readID1;
//					bestRead2 = readID2;
//				}
//			}
//		}
//
//		//merge most similar elements
//		for (ID readID1 = 0; readID1 < tempSimMat.size(); ++readID1) {
//			for (ID readID2 = 0; readID2 < readID1; ++readID2) {
//				if(bestSim < tempSimMat[readID1][readID2]){
//					bestSim = tempSimMat[readID1][readID2];
//					bestRead1 = readID1;
//					bestRead2 = readID2;
//				}
//			}
//		}
//
//
//
//		//create collision ID for these elements
//		//populate groupIndex
//		//populate colliIDs
//		//repeat until at root
//
//		std::ofstream phylip;
//		phylip.open((opt::outputPrefix + ".phylip").c_str());
//
//		phylip << "\t" << m_numRead << endl;
//		for (ID readID1 = 0; readID1 < m_numRead; ++readID1) {
//			phylip << (readID1 + 1);
//			stringstream convert;
//			convert << (readID1 + 1);
//			string id1 = convert.str();
//			for (unsigned i = id1.length(); i < 10; ++i){
//				phylip << " ";
//			}
//			for (ID readID2 = 0; readID2 < m_numRead; ++readID2) {
//				if (readID1 == readID2) {
//					phylip << "\t0";
//				} else {
//					size_t simInd = getSimIdx(std::max(readID1, readID2),
//							std::min(readID1, readID2));
//					double minSize = min(m_indexTbl.getUnique(readID1),
//							m_indexTbl.getUnique(readID2));
//					double simVal = double(m_simMat[simInd]) / minSize;
//					cout << "\t" << simVal;
//				}
//			}
//			phylip << endl;
//		}
//		file.close();

		//for each element compute similarity
		for (size_t readID1 = 1; readID1 < m_numRead; ++readID1) {
			for (unsigned readID2 = 0; readID2 < readID1; ++readID2) {
				size_t simInd = getSimIdx(readID1, readID2);
				double minSize = min(m_indexTbl.getUnique(readID1),
						m_indexTbl.getUnique(readID2));
				double simVal = m_simMat[simInd] / minSize;
				//if element has sufficient similarity
				if (simVal > threshold) {
					boost::shared_ptr<IDMap> &currMap1 = groupIndex[readID1 + 1];
					boost::shared_ptr<IDMap> &currMap2 = groupIndex[readID2 + 1];
					if(currMap1->size() > 0){
						ID candiateID = 0;
						double maxScore = 0;
						//check against all current matching IDs
						for (typename IDMap::const_iterator itr = currMap1->begin();
								itr != currMap1->end(); ++itr) {
							ID maxID = readID2;
							ID minID = itr->first - 1;
							if (maxID != minID) {
								if (maxID < minID) {
									swap(maxID, minID);
								}
								double min = std::min(
										m_indexTbl.getUnique(maxID),
										m_indexTbl.getUnique(minID));
								double score = m_simMat[getSimIdx(maxID, minID)]
										/ min;
								if (score > threshold && maxScore < score) {
									maxScore = score;
									candiateID = itr->second;
								}
							}
						}
						if(maxScore == 0){
							continue;
						}
						else{
							//assign new ID to hash map
							(*currMap1)[readID2 + 1] = candiateID;
							(*currMap2)[readID1 + 1] = candiateID;
							colliIDs[candiateID]->insert(readID2 + 1);
						}
					}
					if(currMap2->size() > 0){
						ID candiateID = 0;
						double maxScore = 0;
						//check against all current matching IDs
						for (typename IDMap::const_iterator itr = currMap2->begin();
								itr != currMap2->end();
								++itr) {
							ID maxID = readID1;
							ID minID = itr->first - 1;
							if (maxID != minID) {
								if (maxID < minID) {
									swap(maxID, minID);
								}
								double min = std::min(
										m_indexTbl.getUnique(maxID),
										m_indexTbl.getUnique(minID));
								double score = m_simMat[getSimIdx(maxID, minID)]
										/ min;
								if (score > threshold && maxScore < score) {
									maxScore = score;
									candiateID = itr->second;
								}
							}
						}
						if(maxScore == 0){
							continue;
						}
						else{
							//assign new ID to hash map
							(*currMap1)[readID2 + 1] = candiateID;
							(*currMap2)[readID1 + 1] = candiateID;
							colliIDs[candiateID]->insert(readID1 + 1);
						}
					}
					//new collision ID
					if(currMap1->find(readID2 + 1) == currMap1->end()){
						(*currMap1)[readID2 + 1] = lastID;
						(*currMap2)[readID1 + 1] = lastID;
						if (colliIDs.find(lastID) == colliIDs.end()) {
							colliIDs[lastID] = boost::shared_ptr<IDSet>(
									new IDSet());
							colliIDs[lastID]->set_empty_key(opt::EMPTY);
						}
						colliIDs[lastID]->insert(readID1 + 1);
						colliIDs[lastID]->insert(readID2 + 1);
						++lastID;
						assert(lastID != opt::COLLI);
					}
				}
			}
		}
		assert(file);
		groupIndex.resize(groupIndex.size()+colliIDs.size());
		for (typename google::dense_hash_map<ID, boost::shared_ptr<IDSet> >::iterator itr =
				colliIDs.begin(); itr != colliIDs.end(); ++itr) {
			IDSet idSet = *(itr->second);
			file << itr->first;
			for (typename IDSet::iterator idItr = idSet.begin();
					idItr != idSet.end(); ++idItr) {
				file << "\t" << *idItr;
			}
			file << endl;
			//insert self (prevent uninitialized values)
			groupIndex[itr->first] = boost::shared_ptr<IDMap>(new IDMap());
			groupIndex[itr->first]->set_empty_key(opt::EMPTY);
			//collision with collision ID yield the collision if ID is part of collision ID
			for (typename IDSet::iterator idItr = idSet.begin();
					idItr != idSet.end(); ++idItr) {
				(*groupIndex[itr->first])[*idItr] = itr->first;
			}
		}
		return (groupIndex);
	}

	/*
	 * For sorting 2 vectors based on first vector
	 */
	inline void sort(vector<unsigned> &counts, vector<ID> &readIDs, int l,
			int r) {
		int i = l - 1, j = r;
		unsigned v = counts[r];
		unsigned tmp_id, tmp_loc;
		if (r <= l)
			return;
		for (;;) {
			while (counts[++i] < v)
				;
			while (v < counts[--j])
				if (j == l)
					break;
			if (i >= j)
				break;
			tmp_id = counts[i];
			counts[i] = counts[j];
			counts[j] = tmp_id;
			tmp_loc = readIDs[i];
			readIDs[i] = readIDs[j];
			readIDs[j] = tmp_loc;
		}
		tmp_id = counts[i];
		counts[i] = counts[r];
		counts[r] = tmp_id;

		tmp_loc = readIDs[i];
		readIDs[i] = readIDs[r];
		readIDs[r] = tmp_loc;

		sort(counts, readIDs, l, i - 1);
		sort(counts, readIDs, i + 1, r);
	}

	virtual ~SimMat() {
	}
private:

	unsigned m_numRead;
	vector<T> m_simMat;
	const SpacedSeedIndex<ID> &m_indexTbl;

	inline void getFname(const char *filename, std::string &bName,
			std::string &eName) {
		std::string fName(filename);
		size_t pos = fName.rfind(".");
		if (pos == std::string::npos) {
			bName = fName;
			eName = "";
			return;
		}
		if (pos == 0) {
			bName = "";
			eName = fName.substr(pos + 1, std::string::npos);
			return;
		}
		bName = fName.substr(0, pos);
		eName = fName.substr(pos + 1, std::string::npos);
	}

	inline void getHeader(const char *aName, vector<unsigned> &lengthSeq) {
		string faSeq, faHead;
		ifstream rFile(aName);
		unsigned readId = 0;
		while (getline(rFile, faHead)) {
			getline(rFile, faSeq);
			lengthSeq.push_back(faSeq.length());
			++readId;
		}
		rFile.close();
	}

	inline void getHeader(const char *aName, vector<string> &headerSeq,
			vector<unsigned> &lengthSeq) {
		string faSeq, faHead;
		ifstream rFile(aName);
		unsigned readId = 0;
		while (getline(rFile, faHead)) {
			getline(rFile, faSeq);
			headerSeq.push_back(faHead.substr(1, string::npos));
			lengthSeq.push_back(faSeq.length());
			++readId;
		}
		rFile.close();
	}

	inline void getsimMat(vector<T> &simMat, const SpacedSeedIndex<ID> &indTable) {
		for (size_t i = 0; i < m_numRead * (m_numRead - 1) / 2; i++)
			simMat[i] = 0;
		size_t i;
#pragma omp parallel for shared(indTable, simMat) private(i) schedule(dynamic)
		for (i = 0; i < indTable.size(); i++) {
			if (indTable.at(i).size > 1) {
				for (ID j = 0; j < (indTable.at(i).size - 1); j++) {
					for (ID k = j + 1; k < indTable.at(i).size; k++) {
						unsigned minInd = indTable.at(i).id[j], maxInd =
								indTable.at(i).id[k];
						if (minInd == maxInd)
							continue;
						if (minInd > maxInd) {
							std::swap(minInd,maxInd);
						}
						size_t cellIdx = getSimIdx(maxInd, minInd);
#pragma omp atomic
						++simMat[cellIdx];
					}
				}
			}
		}
	}

	//readID1 must be > readID2
	inline size_t getSimIdx(ID readID1, ID readID2){
		assert(readID1 > readID2);
		return (readID1 * (readID1 - 1) / 2 + readID2);
	}

	//read1 > read2
	inline void getReadIdx(size_t idx, ID &read1, ID &read2){
		read1 = trunc((sqrt(1 + 4*2*idx)-1)/2);
		read2 = idx - read1*(read1+1)/2;
	}

	inline void outsimMat(const char* aName, unsigned *simMat) {
		string bName, eName;
		getFname(aName, bName, eName);
		ostringstream outss;
		outss << bName << ".smat";
		ofstream outDist(outss.str().c_str());
		for (unsigned i = 1; i < m_numRead; i++) {
			outDist << i;
			for (unsigned j = 0; j < i; j++)
				outDist << "\t"
						<< simMat[i * (i - 1) / 2 + j];
			outDist << "\n";
		}
		outDist.close();
		cerr << "Similarity matrix was written in: " << outss.str() << "\n";
	}

	vector<T> sort_indexes(const vector<double> &v) {
	  // initialize original index locations
	  vector<T> idx(v.size(),0);
	  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

	  // sort indexes based on comparing values in v
	  sort(idx.begin(), idx.end(),
	       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

	  return idx;
	}
};

#endif /* SIMMAT_H_ */
