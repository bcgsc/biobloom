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
			m_indexTbl(index), m_numRead(numRead) {
		cerr << "Size of Similarity Matrix: "
				<< m_numRead * (m_numRead - 1) * sizeof(T) << " bytes"
				<< endl;
		m_simMat = (T *) malloc(
				m_numRead * (m_numRead - 1) / 2 * sizeof(T));
		getsimMat(m_simMat, m_indexTbl);
	}

//	inline void writeSimMat() {
//		outsimMat(m_file, m_fSimMat, m_rSimMat);
//	}

	/*
	 * Note: IDs in matrix are offset by 1 (+1 to get true ID)
	 */
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
		ID lastID = m_numRead + 1;
		//for each element compute similarity
		for (size_t readID1 = 0; readID1 < m_numRead; ++readID1) {
			for (unsigned readID2 = 0; readID2 < readID1; ++readID2) {
				size_t simInd = getSimIdx(readID1, readID2);
				double ikMean = sqrt(
						m_indexTbl.getReadLength(readID1)
								* m_indexTbl.getReadLength(readID2));
				double simVal =	m_simMat[simInd];
//				cout << simVal / ikMean << "\t" << simVal << "\t"
//						<< m_indexTbl.getReadLength(readID1) << "\t"
//						<< m_indexTbl.getReadLength(readID2) << "\t"
//						<< readID1 + 1 << "\t" << readID2 + 1 << endl;
				//if element has sufficient similarity
				if (simVal / ikMean > threshold) {
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
							if (maxID < minID) {
								swap(maxID, minID);
							}
							double mean = sqrt(
									m_indexTbl.getReadLength(maxID)
											* m_indexTbl.getReadLength(minID));
							double score = m_simMat[getSimIdx(maxID, minID)]/mean;
							if (score > threshold && maxScore < score) {
								maxScore = score;
								candiateID = itr->second;
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
								itr != currMap2->end(); ++itr) {
							ID maxID = readID1;
							ID minID = itr->first - 1;
							if (maxID < minID) {
								swap(maxID, minID);
							}
							double mean = sqrt(
									m_indexTbl.getReadLength(maxID)
											* m_indexTbl.getReadLength(minID));
							double score = m_simMat[getSimIdx(maxID, minID)]/mean;
							if (score > threshold && maxScore < score) {
								maxScore = score;
								candiateID = itr->second;
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
			//collision with collision ID yield the collision if ID is part of collision ID
			for (typename IDSet::iterator idItr = idSet.begin();
					idItr != idSet.end(); ++idItr) {
				(*groupIndex[itr->first])[*idItr] = itr->first;
			}
		}
		return (groupIndex);
	}

	inline size_t getSimIdx(ID readID1, ID readID2){
		return (readID1 * (readID1 - 1) / 2 + readID2);
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

//	/*
//	 * Assumes thresholds are sorted from high to low
//	 */
//	inline vector<boost::shared_ptr<google::dense_hash_map<ID, ID> > > getGroupings(
//			const vector<double> &thresholds, std::ofstream &file) {
//		//ID mapping to its associated collision ID
//		google::dense_hash_map<ID, boost::shared_ptr<IDSet> > colliIDs;
//		colliIDs.set_empty_key(opt::EMPTY);
//		vector<boost::shared_ptr<IDMap> > groupIndex(m_numRead + 1);
//		for (size_t readID = 0; readID < m_numRead + 1; ++readID) {
//			groupIndex[readID] = boost::shared_ptr<IDMap>(new IDMap());
//			groupIndex[readID]->set_empty_key(opt::EMPTY);
//		}
//		ID lastID = m_numRead + 1;
//		ID start = lastID;
//		//for each element compute similarity
//		for (unsigned i = 0; i < thresholds.size(); ++i) {
//			cerr << "computing grouping for similarity threshold: "
//					<< thresholds[i] << endl;
//			ID start = lastID;
//			double currentThreshold = thresholds[i];
//			//subgroup index
//			vector<boost::shared_ptr<IDMap> > subGroupIndex(m_numRead + 1);
//			for (size_t readID1 = 0; readID1 < m_numRead; ++readID1) {
//				for (unsigned readID2 = 0; readID2 < readID1; ++readID2) {
//					size_t simInd = getSimIdx(readID1, readID2);
//					double ikMean = sqrt(
//							m_indexTbl.getReadLength(readID1)
//									* m_indexTbl.getReadLength(readID2));
//					double simVal = m_simMat[simInd] / ikMean;
//					//if element has sufficient similarity
//					if (simVal > threshold) {
//						boost::shared_ptr<IDMap> &currMap1 = groupIndex[readID1
//								+ 1];
//						boost::shared_ptr<IDMap> &currMap2 = groupIndex[readID2
//								+ 1];
//						if (currMap1->size() > 0) {
//							ID candiateID = 0;
//							double maxScore = 0;
//							//check against all current matching IDs
//							for (typename IDMap::const_iterator itr =
//									currMap1->begin(); itr != currMap1->end();
//									++itr) {
//								ID maxID = readID2;
//								ID minID = itr->first - 1;
//								if (maxID < minID) {
//									swap(maxID, minID);
//								}
//								double mean = sqrt(
//										m_indexTbl.getReadLength(maxID)
//												* m_indexTbl.getReadLength(
//														minID));
//								double score = m_simMat[getSimIdx(maxID, minID)]
//										/ mean;
//								if (score > threshold && maxScore < score) {
//									maxScore = score;
//									candiateID = itr->second;
//								}
//							}
//							if (maxScore == 0) {
//								continue;
//							} else {
//								//assign new ID to hash map
//								(*currMap1)[readID2 + 1] = candiateID;
//								(*currMap2)[readID1 + 1] = candiateID;
//								colliIDs[candiateID]->insert(readID2 + 1);
//							}
//						}
//						if (currMap2->size() > 0) {
//							ID candiateID = 0;
//							double maxScore = 0;
//							//check against all current matching IDs
//							for (typename IDMap::const_iterator itr =
//									currMap2->begin(); itr != currMap2->end();
//									++itr) {
//								ID maxID = readID1;
//								ID minID = itr->first - 1;
//								if (maxID < minID) {
//									swap(maxID, minID);
//								}
//								double mean = sqrt(
//										m_indexTbl.getReadLength(maxID)
//												* m_indexTbl.getReadLength(
//														minID));
//								double score = m_simMat[getSimIdx(maxID, minID)]
//										/ mean;
//								if (score > threshold && maxScore < score) {
//									maxScore = score;
//									candiateID = itr->second;
//								}
//							}
//							if (maxScore == 0) {
//								continue;
//							} else {
//								//assign new ID to hash map
//								(*currMap1)[readID2 + 1] = candiateID;
//								(*currMap2)[readID1 + 1] = candiateID;
//								colliIDs[candiateID]->insert(readID1 + 1);
//							}
//						}
//						//new collision ID
//						if (currMap1->find(readID2 + 1) == currMap1->end()) {
//							(*currMap1)[readID2 + 1] = lastID;
//							(*currMap2)[readID1 + 1] = lastID;
//							if (colliIDs.find(lastID) == colliIDs.end()) {
//								colliIDs[lastID] = boost::shared_ptr<IDSet>(
//										new IDSet());
//								colliIDs[lastID]->set_empty_key(opt::EMPTY);
//							}
//							colliIDs[lastID]->insert(readID1 + 1);
//							colliIDs[lastID]->insert(readID2 + 1);
//							++lastID;
//							assert(lastID != opt::COLLI);
//						}
//					}
//				}
//			}
//		}
//		assert(file);
//		for (typename google::dense_hash_map<ID, boost::shared_ptr<IDSet> >::iterator itr =
//				colliIDs.begin(); itr != colliIDs.end(); ++itr) {
//			IDSet idSet = *(itr->second);
//			file << itr->first;
//			for (typename IDSet::iterator idItr = idSet.begin();
//					idItr != idSet.end(); ++idItr) {
//				file << "\t" << *idItr;
//			}
//			file << endl;
//		}
//		return (groupIndex);
//	}

//	inline void getCandidates(
//			vector<boost::shared_ptr<vector<Link> > > &candidates) {
//		vector<unsigned> counts(m_numRead, 0);
//		vector<ID> readIDs(m_numRead, 0);
//		//get list of count for each element
////#pragma omp parallel for shared(counts,readIDs) schedule(dynamic)
//		for (size_t i = 0; i < m_numRead; i++) {
//			readIDs[i] = i;
//			candidates[i] = boost::shared_ptr<vector<Link> >(
//					new vector<Link>());
//			for (unsigned j = 0; j < i; j++) {
//				// geometric mean
//				double ikMean = sqrt(
//						m_indexTbl.getReadLength(i)
//								* m_indexTbl.getReadLength(j));
//				size_t simInd = i * (i - 1) / 2 + j;
//				double simVal = m_fSimMat[simInd];
//				if (m_fSimMat[simInd] < m_rSimMat[simInd]) {
//					simVal = m_rSimMat[simInd];
//				}
//				if (simVal / ikMean > opt::thresh) {
////#pragma omp atomic
//					++counts[i];
////#pragma omp atomic
//					++counts[j];
//				}
//			}
//		}
//
//		//TODO Need Parallel sort?
//		//sort counts;
//		sort(counts, readIDs, 0, counts.size());
//
////#pragma omp parallel for shared(candidates, readIDs) schedule(dynamic)
//		for (size_t i = 0; i < readIDs.size(); i++) {
//			size_t readID1 = readIDs[i];
//			for (ID readID2 = 0; readID2 < readID1; ++readID2) {
//				// geometric mean
//				double ikMean = sqrt(
//						m_indexTbl.getReadLength(readID1)
//								* m_indexTbl.getReadLength(readID2));
//				size_t simInd = readID1 * (readID1 - 1) / 2 + readID2;
////#pragma omp critical(simVal)
//				{
//					double simVal =
//							m_fSimMat[simInd] > m_rSimMat[simInd] ?
//									m_fSimMat[simInd] : m_rSimMat[simInd];
//					assert(simVal < ikMean);
//					Direction forward =
//							m_fSimMat[simInd] > m_rSimMat[simInd] ? FW : RV;
//					if (simVal / ikMean > opt::thresh) {
//						//insert values
//						candidates[readID1]->push_back(Link(readID2, forward));
//						//set location to 0 (i.e. marked), both forward and reverse
//						m_fSimMat[simInd] = 0;
//						m_rSimMat[simInd] = 0;
//					}
//				}
//			}
//			for (size_t readID2 = readID1 + 1; readID2 < readIDs.size();
//					++readID2) {
//				// geometric mean
//				double ikMean = sqrt(
//						m_indexTbl.getReadLength(readID1)
//								* m_indexTbl.getReadLength(readID2));
//				size_t simInd = readID2 * (readID2 - 1) / 2 + readID1;
////#pragma omp critical(simVal)
//				{
//					double simVal =
//							m_fSimMat[simInd] > m_rSimMat[simInd] ?
//									m_fSimMat[simInd] : m_rSimMat[simInd];
//					assert(simVal < ikMean);
//					Direction forward =
//							m_fSimMat[simInd] > m_rSimMat[simInd] ? FW : RV;
//					if (simVal / ikMean > opt::thresh) {
//						//insert values
//						candidates[readID1]->push_back(Link(readID2, forward));
//						//set location to 0 (i.e. marked), both forward and reverse
//						m_fSimMat[simInd] = 0;
//						m_rSimMat[simInd] = 0;
//					}
//				}
//			}
//		}
//	}

//	/*
//	 * Obtains candidates counts given a threshold
//	 */
//	inline unsigned getCandidateCounts(vector<double> &normSSs, double thresh) {
//		unsigned count = 0;
//		for (unsigned i = 1; i < m_numRead; i++) {
//			for (unsigned j = 0; j < i; j++) {
//				double ikMean = sqrt(
//						m_indexTbl.getReadLength(i)
//								* m_indexTbl.getReadLength(j)); // geometric mean
//				unsigned simInd = i * (i - 1) / 2 + j;
//				double simVal = m_fSimMat[simInd];
//				if (m_fSimMat[simInd] < m_rSimMat[simInd]) {
//					simVal = m_rSimMat[simInd];
//				}
//				double normSS = simVal / ikMean;
//				if (normSS > thresh) {
//					normSSs.push_back(normSS);
//					++count;
//				}
//			}
//		}
//		return count;
//	}

//	inline void outGraph(const char*aName) {
//		vector<string> headerSeq;
//		vector<unsigned> lengthSeq;
//		getHeader(aName, headerSeq, lengthSeq);
//
//		string bName, eName;
//		getFname(aName, bName, eName);
//		ostringstream outss;
//		outss << bName << "-graph.gv";
//		ofstream outFile(outss.str().c_str());
//
//		outFile << "digraph g{\n";
//		uint64_t *maxVec = (uint64_t *) malloc(
//				((m_numRead + 63) / 64) * sizeof(uint64_t));
//		for (unsigned i = 0; i < m_numRead; i++) {
//			for (unsigned ii = 0; ii < (m_numRead + 63) / 64; ii++)
//				maxVec[ii] = 0;
//			maxVec[i / 64] |= ((uint64_t) 1 << (63 - i % 64));
//			while (true) {
//				unsigned maxInd = 0;
//				double max = 0.0;
//				for (unsigned k = 0; k < m_numRead; k++) {
//					if ((maxVec[k / 64] & ((uint64_t) 1 << (63 - k % 64)))
//							== 0) {
//						double ikMean = sqrt(lengthSeq[i] * lengthSeq[k]); // geometric mean
//						unsigned simInd = std::max(i, k) * (std::max(i, k) - 1)
//								/ 2 + std::min(i, k);
//						double simVal = std::max(m_fSimMat[simInd],
//								m_rSimMat[simInd]);
//						if (simVal / ikMean > max) {
//							maxInd = k;
//							max = simVal / ikMean;
//						}
//					}
//				}
//				if (max < opt::thresh)
//					break;
//				outFile << "\"" << headerSeq[i] << "\"" << "->" << "\""
//						<< headerSeq[maxInd] << "\"\n";
//				maxVec[maxInd / 64] |= ((uint64_t) 1 << (63 - maxInd % 64));
//			}
//		}
//		outFile << "}\n";
//		outFile.close();
//		cerr << "Overlap graph was written in: " << outss.str() << "\n";
//	}

	virtual ~SimMat() {
		free(m_simMat);
	}
private:

	T *m_simMat;
	const SpacedSeedIndex<ID> &m_indexTbl;

	unsigned m_numRead;

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

	inline void getsimMat(T *simMat, const SpacedSeedIndex<ID> &indTable) {
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
						size_t cellIdx = maxInd * (maxInd - 1) / 2 + minInd;
#pragma omp atomic
						++simMat[cellIdx];
					}
				}
			}
		}
	}

//	void rgetsimMat(T *simMat, const SpacedSeedIndex<ID> &indTable) {
//		for (size_t i = 0; i < m_numRead * (m_numRead - 1) / 2; i++)
//			simMat[i] = 0;
//		size_t i;
//#pragma omp parallel for shared(indTable, simMat) private(i) schedule(dynamic)
//		for (i = 0; i < indTable.size(); i++) {
//			if (indTable.at(i).rsize >= 1 && indTable.at(i).fsize >= 1) {
//				unsigned *sid = (unsigned *) calloc(
//						std::min(indTable.at(i).rsize, indTable.at(i).fsize),
//						sizeof(unsigned));
//				size_t sInd = 0;
//				for (size_t j = 0; j < indTable.at(i).rsize; j++) {
//					for (size_t k = 0; k < indTable.at(i).fsize; k++) {
//						size_t minInd = indTable.at(i).rid[j], maxInd =
//								indTable.at(i).fid[k];
//						if (minInd == maxInd) {
//							sid[sInd++] = maxInd;
//							continue;
//						}
//						if (minInd > maxInd) {
//							size_t tmpInd = minInd;
//							minInd = maxInd;
//							maxInd = tmpInd;
//						}
//#pragma omp atomic
//						++simMat[maxInd * (maxInd - 1) / 2 + minInd];
//					}
//				}
//				if (sInd > 1) {
//					for (size_t j = 0; j < sInd - 1; j++) {
//						for (size_t k = j + 1; k < sInd; k++) {
//							size_t minInd = sid[j], maxInd = sid[k];
//							if (minInd > maxInd) {
//								unsigned tmpInd = minInd;
//								minInd = maxInd;
//								maxInd = tmpInd;
//							}
//							size_t cellIdx = maxInd * (maxInd - 1) / 2 + minInd;
//#pragma omp atomic
//							simMat[cellIdx]--;
//						}
//					}
//				}
//				free(sid);
//			}
//		}
//		for (size_t i = 0; i < m_numRead * (m_numRead - 1) / 2; i++)
//			simMat[i] /= 2;
//	}

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
};

#endif /* SIMMAT_H_ */
