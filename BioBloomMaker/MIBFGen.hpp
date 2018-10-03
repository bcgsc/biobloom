/*
 * MIBFGen.hpp
 *
 *  Created on: Dec 19, 2017
 *      Author: cjustin
 */

#ifndef CHROMIUMMAP_MIBFGEN_HPP_
#define CHROMIUMMAP_MIBFGEN_HPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>

#include "btl_bloomfilter/MIBloomFilter.hpp"
#include "btl_bloomfilter/stHashIterator.hpp"
#include "btl_bloomfilter/ntHashIterator.hpp"

#include "btl_bloomfilter/BloomFilter.hpp"

#include "Common/Options.h"

#include <tuple>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>
#include <google/dense_hash_set>
#include <sdsl/int_vector.hpp>

#include <zlib.h>
#include <stdio.h>
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "Common/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class MIBFGen {
public:
	MIBFGen(vector<string> const &filenames, unsigned kmerSize,
			size_t numElements = 0) :
			m_kmerSize(kmerSize), m_expectedEntries(numElements), m_totalEntries(
					0), m_fileNames(filenames), m_failedInsert(0) {
		//Instantiate dense hash map
		m_ids.push_back(""); //first entry is empty
		//dense hash maps take POD, and strings need to live somewhere
		m_nameToID.set_empty_key(m_ids[0]);
		size_t counts = 0;

		if (opt::idByFile) {
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				m_ids.push_back(m_fileNames[i].substr(
						m_fileNames[i].find_last_of("/") + 1));
				m_nameToID[m_ids.back()] = m_ids.size() - 1;
			}
#pragma omp parallel for schedule(dynamic)
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				if (fp == NULL) {
					cerr << "file " << m_fileNames[i] << " cannot be opened"
							<< endl;
					exit(1);
				}
				if (!fexists(m_fileNames[i] + ".rv")) {
					cerr << "Missing file " << m_fileNames[i] << ".rv" << endl;
					exit(1);
				}
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Opening " << m_fileNames[i] << endl;
				}
				kseq_t *seq = kseq_init(fp);
				int l;
				for (;;) {
					l = kseq_read(seq);
					if (l >= 0) {
#pragma omp atomic
						counts += seq->seq.l - m_kmerSize + 1;
					} else {
						kseq_destroy(seq);
						break;
					}
				}
				gzclose(fp);
			}
		} else {
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				if (fp == NULL) {
					cerr << "file " << m_fileNames[i] << " cannot be opened"
							<< endl;
					exit(1);
				}
				if (!fexists(m_fileNames[i] + ".rv")) {
					cerr << "Missing file " << m_fileNames[i] << ".rv" << endl;
					exit(1);
				}
				kseq_t *seq = kseq_init(fp);
				int l;
				for (;;) {
					l = kseq_read(seq);
					if (l >= 0) {
						m_ids.push_back(string(seq->name.s, seq->name.l));
						m_nameToID[m_ids.back()] = m_ids.size() - 1;
						counts += seq->seq.l - m_kmerSize + 1;
					} else {
						kseq_destroy(seq);
						break;
					}
				}
				gzclose(fp);
			}
		}

		//make saturation bit is not exceeded
		assert(m_ids.size() < ID(1 << (sizeof(ID) * 8 - 1)));

		//estimate number of k-mers
		if (m_expectedEntries == 0) {
			m_expectedEntries = counts;
		}
		//assume each file one line per file
		if (opt::verbose) {
			cerr << "Expected number of elements: " << m_expectedEntries
					<< endl;
		}
	}

	MIBloomFilter<ID> * generate(const string &filePrefix, double occ) {
		//keep track of time
		double time = omp_get_wtime();


		MIBloomFilter<ID> * miBFBV;
		if (opt::sseeds.empty()) {
			if (opt::verbose)
				cerr << "K-mer mode detected" << endl;
			miBFBV = generateBV(occ);
		} else {
			if (opt::verbose)
				cerr << "Spaced Seeds Detected" << endl;
			vector<vector<unsigned> > ssVal =
					MIBloomFilter<ID>::parseSeedString(opt::sseeds);
			miBFBV = generateBV(occ, &ssVal);
		}
		if (opt::verbose){
			cerr << "Finishing initial Bit vector construction " <<  omp_get_wtime() - time << "s" << endl;
			time = omp_get_wtime();
			cerr << "Populating values of miBF" << endl;
		}
		//first pass
		unsigned j = 1;
		if (opt::verbose)
			cerr << "Pass " << j << endl;
		if (opt::idByFile) {
#pragma omp parallel for schedule(dynamic)
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
//				if (opt::verbose) {
//#pragma omp critical(stderr)
//					cerr << "Opening: "
//							<< (j % 2 ? m_fileNames[i] : m_fileNames[i] + ".rv")
//							<< endl;
//				}
				fp = j % 2 ?
						gzopen(m_fileNames[i].c_str(), "r") :
						gzopen((m_fileNames[i] + ".rv").c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				int l;
				for (;;) {
					string sequence, name;
					{
						l = kseq_read(seq);
						if (l >= 0) {
							sequence = string(seq->seq.s, seq->seq.l);
							name = m_fileNames[i].substr(
									m_fileNames[i].find_last_of("/") + 1);
						}
					}
					if (l >= 0) {
						loadSeq(*miBFBV, name, sequence, j);
					} else {
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
			}
		} else {
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				if (opt::verbose)
					cerr << "Opening: "
							<< (j % 2 ? m_fileNames[i] : m_fileNames[i] + ".rv")
							<< endl;
				fp = j % 2 ?
						gzopen(m_fileNames[i].c_str(), "r") :
						gzopen((m_fileNames[i] + ".rv").c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				int l;
#pragma omp parallel private(l)
				for (;;) {
					string sequence, name;
#pragma omp critical(seq)
					{
						l = kseq_read(seq);
						if (l >= 0) {
							sequence = string(seq->seq.s, seq->seq.l);
							name = string(seq->name.s, seq->name.l);
						}
					}
					if (l >= 0) {
						loadSeq(*miBFBV, name, sequence, j);
					} else {
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
			}
		}
		if (opt::verbose){
			cerr << "Finishing Pass " << j << " " << omp_get_wtime() - time
					<< "s" << endl;
			time = omp_get_wtime();
		}

		//second pass saturation normalization special
		{
			//record memory before
			size_t memKB = getRSS();
			if (opt::verbose)
				cerr << "Mem usage (kB): " << memKB << endl;

			//use single value Reservoir sampling
			/*
			 * pair<ID,ID> first ID stores the currentID, and ID stores the current observation count
			 * If the second ID exceeds max possible count, the ID is a critical ID
			 * Critical IDs are need for partial hits and always replace existing IDs
			 */
			//TODO since since is known (getPopSaturated()) possible to use sd-bitvector?
			typedef google::sparse_hash_map<size_t, pair<ID,ID>> SatMap;
			SatMap satMap(miBFBV->getPopSaturated());
			const ID criticalCount = m_nameToID.size() + 1;

			if (opt::verbose)
				cerr << "Pass normalize" << endl;
			if (opt::idByFile) {
#pragma omp parallel for schedule(dynamic)
				for (unsigned i = 0; i < m_fileNames.size(); ++i) {
					gzFile fp;
					if (opt::verbose) {
//#pragma omp critical(stderr)
//						cerr << "Opening: "
//								<< (j % 2 ?
//										m_fileNames[i] : m_fileNames[i] + ".rv")
//								<< endl;
					}
					fp = j % 2 ?
							gzopen(m_fileNames[i].c_str(), "r") :
							gzopen((m_fileNames[i] + ".rv").c_str(), "r");
					kseq_t *seq = kseq_init(fp);
					int l;
					for (;;) {
						string sequence, name;
						{
							l = kseq_read(seq);
							if (l >= 0) {
								sequence = string(seq->seq.s, seq->seq.l);
								name = m_fileNames[i].substr(
										m_fileNames[i].find_last_of("/") + 1);
							}
						}
						if (l >= 0) {
							//get positions
							typedef google::dense_hash_set<size_t> SatSet;
							SatSet satVal;
							satVal.set_empty_key(miBFBV->size());
							SatSet critVal;
							critVal.set_empty_key(miBFBV->size());
							ID id = m_nameToID[name];
							recordSaturation(*miBFBV, id, sequence, satVal,
									critVal);
#pragma omp critical(satMap)
							{
								//set critial IDs
								for (SatSet::iterator itr = critVal.begin();
										itr != critVal.end(); itr++) {
									satMap[*itr] = pair<ID, ID>(id,
											criticalCount);
								}
								for (SatSet::iterator itr = satVal.begin();
										itr != satVal.end(); itr++) {
									SatMap::iterator tempItr = satMap.find(
											*itr);
									if (tempItr == satMap.end()) {
										satMap[*itr] = pair<ID, ID>(id, 1);
									} else if (tempItr->second.second
											!= criticalCount) {
										tempItr->second.second++;
										//each position^ID combination is unique
										//is this random enough ont is own?
										size_t randomSeed = *itr ^ id;
										ID randomNum = std::hash<ID> { }(
												randomSeed)
												% tempItr->second.second;
										if (randomNum
												== tempItr->second.second - 1) {
											tempItr->second.first = id;
										}
									}
								}
							}
						} else {
							break;
						}
					}
//					if (opt::verbose) {
//#pragma omp critical(stderr)
//						cerr << "Closing: "
//								<< (j % 2 ?
//										m_fileNames[i] : m_fileNames[i] + ".rv")
//								<< endl;
//					}
					kseq_destroy(seq);
					gzclose(fp);
				}
			} else {
				for (unsigned i = 0; i < m_fileNames.size(); ++i) {
					gzFile fp;
//					if (opt::verbose)
//						cerr << "Opening: "
//								<< (j % 2 ?
//										m_fileNames[i] : m_fileNames[i] + ".rv")
//								<< endl;
					fp = j % 2 ?
							gzopen(m_fileNames[i].c_str(), "r") :
							gzopen((m_fileNames[i] + ".rv").c_str(), "r");
					kseq_t *seq = kseq_init(fp);
					int l;
#pragma omp parallel private(l)
					for (;;) {
						string sequence, name;
#pragma omp critical(seq)
						{
							l = kseq_read(seq);
							if (l >= 0) {
								sequence = string(seq->seq.s, seq->seq.l);
								name = string(seq->name.s, seq->name.l);
							}
						}
						if (l >= 0) {
							//get positions
							typedef google::dense_hash_set<size_t> SatSet;
							SatSet satVal;
							satVal.set_empty_key(miBFBV->size());
							SatSet critVal;
							critVal.set_empty_key(miBFBV->size());
							ID id = m_nameToID[name];
							recordSaturation(*miBFBV, id, sequence, satVal,
									critVal);
#pragma omp critical(satMap)
							{
								//set critial IDs
								for (SatSet::iterator itr = critVal.begin();
										itr != critVal.end(); itr++) {
									satMap[*itr] = pair<ID, ID>(id, criticalCount);
								}
								for (SatSet::iterator itr = satVal.begin();
										itr != satVal.end(); itr++) {
									SatMap::iterator tempItr = satMap.find(
											*itr);
									if (tempItr == satMap.end()) {
										satMap[*itr] = pair<ID, ID>(id, 1);
									} else if (tempItr->second.second
											!= criticalCount) {
										tempItr->second.second++;
										//each position^ID combination is unique
										//is this random enough ont is own?
										size_t randomSeed = *itr ^ id;
										ID randomNum = std::hash<ID> { }(
												randomSeed)
												% tempItr->second.second;
										if (randomNum
												== tempItr->second.second - 1) {
											tempItr->second.first = id;
										}
									}
								}
							}
						} else {
							break;
						}
					}
					kseq_destroy(seq);
					gzclose(fp);
				}
			}

			if (opt::verbose){
				cerr << "Finishing temporary saturation " << omp_get_wtime() - time
						<< "s" << endl;
				time = omp_get_wtime();
			}


			//mutate the saturatedID to random set of saturated IDs
#pragma omp parallel
			for (SatMap::iterator itr = satMap.begin(); itr != satMap.end();
					++itr) {
				//if part of critical list do nothing
				if (itr->second.second != criticalCount) {
					miBFBV->setData(itr->first, itr->second.first | MIBloomFilter<ID>::s_mask);
				}
			}
//			ID check = miBFBV->checkValues(m_nameToID.size());
//			if (m_nameToID.size() != check) {
//				cerr << check << " ID found " << (check
//						& MIBloomFilter<ID>::s_antiMask)
//								<< " (after masking), max possible ID is "
//								<< m_nameToID.size() << endl;
//				exit(1);
//			}
			if (opt::verbose){
				cerr << "Finishing normalization stage " << omp_get_wtime() - time
						<< "s" << endl;
				time = omp_get_wtime();
			}

			//record memory after
			if (opt::verbose)
				cerr << "Mem usage of support datastructure (kB): " << (getRSS() - memKB) << endl;
		}

		//finish the rest
		//j is the number of matches for that iteration possible
		for (j++; j <= opt::hashNum; ++j) {
			if (opt::verbose)
				cerr << "Pass " << j << endl;
			if (opt::idByFile) {
#pragma omp parallel for schedule(dynamic)
				for (unsigned i = 0; i < m_fileNames.size(); ++i) {
					gzFile fp;
					if (opt::verbose) {
#pragma omp critical(stderr)
						cerr << "Opening: "
								<< (j % 2 ?
										m_fileNames[i] : m_fileNames[i] + ".rv")
								<< endl;
					}
					fp = j % 2 ?
							gzopen(m_fileNames[i].c_str(), "r") :
							gzopen((m_fileNames[i] + ".rv").c_str(), "r");
					kseq_t *seq = kseq_init(fp);
					int l;
					for (;;) {
						string sequence, name;
						{
							l = kseq_read(seq);
							if (l >= 0) {
								sequence = string(seq->seq.s, seq->seq.l);
								name = m_fileNames[i].substr(
										m_fileNames[i].find_last_of("/") + 1);
							}
						}
						if (l >= 0) {
							loadSeq(*miBFBV, name, sequence, j);
						} else {
							break;
						}
					}
					kseq_destroy(seq);
					gzclose(fp);
//#pragma omp critical(stderr)
//					if (opt::verbose > 0) {
//						cerr << "Saturation: " << miBFBV->getPopSaturated()
//								<< " popNonZero: " << miBFBV->getPopNonZero()
//								<< endl;
//					}
				}
			} else {
				for (unsigned i = 0; i < m_fileNames.size(); ++i) {
					gzFile fp;
					if (opt::verbose)
						cerr << "Opening: "
								<< (j % 2 ?
										m_fileNames[i] : m_fileNames[i] + ".rv")
								<< endl;
					fp = j % 2 ?
							gzopen(m_fileNames[i].c_str(), "r") :
							gzopen((m_fileNames[i] + ".rv").c_str(), "r");
					kseq_t *seq = kseq_init(fp);
					int l;
#pragma omp parallel private(l)
					for (;;) {
						string sequence, name;
#pragma omp critical(seq)
						{
							l = kseq_read(seq);
							if (l >= 0) {
								sequence = string(seq->seq.s, seq->seq.l);
								name = string(seq->name.s, seq->name.l);
							}
						}
						if (l >= 0) {
							loadSeq(*miBFBV, name, sequence, j);
						} else {
							break;
						}
					}
					kseq_destroy(seq);
					gzclose(fp);
				}
			}
		}
		cerr << "Outputting IDs file: " << filePrefix + "_ids.txt" << endl;
		std::ofstream idFile;
		idFile.open((filePrefix + "_ids.txt").c_str());
		writeIDs(idFile);
		idFile.close();
		cerr << "Failed to insert for set of multi-spaced seeds: "
				<< m_failedInsert << endl;
		cerr << "PopCount: " << miBFBV->getPop() << endl;
		cerr << "PopNonZero: " << miBFBV->getPopNonZero() << endl;
		cerr << "PopSaturated: " << miBFBV->getPopSaturated() << endl;
		cerr << "PopCount Ratio: "
				<< double(miBFBV->getPop()) / double(miBFBV->size()) << endl;
		cerr << "Storing filter" << endl;

		//save filter
		miBFBV->store(filePrefix + ".bf");

		if(opt::verbose > 1){
			vector<size_t> counts(m_ids.size(), 0);
//			cerr << counts.size() << endl;
			miBFBV->getIDCounts(counts);
			size_t count = 0;
			for(vector<size_t>::iterator itr = ++counts.begin(); itr != counts.end(); ++itr){
				cout << ++count << "\t" << *itr << endl;
			}
		}

		return miBFBV;
	}

	vector<string> getIDs() const {
		return m_ids;
	}

private:
	unsigned m_kmerSize;
	size_t m_expectedEntries;
	size_t m_totalEntries;
	vector<string> m_fileNames;
	vector<string> m_ids;
	google::dense_hash_map<string, ID> m_nameToID;
	unsigned m_regionSize;
	size_t m_failedInsert;

	inline MIBloomFilter<ID> * generateBV(double fpr,
			const vector<vector<unsigned> > *ssVal = NULL) {

		size_t filterSize = 0;
		filterSize = MIBloomFilter<ID>::calcOptimalSize(m_expectedEntries,
				opt::hashNum, fpr);

		if (opt::verbose > 0)
			cerr << "Bit vector Size: " << filterSize << endl;
		size_t uniqueCounts = 0;

		sdsl::bit_vector bv(filterSize);

		if (opt::verbose > 0)
			cerr << "Populating initial bit vector" << endl;

		//populate sdsl bitvector (bloomFilter)
		if (opt::idByFile) {
#pragma omp parallel for schedule(dynamic)
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Opening " << m_fileNames[i] << endl;
				}
				fp = gzopen(m_fileNames[i].c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				int l;
				size_t colliCounts = 0;
				size_t totalCount = 0;
				for (;;) {
					string sequence;
					{
						l = kseq_read(seq);
						if (l >= 0) {
							sequence = string(seq->seq.s, seq->seq.l);
						}
					}
					if (l >= 0) {
						if (sequence.length() >= m_kmerSize) {
							//k-merize with rolling hash insert into multi index Bloom filter
							if (ssVal != NULL) {
								colliCounts += loadSeq(bv, sequence, *ssVal);
							} else {
								colliCounts += loadSeq(bv, sequence);
							}
							totalCount += sequence.length() - m_kmerSize + 1;
						}
					} else {
						break;
					}
				}
#pragma omp atomic
				uniqueCounts += totalCount - colliCounts;
				kseq_destroy(seq);
				gzclose(fp);
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Finished processing " << m_fileNames[i] << endl;
				}
			}
		} else {
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				int l;
				size_t colliCounts = 0;
				size_t totalCount = 0;
#pragma omp parallel private(l)
				for (;;) {
					string sequence;
#pragma omp critical(seq)
					{
						l = kseq_read(seq);
						if (l >= 0) {
							sequence = string(seq->seq.s, seq->seq.l);
						}
					}
					if (l >= 0) {
						if (sequence.length() >= m_kmerSize) {
							//k-merize with rolling hash insert into multi index Bloom filter
							if (ssVal != NULL) {
#pragma omp atomic
								colliCounts += loadSeq(bv, sequence, *ssVal);
							} else {
#pragma omp atomic
								colliCounts += loadSeq(bv, sequence);
							}
#pragma omp atomic
							totalCount += sequence.length() - m_kmerSize + 1;
						}
					} else {
						break;
					}
				}
				uniqueCounts += totalCount - colliCounts;
				kseq_destroy(seq);
				gzclose(fp);
			}
		}

		if (opt::verbose > 0) {
			cerr << "Approximate number of unique frames in filter: "
					<< uniqueCounts << endl;
		}
		return new MIBloomFilter<ID>(opt::hashNum, m_kmerSize, bv, opt::sseeds);
	}

//helper methods
	inline void loadSeq(MIBloomFilter<ID> &miBF, const string &name,
			const string &seq, unsigned max = 1) {
		ID id = m_nameToID[name];
		if (miBF.getSeedValues().empty()) {
			ntHashIterator itr(seq, opt::hashNum, m_kmerSize);
			insertTillEnd(miBF, max, itr, id);
		} else {
			stHashIterator itr(seq, miBF.getSeedValues(), opt::hashNum, miBF.getKmerSize());
			insertTillEnd(miBF, max, itr, id);
		}
	}

	template<typename T>
	inline void insertTillEnd(MIBloomFilter<ID> &miBF, unsigned max, T &itr,
			ID id) {
		while (itr != itr.end()) {
			//Last iteration check if value was obliterated
			if (!miBF.insert(*itr, id, max)) {
#pragma omp atomic
				++m_failedInsert;
			}
			++itr;
		}
	}

	/*
	 * Returns set of datavector index where saturation occurs for this ID
	 * Critical values ideally should never be mutated
	 */
	inline void recordSaturation(const MIBloomFilter<ID> &miBF, ID id,
			const string &seq, google::dense_hash_set<size_t> &saturatedValues,
			google::dense_hash_set<size_t> &criticalValues) {
		if (miBF.getSeedValues().empty()) {
			ntHashIterator itr(seq, opt::hashNum, m_kmerSize);
			while (itr != itr.end()) {
				//for each set of hash values, check for saturation
				vector<size_t> rankPos = miBF.getRankPos(*itr);
				vector<ID> results = miBF.getData(rankPos);
				bool saturated = true;
				for (unsigned i = 0; i < opt::hashNum; ++i) {
					ID oldVal = results[i];
					if (oldVal < MIBloomFilter<ID>::s_mask) {
						saturated = false;
						break;
					}
				}
				if (saturated){
					//if completely saturated, record hash location into main set
					for (unsigned i = 0; i < opt::hashNum; ++i) {
						saturatedValues.insert(rankPos[i]);
					}
				}
				else{
					for (unsigned i = 0; i < opt::hashNum; ++i) {
						ID oldVal = results[i];
						if (oldVal > MIBloomFilter<ID>::s_mask) {
							oldVal = oldVal & MIBloomFilter<ID>::s_antiMask;
							//if partially saturated and the same ID, record into critical set
							if(oldVal == id){
								criticalValues.insert(rankPos[i]);
							}
						}
					}
				}
				++itr;
			}
		} else {
			stHashIterator itr(seq, miBF.getSeedValues(), opt::hashNum, miBF.getKmerSize());
			while (itr != itr.end()) {
				//for each set of hash values, check for saturation
				vector<size_t> rankPos = miBF.getRankPos(*itr);
				vector<ID> results = miBF.getData(rankPos);
				bool saturated = true;
				for (unsigned i = 0; i < opt::hashNum; ++i) {
					ID oldVal = results[i];
					if (oldVal < MIBloomFilter<ID>::s_mask) {
						saturated = false;
						break;
					}
				}
				if (saturated){
					//if completely saturated, record hash location into main set
					for (unsigned i = 0; i < opt::hashNum; ++i) {
						saturatedValues.insert(rankPos[i]);
					}
				}
				else{
					for (unsigned i = 0; i < opt::hashNum; ++i) {
						ID oldVal = results[i];
						if (oldVal > MIBloomFilter<ID>::s_mask) {
							oldVal = oldVal & MIBloomFilter<ID>::s_antiMask;
							//if partially saturated and the same ID, record into critical set
							if(oldVal == id){
								criticalValues.insert(rankPos[i]);
							}
						}
					}
				}
				++itr;
			}
		}
	}

	/*
	 * Returns count of collisions (counts unique k-mers)
	 */
	inline size_t loadSeq(sdsl::bit_vector &bv, const string& seq,
			const vector<vector<unsigned> > &seedVal =
					vector<vector<unsigned> >(0)) {
		size_t count = 0;
		if (seq.size() < m_kmerSize)
			return count;
		if (seedVal.empty()) {
			/* init rolling hash state and compute hash values for first k-mer */
			for (ntHashIterator itr(seq, opt::hashNum, m_kmerSize);
					itr != itr.end(); ++itr) {
				unsigned colliCount = 0;
				for (size_t i = 0; i < opt::hashNum; ++i) {
					size_t pos = (*itr)[i] % bv.size();
					uint64_t *dataIndex = bv.data() + (pos >> 6);
					uint64_t bitMaskValue = (uint64_t) 1 << (pos & 0x3F);
					colliCount += __sync_fetch_and_or(dataIndex, bitMaskValue)
							>> (pos & 0x3F) & 1;
				}
				if (colliCount == opt::hashNum) {
					++count;
				}
			}
		} else {
			/* init rolling hash state and compute hash values for first k-mer */
			for (stHashIterator itr(seq, seedVal, opt::hashNum, m_kmerSize);
					itr != itr.end(); ++itr) {
				unsigned colliCount = 0;
				for (size_t i = 0; i < seedVal.size(); ++i) {
					size_t pos = (*itr)[i] % bv.size();
					uint64_t *dataIndex = bv.data() + (pos >> 6);
					uint64_t bitMaskValue = (uint64_t) 1 << (pos & 0x3F);
					colliCount += __sync_fetch_and_or(dataIndex, bitMaskValue)
							>> (pos & 0x3F) & 1;
				}
				if (colliCount == seedVal.size()) {
					++count;
				}
			}
		}
		return count;
	}

	inline void writeIDs(std::ofstream &file) const {
		assert(file);
		for (ID i = 1; i < m_ids.size(); ++i) {
			file << i << "\t" << m_ids[i] << "\n";
			assert(file);
		}
	}

	inline void writeReverseFile(const string &fileName,
			const string &newFileName) const {
		gzFile fp;
		ofstream out = ofstream(newFileName.c_str());
		fp = gzopen(fileName.c_str(), "r");
		kseq_t *seq = kseq_init(fp);
		int l;
		for (;;) {
			l = kseq_read(seq);
			if (l >= 0) {

			} else {
				kseq_destroy(seq);
				break;
			}
		}
		gzclose(fp);
	}

	//TODO move these functions to a common util class?

	/*
	 * Get RSS
	 */
	size_t getRSS(){ //Note: this value is in KB!
	    FILE* file = fopen("/proc/self/status", "r");
	    int result = -1;
	    char line[128];

	    while (fgets(line, 128, file) != NULL){
	        if (strncmp(line, "VmRSS:", 6) == 0){
	            result = parseLine(line);
	            break;
	        }
	    }
	    fclose(file);
	    return result;
	}

	int parseLine(char* line){
	    // This assumes that a digit will be found and the line ends in " Kb".
	    int i = strlen(line);
	    const char* p = line;
	    while (*p <'0' || *p > '9') p++;
	    line[i-3] = '\0';
	    i = atoi(p);
	    return i;
	}

	/*
	 * checks if file exists
	 */
	inline bool fexists(const string &filename) {
		ifstream ifile(filename.c_str());
		return ifile.good();
	}
};

#endif /* CHROMIUMMAP_MIBFGEN_HPP_ */
