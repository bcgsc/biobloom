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
#include "btl_bloomfilter/vendor/stHashIterator.hpp"
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"

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
			m_kmerSize(kmerSize), m_expectedEntries(numElements), m_fileNames(filenames) {
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

		//record memory before
		size_t memKB = getRSS();
		if (opt::verbose)
			cerr << "Mem usage (kB): " << memKB << endl;

		//use single value Reservoir sampling
		/*
		 * pair<ID,ID> first ID stores the currentID, and ID stores the current observation count
		 * If the second ID exceeds max possible count, the ID is a critical ID
		 * Critical IDs are needed for partial hits and always replace existing IDs
		 */
		//TODO since size is known (getPopSaturated()) possible to use sd-bitvector?
		vector<ID> counts(miBFBV->getPop(), 0);
		if (opt::verbose)
			cerr << "Pass normalize" << endl;
		if (opt::idByFile) {
#pragma omp parallel for schedule(dynamic)
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Opening " << m_fileNames[i] << endl;
				}
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
						typedef google::dense_hash_set<uint64_t> SatSet;
						SatSet satVal;
						satVal.set_empty_key(miBFBV->size());
						ID id = m_nameToID[name];
						if (opt::verbose) {
#pragma omp critical(stderr)
							cerr << "Recording Hash Positions "
									<< m_fileNames[i] << endl;
						}

						recordPos(*miBFBV, sequence, satVal);

						if (opt::verbose) {
#pragma omp critical(stderr)
							cerr << "Reservoir sampling positions "
									<< m_fileNames[i] << endl;
						}
						for (SatSet::iterator itr = satVal.begin();
								itr != satVal.end(); itr++) {
							uint64_t randomSeed = *itr ^ id;
							uint64_t rank = miBFBV->getRankPos(*itr);
							ID count = __sync_add_and_fetch(&counts[rank], 1);
							ID randomNum = std::hash<ID> { }(randomSeed)
									% count;
							if (randomNum == count - 1) {
								miBFBV->setData(rank, id);
							}
						}
					} else {
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Finished processing " << m_fileNames[i] << endl;
				}
			}
			//apply saturation
			if(opt::verbose){
				cerr << "Applying saturation" << endl;
			}
#pragma omp parallel for schedule(dynamic)
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Opening " << m_fileNames[i] << endl;
				}
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
						typedef google::dense_hash_set<uint64_t> SatSet;
						SatSet satVal;
						satVal.set_empty_key(miBFBV->size());
						ID id = m_nameToID[name];
						setSatIfMissing(*miBFBV, id, sequence, counts);
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
				fp = gzopen(m_fileNames[i].c_str(), "r");
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
						typedef google::dense_hash_set<uint64_t> SatSet;
						SatSet satVal;
						satVal.set_empty_key(miBFBV->size());
						ID id = m_nameToID[name];

						if (opt::verbose) {
#pragma omp critical(stderr)
							cerr << "Recording Hash Positions "
									<< m_fileNames[i] << endl;
						}

						recordPos(*miBFBV, sequence, satVal);

						if (opt::verbose) {
#pragma omp critical(stderr)
							cerr << "Reservoir sampling positions "
									<< m_fileNames[i] << endl;
						}

						for (SatSet::iterator itr = satVal.begin();
								itr != satVal.end(); itr++) {
							uint64_t randomSeed = *itr ^ id;
							uint64_t rank = miBFBV->getRankPos(*itr);
							ID count = __sync_add_and_fetch(&counts[rank], 1);
							ID randomNum = std::hash<ID> { }(randomSeed)
									% count;
							if (randomNum == count - 1) {
								miBFBV->setData(rank, id);
							}
						}
					} else {
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
			}
			//apply saturation
			if(opt::verbose){
				cerr << "Applying saturation" << endl;
			}
			//another pass through references
			//if target frame does not have a single representative, mark frame as saturated
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
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
						typedef google::dense_hash_set<uint64_t> SatSet;
						SatSet satVal;
						satVal.set_empty_key(miBFBV->size());
						ID id = m_nameToID[name];
						setSatIfMissing(*miBFBV, id, sequence, counts);
					} else {
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
			}
		}

		cerr << "Outputting IDs file: " << filePrefix + "_ids.txt" << endl;
		std::ofstream idFile;
		idFile.open((filePrefix + "_ids.txt").c_str());
		writeIDs(idFile);
		idFile.close();
//		cerr << "Failed to insert for set of multi-spaced seeds: "
//				<< m_failedInsert << endl;
		cerr << "PopCount: " << miBFBV->getPop() << endl;
//		cerr << "PopNonZero: " << miBFBV->getPopNonZero() << endl;
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
//	size_t m_totalEntries;
	vector<string> m_fileNames;
	vector<string> m_ids;
	google::dense_hash_map<string, ID> m_nameToID;
//	unsigned m_regionSize;
//	size_t m_failedInsert;

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
//	inline void loadSeq(MIBloomFilter<ID> &miBF, const string &name,
//			const string &seq, unsigned max = 1) {
//		ID id = m_nameToID[name];
//		if (miBF.getSeedValues().empty()) {
//			ntHashIterator itr(seq, opt::hashNum, m_kmerSize);
//			insertTillEnd(miBF, max, itr, id);
//		} else {
//			stHashIterator itr(seq, miBF.getSeedValues(), opt::hashNum, miBF.getKmerSize());
//			insertTillEnd(miBF, max, itr, id);
//		}
//	}

//	template<typename T>
//	inline void insertTillEnd(MIBloomFilter<ID> &miBF, unsigned max, T &itr,
//			ID id) {
//		while (itr != itr.end()) {
//			//Last iteration check if value was obliterated
//			if (!miBF.insert(*itr, id, max)) {
//#pragma omp atomic
//				++m_failedInsert;
//			}
//			++itr;
//		}
//	}

	/*
	 * Returns set of datavector index where saturation occurs for this ID
	 * Critical values ideally should never be mutated
	 */
	void recordPos(const MIBloomFilter<ID> &miBF, const string &seq,
			google::dense_hash_set<uint64_t> &values) {
		if (miBF.getSeedValues().empty()) {
			ntHashIterator itr(seq, opt::hashNum, m_kmerSize);
			while (itr != itr.end()) {
				for (unsigned i = 0; i < opt::hashNum; ++i) {
					values.insert((*itr)[i]);
				}
				++itr;
			}
		} else {
			stHashIterator itr(seq, miBF.getSeedValues(), opt::hashNum, miBF.getKmerSize());
			while (itr != itr.end()) {
				for (unsigned i = 0; i < opt::hashNum; ++i) {
					values.insert((*itr)[i]);
				}
				++itr;
			}
		}
	}

	inline void setSatIfMissing(MIBloomFilter<ID> &miBF, ID id,
			const string &seq, vector<ID> &counts) {
		if (miBF.getSeedValues().empty()) {
			ntHashIterator itr(seq, opt::hashNum, m_kmerSize);
			while (itr != itr.end()) {
				//for each set of hash values, check for saturation
				vector<uint64_t> rankPos = miBF.getRankPos(*itr);
				vector<ID> results = miBF.getData(rankPos);
				vector<ID> replacementIDs(opt::hashNum);
				bool valueFound = false;
				vector<ID> seenSet(opt::hashNum);
				for (unsigned i = 0; i < opt::hashNum; ++i) {
					ID currentResult = results[i] & MIBloomFilter<ID>::s_antiMask;
					if (currentResult == id) {
						valueFound = true;
						break;
					}
					if (find(seenSet.begin(), seenSet.end(), currentResult)
							== seenSet.end()) {
						seenSet.push_back(currentResult);
					}
					else{
						replacementIDs.push_back(currentResult);
					}
				}
				if (!valueFound) {
					uint64_t replacementPos = counts.size();
					ID minCount = numeric_limits<ID>::min();
					for (unsigned i = 0; i < opt::hashNum; ++i) {
						ID currentResult = results[i] & MIBloomFilter<ID>::s_antiMask;
						if (find(replacementIDs.begin(), replacementIDs.end(),
								currentResult) != replacementIDs.end()) {
							if(minCount < counts[rankPos[i]]){
								minCount = counts[rankPos[i]];
								replacementPos = rankPos[i];
							}
						}
					}
					//mutate if possible
					if (replacementPos != counts.size()) {
						miBF.setData(replacementPos, id);
#pragma omp atomic update
						++counts[replacementPos];
					}
					else{
						miBF.saturate(*itr);
					}
				}
				++itr;
			}
		} else {
			stHashIterator itr(seq, miBF.getSeedValues(), opt::hashNum, miBF.getKmerSize());
			while (itr != itr.end()) {
				//for each set of hash values, check for saturation
				vector<uint64_t> rankPos = miBF.getRankPos(*itr);
				vector<ID> results = miBF.getData(rankPos);
				vector<ID> replacementIDs(opt::hashNum);
				bool valueFound = false;
				vector<ID> seenSet(opt::hashNum);
				for (unsigned i = 0; i < opt::hashNum; ++i) {
					ID currentResult = results[i] & MIBloomFilter<ID>::s_antiMask;
					if (currentResult == id) {
						valueFound = true;
						break;
					}
					if (find(seenSet.begin(), seenSet.end(), currentResult)
							== seenSet.end()) {
						seenSet.push_back(currentResult);
					}
					else{
						replacementIDs.push_back(currentResult);
					}
				}
				if (!valueFound) {
					uint64_t replacementPos = counts.size();
					ID minCount = numeric_limits<ID>::min();
					for (unsigned i = 0; i < opt::hashNum; ++i) {
						ID currentResult = results[i] & MIBloomFilter<ID>::s_antiMask;
						if (find(replacementIDs.begin(), replacementIDs.end(),
								currentResult) != replacementIDs.end()) {
							if(minCount < counts[rankPos[i]]){
								minCount = counts[rankPos[i]];
								replacementPos = rankPos[i];
							}
						}
					}
					//mutate if possible
					if (replacementPos != counts.size()) {
						miBF.setData(replacementPos, id);
#pragma omp atomic update
						++counts[replacementPos];
					}
					else{
						miBF.saturate(*itr);
					}
				}
				++itr;
			}
		}
	}

	/*
	 * Returns set of datavector index where saturation occurs for this ID
	 * Critical values ideally should never be mutated
	 */
	inline void recordSaturation(const MIBloomFilter<ID> &miBF, ID id,
			const string &seq, google::dense_hash_set<uint64_t> &saturatedValues,
			google::dense_hash_set<uint64_t> &criticalValues) {
		if (miBF.getSeedValues().empty()) {
			ntHashIterator itr(seq, opt::hashNum, m_kmerSize);
			while (itr != itr.end()) {
				//for each set of hash values, check for saturation
				vector<uint64_t> rankPos = miBF.getRankPos(*itr);
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
				vector<uint64_t> rankPos = miBF.getRankPos(*itr);
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
				for (unsigned i = 0; i < opt::hashNum; ++i) {
					uint64_t pos = (*itr)[i] % bv.size();
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
				for (unsigned i = 0; i < seedVal.size(); ++i) {
					uint64_t pos = (*itr)[i] % bv.size();
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
