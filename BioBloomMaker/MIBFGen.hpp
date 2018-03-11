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

#include "bloomfilter/MIBloomFilter.hpp"
#include "bloomfilter/RollingHashIterator.h"
#include "btl_bloomfilter/BloomFilter.hpp"
#include "btl_bloomfilter/ntHashIterator.hpp"

#include "Common/Options.h"

#include <tuple>
#include <google/dense_hash_map>
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
				m_ids.push_back(m_fileNames[i].substr(
						m_fileNames[i].find_last_of("/") + 1));
				m_nameToID[m_ids.back()] = m_ids.size() - 1;
				int l;
				for (;;) {
					l = kseq_read(seq);
					if (l >= 0) {
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
		cerr << "Populating values of miBF" << endl;
		//j is the number of matches for that iteration possible
		for (unsigned j = 1; j <= opt::hashNum; ++j) {
			if (opt::verbose)
				cerr << "Pass " << j << endl;
			if (opt::idByFile) {
#pragma omp parallel for
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
					if (opt::verbose > 0) {
						cerr << "Saturation: " << miBFBV->getPopSaturated()
								<< " popNonZero: " << miBFBV->getPopNonZero()
								<< endl;
					}
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
					if (opt::verbose > 0) {
						cerr << "Saturation: " << miBFBV->getPopSaturated()
								<< " popNonZero: " << miBFBV->getPopNonZero()
								<< endl;
					}
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

		if (opt::verbose > 0) {
			cerr << "Approximate number of unique entries in filter: "
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
#pragma omp atomic update
			m_failedInsert += miBF.insert(itr, max, id);
		} else {
			RollingHashIterator itr(seq, m_kmerSize, miBF.getSeedValues());
#pragma omp atomic update
			m_failedInsert += miBF.insert(itr, max, id);
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
			for (RollingHashIterator itr(seq, m_kmerSize, seedVal);
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

	/*
	 * checks if file exists
	 */
	inline bool fexists(const string &filename) {
		ifstream ifile(filename.c_str());
		return ifile.good();
	}
};

#endif /* CHROMIUMMAP_MIBFGEN_HPP_ */
