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
#include "btl_bloomfilter/MIBFConstructSupport.hpp"
#include "btl_bloomfilter/vendor/stHashIterator.hpp"
#include "Common/sntHashIterator.hpp"

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

		if(m_ids.size() >= ID(1 << (sizeof(ID) * 8 - 1))){
			cerr << "Max Size: " << ID(1 << (sizeof(ID) * 8 - 1)) << " Required Size: " << m_ids.size()  << endl;
			cerr << "Recompile with a larger integer for storing IDs." << endl;
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

	template<typename H>
	void generate(const string &filePrefix, double occ) {
		//keep track of time
		double time = omp_get_wtime();

		MIBFConstructSupport<ID, H> miBFCS(m_expectedEntries, m_kmerSize,
				opt::hashNum, occ, opt::sseeds);
		vector<vector<unsigned> > ssVal;
		if (!opt::sseeds.empty()) {
			ssVal =	MIBloomFilter<ID>::parseSeedString(opt::sseeds);
		}
		generateBV(miBFCS, ssVal);

		if (opt::verbose){
			cerr << "Finishing initial Bit vector construction " <<  omp_get_wtime() - time << "s" << endl;
			time = omp_get_wtime();
			cerr << "Populating values of miBF" << endl;
		}
		MIBloomFilter<ID> *miBF = miBFCS.getEmptyMIBF();

		//record memory before
		size_t memKB = getRSS();
		if (opt::verbose)
			cerr << "Mem usage (kB): " << memKB << endl;

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
						H itr = hashIterator<H>(sequence, ssVal);
						miBFCS.insertMIBF(*miBF, itr, m_nameToID[name]);
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
						H itr = hashIterator<H>(sequence, ssVal);
						miBFCS.insertSaturation(*miBF, itr, m_nameToID[name]);
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
						H itr = hashIterator<H>(sequence, ssVal);
						miBFCS.insertMIBF(*miBF, itr, m_nameToID[name]);
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
						H itr = hashIterator<H>(sequence, ssVal);
						miBFCS.insertSaturation(*miBF, itr, m_nameToID[name]);
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
		cerr << "PopCount: " << miBF->getPop() << endl;
		cerr << "PopSaturated: " << miBF->getPopSaturated() << endl;
		cerr << "PopCount Ratio: "
				<< double(miBF->getPop()) / double(miBF->size()) << endl;
		cerr << "Storing filter" << endl;

		//save filter
		miBF->store(filePrefix + ".bf");

		if(opt::verbose > 1){
			vector<size_t> counts(m_ids.size(), 0);
			miBF->getIDCounts(counts);
			size_t count = 0;
			for(vector<size_t>::iterator itr = ++counts.begin(); itr != counts.end(); ++itr){
				cout << ++count << "\t" << *itr << endl;
			}
		}
		delete(miBF);
	}

	vector<string> getIDs() const {
		return m_ids;
	}

private:
	unsigned m_kmerSize;
	size_t m_expectedEntries;
	vector<string> m_fileNames;
	vector<string> m_ids;
	google::dense_hash_map<string, ID> m_nameToID;

	template<typename H>
	void generateBV(MIBFConstructSupport<ID, H> &miBFCS,
			const vector<vector<unsigned>> &ssVal) {
		if (opt::verbose > 0)
			cerr << "Bit vector Size: " << miBFCS.getFilterSize() << endl;

		size_t uniqueCounts = 0;
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
							H itr = hashIterator<H>(sequence, ssVal);
							colliCounts += miBFCS.insertBVColli(itr);
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
				if(opt::verbose){
#pragma omp critical(stderr)
					cerr << "Opening " << m_fileNames[i] << endl;
				}
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
							H itr = hashIterator<H>(sequence, ssVal);
							colliCounts += miBFCS.insertBVColli(itr);
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
	}


	template<typename H>
	H hashIterator(const string &seq,
			const vector<vector<unsigned> > &seedVal) {
		return H(seq, seedVal, opt::hashNum, 1, m_kmerSize);
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
