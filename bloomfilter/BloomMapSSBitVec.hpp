/*
 * BloomMap.hpp
 *
 *  Created on: Dec 17, 2015
 *      Author: gjahesh
 */

#ifndef BLOOMMAPSSBITVECT_HPP_
#define BLOOMMAPSSBITVECT_HPP_

#include <string>
#include <vector>
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <limits>
#include <google/dense_hash_map>
#include <sdsl/bit_vector_il.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include "BloomMapSS.hpp"

using namespace std;
template<typename T>
class BloomMapSS;

static const unsigned BLOCKSIZE = 512;

template<typename T>
class BloomMapSSBitVec {
public:

#pragma pack(1) //to maintain consistent values across platforms
	struct FileHeader {
		char magic[8];
		uint32_t hlen;	//header length (including spaced seeds)
		uint64_t size;
		uint64_t nEntry;
		uint64_t tEntry;
		double dFPR;
		uint32_t nhash;
		uint32_t kmer;
	};

	/*
	 * Constructor using an existing bloomMap
	 */
	BloomMapSSBitVec<T>(BloomMapSS<T> &bloomMap) :
			m_dSize(bloomMap.getPop()), m_dFPR(bloomMap.getDesiredFPR()), m_tEntry(
					bloomMap.getTotalEntries()), m_sseeds(
					bloomMap.getSeedStrings()), m_kmerSize(
					bloomMap.getKmerSize()), m_ssVal(bloomMap.getSeedValues()) {
		m_data = new T[m_dSize]();
		sdsl::bit_vector tmpBV(bloomMap.size());
		size_t idx = 0;
		for (size_t i = 0; i < tmpBV.size(); ++i) {
			if (bloomMap.m_array[i] != 0) {
				tmpBV[i] = true;
				m_data[idx++] = bloomMap.m_array[i];
			}
		}
		m_bv = sdsl::bit_vector_il<BLOCKSIZE>(tmpBV);
		m_rankSupport = sdsl::rank_support_il<1>(&m_bv);
	}

	/*
	 * Constructor using a prebuilt bitvector
	 */
	BloomMapSSBitVec<T>(size_t expectedElemNum, double fpr,
			vector<string> seeds, const sdsl::bit_vector &bv) :
			m_dSize(0), m_dFPR(fpr), m_tEntry(expectedElemNum), m_sseeds(
					seeds), m_kmerSize(m_sseeds[0].size()), m_ssVal(
					parseSeedString(m_sseeds)) {
		for (vector<string>::const_iterator itr = m_sseeds.begin();
				itr != m_sseeds.end(); ++itr) {
			//check if spaced seeds are all the same length
			assert(m_kmerSize == itr->size());
		}
		m_bv = sdsl::bit_vector_il<BLOCKSIZE>(bv);
		m_rankSupport = sdsl::rank_support_il<1>(&m_bv);
		m_dSize = getPop();
		m_data = new T[m_dSize]();
	}

	BloomMapSSBitVec<T>(const string &filterFilePath) {
#pragma omp parallel for
		for (unsigned i = 0; i < 2; ++i) {
			if (i == 0) {
				FILE *file = fopen(filterFilePath.c_str(), "rb");
				if (file == NULL) {
					cerr << "file \"" << filterFilePath << "\" could not be read."
							<< endl;
					exit(1);
				}

				FileHeader header;
				if (fread(&header, sizeof(struct FileHeader), 1, file) == 1) {
					cerr << "Loading header..." << endl;
				} else {
					cerr << "Failed to Load header" << endl;
					exit(1);
				}
				char magic[9];
				strncpy(magic, header.magic, 8);
				magic[8] = '\0';

				cerr << "Loaded header... magic: " << magic << " hlen: "
						<< header.hlen << " size: " << header.size << " nhash: "
						<< header.nhash << " kmer: " << header.kmer << " dFPR: "
						<< header.dFPR << " nEntry: " << header.nEntry
						<< " tEntry: " << header.tEntry << endl;

				assert(strcmp(MAGIC, magic) == 0);

				m_sseeds = vector<string>(header.nhash);

				//load seeds
				for (unsigned i = 0; i < header.nhash; ++i) {
					char temp[header.kmer];

					if (fread(temp, header.kmer, 1, file) != 1) {
						cerr << "Failed to load spaced seed string" << endl;
						exit(1);
					} else {
						cerr << "Spaced Seed " << i << ": "
								<< string(temp, header.kmer) << endl;
					}
					m_sseeds[i] = string(temp, header.kmer);
				}
				m_dFPR = header.dFPR;
				m_nEntry = header.nEntry;
				m_tEntry = header.tEntry;
				m_kmerSize = header.kmer;
				m_dSize = header.size;
				m_data = new T[m_dSize]();

				m_ssVal = parseSeedString(m_sseeds);

				cerr << "Loading data vector" << endl;

				long int lCurPos = ftell(file);
				fseek(file, 0, 2);
				size_t fileSize = ftell(file) - header.hlen;
				fseek(file, lCurPos, 0);
				if (fileSize != m_dSize * sizeof(T)) {
					cerr << "Error: " << filterFilePath
							<< " does not match size given by its header. Size: "
							<< fileSize << " vs " << m_dSize * sizeof(T)
							<< " bytes." << endl;
					exit(1);
				}

				size_t countRead = fread(m_data, fileSize, 1, file);
				if (countRead != 1 && fclose(file) != 0) {
					cerr << "file \"" << filterFilePath
							<< "\" could not be read." << endl;
					exit(1);
				}
			}
			else{
				string bvFilename = filterFilePath + ".sdsl";
				cerr << "Loading sdsl interleaved bit vector from: "
						<< bvFilename << endl;
				load_from_file(m_bv, bvFilename);
				m_rankSupport = sdsl::rank_support_il<1>(&m_bv);
			}
		}

		cerr << "FPR based on PopCount: " << getFPR() << endl;
	}

	/*
	 * Stores the filter as a binary file to the path specified
	 * Stores uncompressed because the random data tends to
	 * compress poorly anyway
	 */
	void store(string const &filterFilePath) const {

#pragma omp parallel for
		for (unsigned i = 0; i < 2; ++i) {
			if (i == 0) {
				ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

				assert(myFile);
				writeHeader(myFile);

				cerr << "Storing filter. Filter is " << m_dSize * sizeof(T)
						<< "bytes." << endl;

				//write out each block
				myFile.write(reinterpret_cast<char*>(m_data),
						m_dSize * sizeof(T));

				myFile.close();
				assert(myFile);

				FILE *file = fopen(filterFilePath.c_str(), "rb");
				if (file == NULL) {
					cerr << "file \"" << filterFilePath
							<< "\" could not be read." << endl;
					exit(1);
				}
			} else {
				string bvFilename = filterFilePath + ".sdsl";

				cerr << "Storing sdsl interleaved bit vector to: " << bvFilename
						<< endl;
				store_to_file(m_bv, bvFilename);
				cerr << "Number of bit vector buckets is " << m_bv.size()
						<< endl;
				cerr << "Uncompressed bit vector size is "
						<< (m_bv.size() + m_bv.size() * 64 / BLOCKSIZE) / 8
						<< "bytes" << endl;
			}
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 * ONLY REPLACE VALUE if it is larger than current value (for thread safety)
	 */
	void insert(std::vector<size_t> const &hashes, T value) {
		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < hashes.size(); ++i) {
			size_t pos = m_rankSupport(hashes.at(i) % m_bv.size());
			setIfGreater(&m_data[pos], value);
		}
	}

	/*
	 * Mutates value vector to contain values in bloom map
	 * Assumes vector is the size of hashes
	 */
	void query(std::vector<size_t> const &hashes, vector<T> &values) const {
		for (unsigned i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_bv.size();
			if (m_bv[pos] == 0) {
				values[i] = 0;
				continue;
			}
			values[i] = m_data[m_rankSupport(pos)];
		}
	}

	/*
	 * Returns unambiguous hit to object
	 * Returns best hit on ambiguous collision
	 * Returns numeric_limits<T>::max() on completely ambiguous collision
	 * Returns 0 on if missing element
	 */
	T atBest(std::vector<size_t> const &hashes, unsigned missMin) const {
		google::dense_hash_map<T, unsigned> tmpHash;
		tmpHash.set_empty_key(0);
		unsigned maxCount = 0;
		unsigned miss = 0;
		T value = 0;

		for (unsigned i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_bv.size();
			if(m_bv[pos] == 0){
				++miss;
				continue;
			}
			size_t rankPos = m_rankSupport(pos);
			if (tmpHash.find(m_data[rankPos]) != tmpHash.end()) {
				++tmpHash[m_data[rankPos]];
			} else {
				tmpHash[m_data[rankPos]] = 1;
			}
			if(maxCount == tmpHash[m_data[rankPos]]) {
				value = m_data[rankPos] < value ? m_data[rankPos] : value;
			}
			else if ( maxCount < tmpHash[m_data[rankPos]] ){
				value = m_data[rankPos];
				maxCount = tmpHash[m_data[rankPos]];
			}
		}
		if (missMin < miss) {
			return 0;
		} else {
			return value;
		}
	}

	const vector< vector <unsigned > > &getSeedValues() const {
		return m_ssVal;
	}

	unsigned getKmerSize(){
		return m_kmerSize;
	}

	/*
	 * Return FPR based on popcount
	 */
	double getFPR() const {
		return pow(double(getPop())/double(m_bv.size()), double(m_ssVal.size()));
	}

	/*
	 * Return FPR based on popcount and minimum number of matches for a hit
	 */
	double getFPR(unsigned minNum) const {
		assert(minNum <= m_ssVal.size());
		double cumulativeProb= 0;
		double popCount = getPop();
		double p = popCount/double(m_bv.size());
		for (unsigned i = 0; i < m_ssVal.size() - minNum; ++i) {
			cumulativeProb += double(nChoosek(m_ssVal.size(), i)) * pow(p, i)
					* pow(1.0 - p, (m_ssVal.size() - i));
		}
		return pow(double(getPop())/double(m_bv.size()), double(minNum));
	}

	/*
	 * Return FPR based on number of inserted elements
	 */
	double getFPR_numEle() const {
		assert(m_nEntry > 0);
		return calcFPR_numInserted(m_nEntry);
	}

	size_t getPop() const {
		size_t count = 0;
		size_t index = m_bv.size();
		while(count == 0){
			count = m_rankSupport(--index);
		}
		return count;
	}

	size_t getUniqueEntries() const {
		return m_nEntry;
	}

	void setUnique(size_t count){
		m_nEntry = count;
	}

	~BloomMapSSBitVec() {
		delete[] m_data;
	}

private:

	/*
	 * Helper function for header storage
	 */
	void writeHeader(ofstream &out) const {
		FileHeader header;
		strncpy(header.magic, MAGIC, 8);
		char magic[9];
		strncpy(magic, header.magic, 8);
		magic[8] = '\0';

		header.hlen = sizeof(struct FileHeader) + m_kmerSize * m_sseeds.size();
		header.kmer = m_kmerSize;
		header.size = m_dSize;
		header.nhash = m_sseeds.size();
		header.dFPR = m_dFPR;
		header.nEntry = m_nEntry;
		header.tEntry = m_tEntry;

		cerr << "Writing header... magic: " << magic << " hlen: " << header.hlen
				<< " size: " << header.size << " dFPR: " << header.dFPR
				<< " nEntry: " << header.nEntry << " tEntry: " << header.tEntry
				<< endl;

		out.write(reinterpret_cast<char*>(&header), sizeof(struct FileHeader));

		for (vector<string>::const_iterator itr = m_sseeds.begin();
				itr != m_sseeds.end(); ++itr) {
			out.write(itr->c_str(), m_kmerSize);
		}
	}

	/*
	 * Parses spaced seed string (string consisting of 1s and 0s) to vector
	 */
	inline vector< vector<unsigned> > parseSeedString(const vector<string> &spacedSeeds) {
		SeedVal seeds(spacedSeeds.size());
		for(unsigned i = 0; i < spacedSeeds.size(); ++i){
			const string ss = spacedSeeds.at(i);
			for(unsigned j = 0; j < ss.size(); ++j){
				if(ss.at(j) == '0'){
					seeds[i].push_back(j);
				}
			}
		}
		return seeds;
	}

	/*
	 * Calculates the optimal number of hash function to use
	 * Calculation assumes optimal ratio of bytes per entry given a fpr
	 */
	static unsigned calcOptiHashNum(double fpr) {
		return unsigned(-log(fpr) / log(2));
	}

	/*
	 * Calculate FPR based on hash functions, size and number of entries
	 * see http://en.wikipedia.org/wiki/Bloom_filter
	 */
	double calcFPR_numInserted(size_t numEntr) const {
		return pow( 1.0 - pow(1.0 - 1.0 / double(m_bv.size()),
					double(numEntr) * double(m_ssVal.size())), double(m_ssVal.size()));
	}

	/*
	 * Calculates the optimal FPR to use based on hash functions
	 */
	double calcFPR_hashNum(int hashFunctNum) const {
		return pow(2.0, -hashFunctNum);
	}

	/*
	 * Lock free cas value setting for larger element
	 */
	void setIfGreater(T *val, T newVal) {
//		assert(newVal);
//		__sync_or_and_fetch(val, 1);
		T oldValue;
		do {
			oldValue = *val;
			if (oldValue >= newVal)
				break;
		} while (!__sync_bool_compare_and_swap(val, oldValue, newVal));
	}

	unsigned nChoosek( unsigned n, unsigned k ) const
	{
	    if (k > n) return 0;
	    if (k * 2 > n) k = n-k;
	    if (k == 0) return 1;

	    int result = n;
	    for( unsigned i = 2; i <= k; ++i ) {
	        result *= (n-i+1);
	        result /= i;
	    }
	    return result;
	}

	//size of bitvector
	size_t m_dSize;

	sdsl::bit_vector_il<BLOCKSIZE> m_bv;
	T* m_data;
	sdsl::rank_support_il<1> m_rankSupport;

	double m_dFPR;
	uint64_t m_nEntry;
	uint64_t m_tEntry;
	vector<string> m_sseeds;
	unsigned m_kmerSize;

	typedef vector<vector<unsigned> > SeedVal;
	SeedVal m_ssVal;
	const char* MAGIC = "BLOOMMBV";
};

#endif /* BLOOMMAPSSBITVECT_HPP_ */
