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
#include <sdsl/rank_support.hpp>
#include <BloomMapSS.hpp>

using namespace std;

static const char* MAGIC = "BLOOMMSS";

template<typename T>
class BloomMapSSBitVec {
friend class BloomMapSS;

public:

	struct FileHeader {
		char magic[8];
		uint32_t hlen;	//header length (including spaced seeds)
		uint64_t size;
		uint32_t nhash;
		uint32_t kmer;
		double dFPR;
		uint64_t nEntry;
		uint64_t tEntry;
	};

	BloomMapSSBitVec<T>(BloomMapSS<T> &bloomMap) :
			m_size(bloomMap.m_size), m_bv(sdsl::bit_vector_il(m_size)), m_data(), m_dFPR(
					bloomMap.m_dFPR), m_nEntry(bloomMap.m_nEntry), m_tEntry(
					bloomMap.m_tEntry), m_sseeds(bloomMap.m_seeds), m_kmerSize(
					bloomMap.m_kmerSize), m_ssVal(bloomMap.m_sseeds) {
		m_data = new T[bloomMap.getPop()]();
		size_t idx = 0;
		for(size_t i = 0; i< m_size; ++i){
			if(bloomMap.m_array[i] != 0){
				m_bv[i] = true;
				bloomMap[idx++] = bloomMap.m_array[i];
			}
		}
		m_rankSupport(&m_bv);
	}

	BloomMapSSBitVec<T>(const string &filterFilePath) {
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
			cerr << "Failed to header" << endl;
			exit(1);
		}
		char magic[9];
		strncpy(magic, header.magic, 8);
		magic[8] = '\0';

		cerr << "Loaded header... magic: " << magic << " hlen: " << header.hlen
				<< " size: " << header.size << " nhash: " << header.nhash
				<< " kmer: " << header.kmer << " dFPR: " << header.dFPR
				<< " nEntry: " << header.nEntry << " tEntry: " << header.tEntry
				<< endl;

		assert(strcmp(MAGIC, magic) == 0);

		m_sseeds = vector<string>(header.nhash);

		//load seeds
		for (unsigned i = 0; i < header.nhash; ++i) {
			char temp[header.kmer];

			if (fread(temp, header.kmer, 1, file) != 1) {
				cerr << "Failed to load spaced seed string" << endl;
				exit(1);
			} else {
				cerr << "Spaced Seed " << i << ": " << string(temp, header.kmer)
						<< endl;
			}
			m_sseeds[i] = string(temp, header.kmer);
		}
		m_dFPR = header.dFPR;
		m_nEntry = header.nEntry;
		m_tEntry = header.tEntry;
		m_kmerSize = header.kmer;
		m_size = header.size;
		m_data = new T[m_size]();

		m_ssVal = parseSeedString(m_sseeds);

		long int lCurPos = ftell(file);
		fseek(file, 0, 2);
		size_t fileSize = ftell(file) - header.hlen;
		fseek(file, lCurPos, 0);
		if (fileSize != m_size * sizeof(T)) {
			cerr << "Error: " << filterFilePath
					<< " does not match size given by its header. Size: "
					<< fileSize << " vs " << m_size * sizeof(T) << " bytes."
					<< endl;
			exit(1);
		}

		size_t countRead = fread(m_data, fileSize, 1, file);
		if (countRead != 1 && fclose(file) != 0) {
			cerr << "file \"" << filterFilePath << "\" could not be read."
					<< endl;
			exit(1);
		}

		string bvFilename = filterFilePath + ".sdsl";
		cerr << "Loading sdsl interleaved bit vector from: " << bvFilename << endl;
		load_from_file(m_bv, bvFilename);

		cerr << "FPR based on PopCount: " << getFPR() << endl;
		cerr << "FPR based on Unique Elements: " << getFPR_numEle() << endl;
	}

	/*
	 * Stores the filter as a binary file to the path specified
	 * Stores uncompressed because the random data tends to
	 * compress poorly anyway
	 */
	void store(string const &filterFilePath) const {
		ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

		cerr << "Storing filter. Filter is " << m_size * sizeof(T) << "bytes."
				<< endl;

		assert(myFile);
		writeHeader(myFile);

		//write out each block
		myFile.write(reinterpret_cast<char*>(m_data), m_size * sizeof(T));

		myFile.close();
		assert(myFile);

		string bvFilename = filterFilePath + ".sdsl";

		cerr << "Storing sdsl interleaved bit vector to: " << bvFilename << endl;
		store_to_file(m_bv, bvFilename);
	}

	/*
	 * Mutates value vector to contain values in bloom map
	 * Assumes vector is the size of hashes
	 */
	void query(std::vector<size_t> const &hashes, vector<T> &values) const {
		for (unsigned i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_size;
			values[i] = m_data[m_rankSupport(pos)];
		}
	}

	/*
	 * Returns unambiguous hit to object
	 * Returns best hit on ambiguous collision
	 * Returns numeric_limits<T>::max() on completely ambiguous collision
	 * Returns 0 on if missing element
	 */
	//TODO investigate more efficient ways to do this
	T atBest(std::vector<size_t> const &hashes) const {
		google::dense_hash_map<T, unsigned> tmpHash;
		tmpHash.set_empty_key(0);
		unsigned maxCount = 0;
		T value = 0;

		for (unsigned i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_size;
			if(m_data[m_rankSupport(pos)] == 0){
				return 0;
			}
			if (tmpHash.find(m_data[m_rankSupport(pos)]) != tmpHash.end()) {
				++tmpHash[m_data[m_rankSupport(pos)]];
			} else {
				tmpHash[m_data[m_rankSupport(pos)]] = 1;
			}
			if(maxCount == tmpHash[m_data[m_rankSupport(pos)]]) {
				value = numeric_limits<T>::max();
			}
			else if ( maxCount < tmpHash[m_data[m_rankSupport(pos)]] ){
				value = m_data[m_rankSupport(pos)];
				maxCount = tmpHash[m_data[m_rankSupport(pos)]];
			}
		}
		return value;
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
			size_t pos = hashes.at(i) % m_size;
			if(m_data[m_rankSupport(pos)] == 0){
				++miss;
				continue;
			}
			if (tmpHash.find(m_data[m_rankSupport(pos)]) != tmpHash.end()) {
				++tmpHash[m_data[m_rankSupport(pos)]];
			} else {
				tmpHash[m_data[m_rankSupport(pos)]] = 1;
			}
			if(maxCount == tmpHash[m_data[m_rankSupport(pos)]]) {
				value = numeric_limits<T>::max();
			}
			else if ( maxCount < tmpHash[m_data[m_rankSupport(pos)]] ){
				value = m_data[m_rankSupport(pos)];
				maxCount = tmpHash[m_data[m_rankSupport(pos)]];
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
		return pow(double(getPop())/double(m_size), double(m_ssVal.size()));
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
		size_t index = m_size;
		while(count == 0){
			count = m_rankSupport(--index);
		}
		return count;
	}

	size_t getSize() const {
		return m_size;
	}

	void setUnique(size_t count){
		m_nEntry = count;
	}

	~BloomMapSSBitVec() {
		delete[] m_data;
	}

private:

	/*
	 * Helper functions for filter storage
	 */

	void writeHeader(ofstream &out) const {
		FileHeader header;
		strncpy(header.magic, MAGIC, 8);
		char magic[9];
		strncpy(magic, header.magic, 8);
		magic[8] = '\0';

		header.hlen = sizeof(struct FileHeader) + m_kmerSize * m_sseeds.size();
		header.kmer = m_kmerSize;
		header.size = m_size;
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
	 * Only returns multiples of 64 for filter building purposes
	 * Is an estimated size using approximations of FPR formula
	 * given the number of hash functions
	 */
	size_t calcOptimalSize(size_t entries, double fpr) const {
		size_t non64ApproxVal = size_t(
				-double(entries) * double(m_ssVal.size())
						/ log(1.0 - pow(fpr, double(1 / double(m_ssVal.size())))));

		return non64ApproxVal + (64 - non64ApproxVal % 64);
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
		return pow( 1.0 - pow(1.0 - 1.0 / double(m_size),
					double(numEntr) * double(m_ssVal.size())), double(m_ssVal.size()));
	}

	/*
	 * Calculates the optimal FPR to use based on hash functions
	 */
	double calcFPR_hashNum(int hashFunctNum) const {
		return pow(2.0, -hashFunctNum);
	}

	size_t m_size;
	sdsl::bit_vector_il m_bv;
	T* m_data;
	sdsl::rank_support_il<1> m_rankSupport;

	double m_dFPR;
	uint64_t m_nEntry;
	uint64_t m_tEntry;
	vector<string> m_sseeds;
	unsigned m_kmerSize;

	typedef vector< vector<unsigned> > SeedVal;
	SeedVal m_ssVal;
};

#endif /* BLOOMMAPSSBITVECT_HPP_ */
