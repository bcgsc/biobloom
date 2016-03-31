/*
 * BloomMap.hpp
 *
 *  Created on: Dec 17, 2015
 *      Author: gjahesh
 */

#ifndef BLOOMMAP_HPP_
#define BLOOMMAP_HPP_

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
//#include "rolling.h"
#include <limits>
#include <google/dense_hash_map>

using namespace std;

static const char* MAGIC = "BlOOMMSS";

template<typename T>
class BloomMap {

public:

	struct FileHeader {
		char magic[8];
		uint32_t hlen;
		uint64_t size;
		double dFPR;
		uint64_t nEntry;
		uint64_t tEntry;
	};

	BloomMap<T>(size_t filterSize, vector<string> seeds) :
			m_size(filterSize), m_dFPR(0), m_nEntry(0), m_tEntry(0) {
		m_ssVal = parseSeedString(seeds);
		m_array = new T[m_size]();
	}

	/* De novo filter constructor.
	 * Allocates a filter size based on the number of expected elements and FPR
	 *
	 * If hashNum is set to 0, an optimal value is computed based on the FPR
	 */
	BloomMap<T>(size_t expectedElemNum, double fpr, vector<string> seeds) :
			m_size(0), m_dFPR(fpr), m_nEntry(
					0), m_tEntry(0) {
		m_ssVal = parseSeedString(seeds);
		if (m_size == 0) {
			m_size = calcOptimalSize(expectedElemNum, m_dFPR);
		}
		m_array = new T[m_size]();
	}

	~BloomMap() {
		delete[] m_array;
	}

	BloomMap<T>(const string &filterFilePath) {
		FILE *file = fopen(filterFilePath.c_str(), "rb");
		if (file == NULL) {
			cerr << "file \"" << filterFilePath << "\" could not be read."
					<< endl;
			exit(1);
		}

		loadHeader(file);

		long int lCurPos = ftell(file);
		fseek(file, 0, 2);
		size_t fileSize = ftell(file) - sizeof(struct FileHeader);
		fseek(file, lCurPos, 0);
		if (fileSize != m_size * sizeof(T)) {
			cerr << "Error: " << filterFilePath
					<< " does not match size given by its header. Size: "
					<< fileSize << " vs " << m_size * sizeof(T) << " bytes."
					<< endl;
			exit(1);
		}

		size_t countRead = fread(m_array, fileSize, 1, file);
		if (countRead != 1 && fclose(file) != 0) {
			cerr << "file \"" << filterFilePath << "\" could not be read."
					<< endl;
			exit(1);
		}
	}

	void loadHeader(FILE *file) {

		FileHeader header;
		if (fread(&header, sizeof(struct FileHeader), 1, file) == 1) {
			cerr << "Loading header..." << endl;
		} else {
			cerr << "Failed to header" << endl;
		}
		char magic[9];
		strncpy(magic, header.magic, 8);
		magic[8] = '\0';

        cerr << "Loaded header... magic: " <<
            magic << " hlen: " <<
            header.hlen << " size: " <<
            header.size << " nhash: " <<
            header.dFPR << " aFPR: " <<
            header.nEntry << " tEntry: " <<
            header.tEntry << endl;

        assert(MAGIC == magic);

		m_size = header.size;
		m_array = new T[m_size]();
	}

	void insert(std::vector<size_t> const &hashes, std::vector<T> &values) {
		//iterates through hashed values adding it to the filter
		for (unsigned i = 0; i < values.size(); ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
			m_array[pos] = values[i];
		}
	}

	void insert(std::vector<size_t> const &hashes, T value) {
		//add same value to the filter
		for (unsigned i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
			m_array[pos] = value;
		}
	}

	std::vector<T> query(std::vector<size_t> const &hashes) const {
		std::vector<T> values;
		for (unsigned i = 0; i < values.size(); ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
			values.push_back(m_array[pos]);
		}
		return values;
	}

	/*
	 * Returns unambiguous hit to object
	 * Returns numeric_limits<T>::max() on ambiguous collision
	 * Returns 0 on if missing element
	 */
	T at(std::vector<size_t> const &hashes) const {
		T value = 0;
		for (unsigned i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
			if(m_array[pos] == 0){
				return 0;
			}
			if (value != 0 && value != m_array[pos]) {
				value = numeric_limits<T>::max();
			} else {
				value = m_array[pos];
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
	//TODO investigate more efficient ways to do this
	T atBest(std::vector<size_t> const &hashes) const {
		google::dense_hash_map<T, unsigned> tmpHash;
		tmpHash.set_empty_key(0);
		unsigned maxCount = 0;
		T value = 0;

		for (unsigned i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
			if(m_array[pos] == 0){
				return 0;
			}
			if(tmpHash.find(m_array[pos]) != tmpHash.end()){
				++tmpHash[m_array[pos]];
			}
			if(maxCount == tmpHash[m_array[pos]]) {
				value = numeric_limits<T>::max();
			}
			else if ( maxCount < tmpHash[m_array[pos]] ){
				value = m_array[pos];
				maxCount = tmpHash[m_array[pos]];
			}
		}
		return value;
	}

	void writeHeader(ofstream &out) const {
		FileHeader header;
		strncpy(header.magic, MAGIC, 8);
		char magic[9];
		strncpy(magic, header.magic, 8);
		magic[8] = '\0';

		header.hlen = sizeof(struct FileHeader);
		header.size = m_size;
		header.dFPR = m_dFPR;
		header.nEntry = m_nEntry;
		header.tEntry = m_tEntry;

		cerr << "Writing header... magic: " << magic << " hlen: " << header.hlen
				<< " size: " << header.size << " dFPR: " << header.dFPR
				<< " nEntry: " << header.nEntry << " tEntry: " << header.tEntry
				<< endl;

		out.write(reinterpret_cast<char*>(&header), sizeof(struct FileHeader));
	}

	/*
	 * Stores the filter as a binary file to the path specified
	 * Stores uncompressed because the random data tends to
	 * compress poorly anyway
	 */
	void storeFilter(string const &filterFilePath) const {
		ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

		cerr << "Storing filter. Filter is " << m_size * sizeof(T) << "bytes."
				<< endl;

		assert(myFile);
		writeHeader(myFile);

		//write out each block
		myFile.write(reinterpret_cast<char*>(m_array), m_size * sizeof(T));

		myFile.close();
		assert(myFile);
	}

	const vector< vector <unsigned > > getSeedValues() const {
		return m_ssVal;
	}

private:

	/*
	 * Parses spaced seed string (string consisting of 1s and 0s) to vector
	 */
	inline vector< vector<unsigned> > parseSeedString(const vector<string> spacedSeeds) {
		SeedVal seeds(spacedSeeds.size());
		for(unsigned i = 0; i < spacedSeeds.size(); ++i){
			const string ss = spacedSeeds.at(i);
			for(unsigned j = 0; j < ss.size(); ++j){
				if(ss.at(i) != 0){
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
					double(numEntr) * m_ssVal.size()), double(m_ssVal.size()));
	}

	/*
	 * Calculates the optimal FPR to use based on hash functions
	 */
	double calcFPR_hashNum(int hashFunctNum) const {
		return pow(2.0, -hashFunctNum);
	}

	size_t m_size;
	T* m_array;

	double m_dFPR;
	uint64_t m_nEntry;
	uint64_t m_tEntry;

	typedef vector< vector<unsigned> > SeedVal;
	SeedVal m_ssVal;

};

#endif /* BLOOMMAP_HPP_ */
