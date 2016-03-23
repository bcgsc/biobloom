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
#include "rolling.h"

using namespace std;

template<typename T>
class BloomMap {

public:

	struct FileHeader {
		char magic[8];
		uint32_t hlen;
		uint64_t size;
		uint32_t nhash;
		uint32_t kmer;
		double dFPR;
		uint64_t nEntry;
		uint64_t tEntry;
	};

	BloomMap<T>(size_t filterSize, unsigned hashNum, unsigned kmerSize) :
			m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize), m_dFPR(
					0), m_nEntry(0), m_tEntry(0) {
		m_array = new T[m_size]();
	}

	/* De novo filter constructor.
	 * Allocates a filter size based on the number of expected elements and FPR
	 *
	 * If hashNum is set to 0, an optimal value is computed based on the FPR
	 */
	BloomMap<T>(size_t expectedElemNum, double fpr, unsigned hashNum,
			unsigned kmerSize) :
			m_size(0), m_hashNum(hashNum), m_kmerSize(kmerSize), m_dFPR(fpr), m_nEntry(
					0), m_tEntry(0) {
		if (m_hashNum == 0) {
			m_hashNum = calcOptiHashNum(m_dFPR);
		}
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
					<< " does not match size given by its information file. Size: "
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

        cerr << "Loading header... magic: " <<
            magic << " hlen: " <<
            header.hlen << " size: " <<
            header.size << " nhash: " <<
            header.nhash << " kmer: " <<
           header.kmer << " dFPR: " <<
            header.dFPR << " aFPR: " <<
            header.nEntry << " tEntry: " <<
            header.tEntry << endl;

		m_size = header.size;
		m_array = new T[m_size]();
		m_hashNum = header.nhash;
		m_kmerSize = header.kmer;
	}

	T& operator[](size_t i) {
		assert(i < m_size);
		return m_array[i];
	}

	const T& operator[](size_t i) const {
		assert(i < m_size);
		return m_array[i];
	}

	void insert(std::vector<size_t> const &hashes, std::vector<T> &values) {
		assert(hashes.size() == m_hashNum);
		//iterates through hashed values adding it to the filter
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
			m_array[pos] = values[i];
		}
	}

	void insert(std::vector<size_t> const &hashes, T value) {
		assert(hashes.size() == m_hashNum);
		//add same value to the filter
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
			m_array[pos] = value;
		}
	}

	std::vector<T> query(std::vector<size_t> const &hashes) {
		assert(hashes.size() == m_hashNum);
		std::vector<T> values;

		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
			values.push_back(m_array[pos]);

		}
		return values;
	}

	void writeHeader(ofstream &out) const {
		FileHeader header;
		strncpy(header.magic, "BlOOMFXX", 8);
		char magic[9];
		strncpy(magic, header.magic, 8);
		magic[8] = '\0';

		header.hlen = sizeof(struct FileHeader);
		header.size = m_size;
		header.nhash = m_hashNum;
		header.kmer = m_kmerSize;
		header.dFPR = m_dFPR;
		header.nEntry = m_nEntry;
		header.tEntry = m_tEntry;

        cerr << "Writing header... magic: "
            << magic << " hlen: "
            << header.hlen << " size: "
            << header.size << " nhash: "
            << header.nhash << " kmer: "
           << header.kmer << " dFPR: "
            << header.dFPR << " aFPR: "
           << header.nEntry << " tEntry: "
            << header.tEntry << endl;

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

private:

	/*
	 * Only returns multiples of 64 for filter building purposes
	 * Is an estimated size using approximations of FPR formula
	 * given the number of hash functions
	 */
	size_t calcOptimalSize(size_t entries, double fpr) const {
		size_t non64ApproxVal = size_t(
				-double(entries) * double(m_hashNum)
						/ log(1.0 - pow(fpr, double(1 / double(m_hashNum)))));

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
		return pow(
				1.0
						- pow(1.0 - 1.0 / double(m_size),
								double(numEntr) * m_hashNum), double(m_hashNum));
	}

	/*
	 * Calculates the optimal FPR to use based on hash functions
	 */
	double calcFPR_hashNum(unsigned hashFunctNum) const {
		return pow(2, -hashFunctNum);
	}

	size_t m_size;
	unsigned m_hashNum;
	T* m_array;

	unsigned m_kmerSize;
	double m_dFPR;
	uint64_t m_nEntry;
	uint64_t m_tEntry;

};

#endif /* BLOOMMAP_HPP_ */
