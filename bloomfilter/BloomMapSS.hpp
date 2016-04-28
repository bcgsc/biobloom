/*
 * BloomMap.hpp
 *
 *  Created on: Dec 17, 2015
 *      Author: gjahesh
 */

#ifndef BLOOMMAPSS_HPP_
#define BLOOMMAPSS_HPP_

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
#include "BloomMapSSBitVec.hpp"

using namespace std;
template<typename T>
class BloomMapSSBitVec;

/*
 * Returns an a filter size large enough to maintain an occupancy specified
 */
inline size_t calcOptimalSize(size_t entries, unsigned hashNum, double occupancy) {
	assert(hashNum > 0);
	return size_t(-double(entries) * double(hashNum) / log(1.0 - occupancy));
}

/*
 * Parses spaced seed string (string consisting of 1s and 0s) to vector
 */
inline vector<vector<unsigned> > parseSeedString(const vector<string> &spacedSeeds) {
	vector<vector<unsigned> > seeds(spacedSeeds.size());
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

template<typename T>
class BloomMapSS {
	const char* MAGIC = "BLOOMMSS";
	friend class BloomMapSSBitVec<T>;

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

	/* De novo filter constructor.
	 * Allocates a filter size based on the number of expected elements and FPR
	 *
	 * If hashNum is set to 0, an optimal value is computed based on the FPR
	 */
	BloomMapSS<T>(size_t expectedElemNum, double fpr, vector<string> seeds) :
			m_size(0), m_dFPR(fpr), m_nEntry(0), m_tEntry(expectedElemNum), m_sseeds(seeds), m_kmerSize(
					m_sseeds[0].size()), m_ssVal(parseSeedString(m_sseeds)) {
		for (vector<string>::const_iterator itr = m_sseeds.begin();
				itr != m_sseeds.end(); ++itr) {
			//check if spaced seeds are all the same length
			assert(m_kmerSize == itr->size());
		}
		if (m_size == 0) {
			m_size = calcOptimalSize(expectedElemNum, m_dFPR, m_sseeds.size());
		}
		m_array = new T[m_size]();
	}

	BloomMapSS<T>(const string &filterFilePath) {
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
		m_array = new T[m_size]();

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

		size_t countRead = fread(m_array, fileSize, 1, file);
		if (countRead != 1 && fclose(file) != 0) {
			cerr << "file \"" << filterFilePath << "\" could not be read."
					<< endl;
			exit(1);
		}
		cerr << "FPR based on PopCount: " << getFPR() << endl;
		cerr << "FPR based on Unique Elements: " << getFPR_numEle() << endl;
	}

	/*
	 * Stores the filter as a binary file to the path specified
	 * Stores uncompressed because the random data tends to
	 * compress poorly anyway
	 */
	void storeFilter(string const &filterFilePath) const {
		ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

		assert(myFile);
		writeHeader(myFile);

		cerr << "Storing filter. Filter is " << m_size * sizeof(T) << "bytes."
				<< endl;

		//write out each block
		myFile.write(reinterpret_cast<char*>(m_array), m_size * sizeof(T));

		myFile.close();
		assert(myFile);
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

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 * Returns if already inserted
	 * Inserts ambiguous ID when collision occurs
	 */
	bool insertAndCheck(std::vector<size_t> const &hashes, T value, size_t &collisionCount) {
		//iterates through hashed values adding it to the filter
		bool found = true;
		for (size_t i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_size;
			if(m_array[pos] == 0){
				found = false;
			}
			else if( m_array[pos] != value){
//				m_array[pos] = numeric_limits<T>::max();
				++collisionCount;
			}
			m_array[pos] = value;
		}
		return found;
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
	vector<T> at(std::vector<size_t> const &hashes) const {
		vector<T> values(hashes.size());
		for (unsigned i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_size;
			values[i] = m_array[pos];
		}
		return values;
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
			if(m_array[pos] == 0){
				return 0;
			}
			if (tmpHash.find(m_array[pos]) != tmpHash.end()) {
				++tmpHash[m_array[pos]];
			} else {
				tmpHash[m_array[pos]] = 1;
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
			if(m_array[pos] == 0){
				++miss;
				continue;
			}
			if (tmpHash.find(m_array[pos]) != tmpHash.end()) {
				++tmpHash[m_array[pos]];
			} else {
				tmpHash[m_array[pos]] = 1;
			}
			if(maxCount == tmpHash[m_array[pos]]) {
				value = numeric_limits<T>::max();
			}
			else if ( maxCount < tmpHash[m_array[pos]] ){
				value = m_array[pos];
				maxCount = tmpHash[m_array[pos]];
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

	const vector< string > &getSeedStrings() const {
		return m_sseeds;
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
	 * Return FPR based on popcount and minimum number of matches for a hit
	 */
	double getFPR(unsigned minNum) const {
		return pow(double(getPop())/double(m_size), double(minNum));
	}

	/*
	 * Return FPR based on number of inserted elements
	 */
	double getFPR_numEle() const {
		assert(m_nEntry > 0);
		return calcFPR_numInserted(m_nEntry);
	}

	size_t getPop() const {
		size_t i, popBF = 0;
//#pragma omp parallel for reduction(+:popBF)
		for (i = 0; i < m_size; i++)
			popBF = popBF + (m_array[i] != 0);
		return popBF;
	}

	size_t size() const {
		return m_size;
	}

	double getDesiredFPR() const {
		return m_dFPR;
	}

	size_t getUniqueEntries() const {
		return m_nEntry;
	}

	size_t getTotalEntries() const {
		return m_tEntry;
	}

	void setUnique(size_t count){
		m_nEntry = count;
	}

	~BloomMapSS() {
		delete[] m_array;
	}

private:

	size_t m_size;
	T* m_array;

	double m_dFPR;
	uint64_t m_nEntry;
	uint64_t m_tEntry;
	vector<string> m_sseeds;
	unsigned m_kmerSize;

	typedef vector<vector<unsigned> > SeedVal;
	SeedVal m_ssVal;

//	void loadHeader(FILE *file) {
//
//		FileHeader header;
//		if (fread(&header, sizeof(struct FileHeader), 1, file) == 1) {
//			cerr << "Loading header..." << endl;
//		} else {
//			cerr << "Failed to header" << endl;
//			exit(1);
//		}
//		char magic[9];
//		strncpy(magic, header.magic, 8);
//		assert(string(MAGIC) == string(header.magic));
//		magic[8] = '\0';
//
//		cerr << "Loaded header... magic: " << magic << " hlen: " << header.hlen
//				<< " size: " << header.size << " nhash: " << header.nhash
//				<< " kmer: " << header.kmer << " dFPR: " << header.dFPR
//				<< " aFPR: " << header.nEntry << " tEntry: " << header.tEntry
//				<< endl;
//
//		m_sseeds = vector<string>(header.size);
//
//        //load seeds
//		for (unsigned i = 0; i < header.nhash; ++i) {
//			char temp[header.kmer];
//
//			if (fread(temp, header.kmer, 1, file) != 1) {
//				cerr << "Failed to load spaced seed string" << endl;
//				exit(1);
//			}
//			else{
//				cerr << "Spaced Seed " << i <<": " << string(temp, header.kmer) << endl;
//			}
//			m_sseeds[i] = string(temp, header.kmer);
//		}
//		m_dFPR = header.dFPR;
//		m_nEntry = header.nEntry;
//		m_tEntry = header.tEntry;
//		m_kmerSize = header.kmer;
//		m_size = header.size;
//		m_array = new T[m_size]();
//	}

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
};

#endif /* BLOOMMAPSS_HPP_ */
