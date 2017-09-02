/*
 * MIBloomFilter.hpp
 *
 * Agnostic of hash function used -> cannot call contains without an array of hash values
 *
 *  Created on: Jan 14, 2016
 *      Author: cjustin
 */

#ifndef MIBLOOMFILTER_HPP_
#define MIBLOOMFILTER_HPP_

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
#include <boost/shared_ptr.hpp>
#include <omp.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
static const unsigned BLOCKSIZE = 512;

/*
 * Parses spaced seed string (string consisting of 1s and 0s) to vector
 */
static vector<vector<unsigned> > parseSeedString(
		const vector<string> &spacedSeeds) {
	vector<vector<unsigned> > seeds(spacedSeeds.size(), vector<unsigned>());
	for (unsigned i = 0; i < spacedSeeds.size(); ++i) {
		const string ss = spacedSeeds.at(i);
		for (unsigned j = 0; j < ss.size(); ++j) {
			if (ss.at(j) == '0') {
				seeds[i].push_back(j);
			}
		}
	}
	return seeds;
}

///*
// * Returns an a filter size large enough to maintain an occupancy specified
// */
//inline size_t calcOptimalSize(size_t entries, unsigned hashNum,
//		double occupancy) {
//	assert(hashNum > 0);
//	return size_t(-double(entries) * double(hashNum) / log(1.0 - occupancy));
//}

/*
 * Only returns multiples of 64 for filter building purposes
 * Is an estimated size using approximations of FPR formula
 * given the number of hash functions
 * see http://en.wikipedia.org/wiki/Bloom_filter
 */
inline size_t calcOptimalSize(size_t entries, unsigned hashNum, double fpr)
{
	size_t non64ApproxVal = size_t(
			-double(entries) * double(hashNum)
					/ log(1.0 - pow(fpr, double(1 / (double(hashNum))))));

	return non64ApproxVal + (64 - non64ApproxVal % 64);
}

template<typename T>
class MIBloomFilter {
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
	 * Constructor using a prebuilt bitvector
	 */
	MIBloomFilter<T>(size_t expectedElemNum, double fpr, unsigned hashNum,
			unsigned kmerSize, sdsl::bit_vector &bv, size_t unique,
			const vector<string> seeds = vector<string>(0)) :
			m_dSize(0), m_dFPR(fpr), m_nEntry(unique), m_tEntry(
					expectedElemNum), m_hashNum(hashNum), m_kmerSize(kmerSize), m_sseeds(
					seeds) {
		cerr << "Converting bit vector to rank interleaved form" << endl;
		double start_time = omp_get_wtime();
		m_bv = sdsl::bit_vector_il<BLOCKSIZE>(bv);
		bv = sdsl::bit_vector();
		double time = omp_get_wtime() - start_time;
		cerr << "Converted bit vector to rank interleaved form " << time << "s"
				<< endl;
		if (!seeds.empty()) {
			m_ssVal = parseSeedString(m_sseeds);
			assert(m_sseeds[0].size() == kmerSize);
			for (vector<string>::const_iterator itr = m_sseeds.begin();
					itr != m_sseeds.end(); ++itr) {
//				cerr << *itr << endl;
				//check if spaced seeds are all the same length
				assert(m_kmerSize == itr->size());
			}
		}
		m_rankSupport = sdsl::rank_support_il<1>(&m_bv);
		m_dSize = getPop();
		m_data = new T[m_dSize]();
	}

	MIBloomFilter<T>(const string &filterFilePath) {
#pragma omp parallel for
		for (unsigned i = 0; i < 2; ++i) {
			if (i == 0) {
				FILE *file = fopen(filterFilePath.c_str(), "rb");
				if (file == NULL) {
#pragma omp critical(stderr)
					cerr << "file \"" << filterFilePath
							<< "\" could not be read." << endl;
					exit(1);
				}

				FileHeader header;
				if (fread(&header, sizeof(struct FileHeader), 1, file) == 1) {
#pragma omp critical(stderr)
					cerr << "Loading header..." << endl;
				} else {
#pragma omp critical(stderr)
					cerr << "Failed to Load header" << endl;
					exit(1);
				}
				char magic[9];
				strncpy(magic, header.magic, 8);
				magic[8] = '\0';
#pragma omp critical(stderr)
				cerr << "Loaded header... magic: " << magic << " hlen: "
						<< header.hlen << " size: " << header.size << " nhash: "
						<< header.nhash << " kmer: " << header.kmer << " dFPR: "
						<< header.dFPR << " nEntry: " << header.nEntry
						<< " tEntry: " << header.tEntry << endl;

				assert(strcmp(MAGIC, magic) == 0);

				m_dFPR = header.dFPR;
				m_nEntry = header.nEntry;
				m_hashNum = header.nhash;
				m_tEntry = header.tEntry;
				m_kmerSize = header.kmer;
				m_dSize = header.size;
				m_data = new T[m_dSize]();

				if (!m_sseeds.empty()) {
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
						m_sseeds.push_back(string(temp, header.kmer));
					}

					m_ssVal = parseSeedString(m_sseeds);
					assert(m_sseeds[0].size() == m_kmerSize);
					for (vector<string>::const_iterator itr = m_sseeds.begin();
							itr != m_sseeds.end(); ++itr) {
						//check if spaced seeds are all the same length
						assert(m_kmerSize == itr->size());
					}
				}

#pragma omp critical(stderr)
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
			} else {
				string bvFilename = filterFilePath + ".sdsl";
#pragma omp critical(stderr)
				cerr << "Loading sdsl interleaved bit vector from: "
						<< bvFilename << endl;
				load_from_file(m_bv, bvFilename);
				m_rankSupport = sdsl::rank_support_il<1>(&m_bv);
			}
		}

		cerr << "Bit Vector Size: " << m_bv.size() << endl;
		cerr << "Popcount: " << getPop() << endl;
	}

	/*
	 * Stores the filter as a binary file to the path specified
	 * Stores uncompressed because the random data tends to
	 * compress poorly anyway
	 */
	inline void store(string const &filterFilePath) const {

#pragma omp parallel for
		for (unsigned i = 0; i < 2; ++i) {
			if (i == 0) {
				ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

				assert(myFile);
				writeHeader(myFile);

				cerr << "Storing filter. Filter is " << m_dSize * sizeof(T)
						<< " bytes." << endl;

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
						<< " bytes" << endl;
			}
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 * ONLY REPLACE VALUE if it is larger than current value (deterministic)
	 */
	inline void insert(std::vector<size_t> const &hashes, T value) {
		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < hashes.size(); ++i) {
			size_t pos = m_rankSupport(hashes.at(i) % m_bv.size());
			setIfGreater(&m_data[pos], value);
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 * ONLY REPLACE VALUE if it is larger than current value (deterministic)
	 */
	inline void insert(const size_t *hashes, T value) {
		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t pos = m_rankSupport(hashes[i] % m_bv.size());
			setIfGreater(&m_data[pos], value);
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 * ALWAYS SETS VALUE
	 * NOT DETERMINSTIC
	 */
	inline void insert(std::vector<size_t> const &hashes, T value,
			boost::numeric::ublas::matrix<unsigned> &mat) {
		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < hashes.size(); ++i) {
			size_t pos = m_rankSupport(hashes.at(i) % m_bv.size());
			setVal(&m_data[pos], value, mat);
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 * ALWAYS SETS VALUE
	 * NOT DETERMINSTIC
	 */
	inline void insert(const size_t *hashes, T value,
			boost::numeric::ublas::matrix<unsigned> &mat) {
		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t pos = m_rankSupport(hashes[i] % m_bv.size());
			setVal(&m_data[pos], value, mat);
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 * Replaces values according to collisionID hashtable (for thread safety)
	 */
	inline void insert(std::vector<size_t> const &hashes, T value,
			const vector<boost::shared_ptr<google::dense_hash_map<T, T> > > &colliIDs) {
		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < hashes.size(); ++i) {
			size_t pos = m_rankSupport(hashes.at(i) % m_bv.size());
			setVal(&m_data[pos], value, colliIDs);
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 * Replaces values according to collisionID hashtable (for thread safety)
	 */
	inline void insert(const size_t *hashes, T value,
			const vector<boost::shared_ptr<google::dense_hash_map<T, T> > > &colliIDs) {
		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t pos = m_rankSupport(hashes[i] % m_bv.size());
			setVal(&m_data[pos], value, colliIDs);
		}
	}

	/*
	 * Mutates value vector to contain values in bloom map
	 * Assumes vector is the size of hashes
	 */
	inline void query(std::vector<size_t> const &hashes,
			vector<T> &values) const {
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
	 * Mutates value vector to contain values in bloom map
	 * Assumes vector is the size of hashes
	 */
	inline void query(const size_t *hashes, vector<T> &values) const {
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_bv.size();
			if (m_bv[pos] == 0) {
				values[i] = 0;
				continue;
			}
			values[i] = m_data[m_rankSupport(pos)];
		}
	}

	//TODO optimize
	inline vector<T> at(const size_t *hashes,
			const vector<boost::shared_ptr<vector<T> > > &colliIDs,
			unsigned maxMiss, unsigned &misses) {
		vector<T> results;
		results.reserve(m_hashNum);
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_bv.size();
			if (m_bv[pos] == 0) {
				++misses;
				if (misses > maxMiss) {
					return results;
				}
			}
		}

		google::dense_hash_map<T, unsigned> tmpHash;
		tmpHash.set_empty_key(0);
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_bv.size();
			if (m_bv[pos] != 0) {
				size_t rankPos = m_rankSupport(pos);
				T currID = m_data[rankPos];
				if (currID
						!= std::numeric_limits<T>::max() && colliIDs[currID] != NULL) {
					if (tmpHash.find(currID) != tmpHash.end()) {
						++tmpHash[currID];
					} else {
						tmpHash[currID] = 1;
					}
				}
			}
		}

		google::dense_hash_map<T, unsigned> tmpHash2;
		tmpHash2.set_empty_key(0);
		unsigned maxCount = 0;
		for (typename google::dense_hash_map<T, unsigned>::iterator itr =
				tmpHash.begin(); itr != tmpHash.end(); ++itr) {
			for (unsigned i = 0; i < colliIDs[itr->first]->size(); ++i) {
				T id = colliIDs[itr->first]->at(i);
				if (tmpHash2.find(id) != tmpHash2.end()) {
					tmpHash2[id] += itr->second;
					if (maxCount < tmpHash2[id]) {
						maxCount = tmpHash2[id];
					}
				} else {
					tmpHash2[id] = itr->second;
					if (maxCount < tmpHash2[id]) {
						maxCount = tmpHash2[id];
					}
				}
			}
		}
		for (typename google::dense_hash_map<T, unsigned>::iterator itr =
				tmpHash2.begin(); itr != tmpHash2.end(); ++itr) {
			if (maxCount == itr->second) {
				results.push_back(itr->first);
			}
		}
		return results;
	}

	/*
	 * Returns unambiguous hit to object
	 * Returns best hit on ambiguous collisions
	 * Returns numeric_limits<T>::max() on completely ambiguous collision
	 * Returns 0 on if missing element
	 * Uses colliIDs to resolve ambiguities
	 */
	//TODO optimize
	inline vector<T> at(std::vector<size_t> const &hashes,
			const vector<boost::shared_ptr<vector<T> > > &colliIDs,
			unsigned maxMiss, unsigned &misses) const {
		vector<T> results;
		results.reserve(hashes.size());
		for (unsigned i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_bv.size();
			if (m_bv[pos] == 0) {
				++misses;
				if (misses > maxMiss) {
					return results;
				}
			}
		}

		google::dense_hash_map<T, unsigned> tmpHash;
		tmpHash.set_empty_key(0);
		for (unsigned i = 0; i < hashes.size(); ++i) {
			size_t pos = hashes.at(i) % m_bv.size();
			if (m_bv[pos] != 0) {
				size_t rankPos = m_rankSupport(pos);
				T currID = m_data[rankPos];
				if (currID
						!= std::numeric_limits<T>::max() && colliIDs[currID] != NULL) {
					if (tmpHash.find(currID) != tmpHash.end()) {
						++tmpHash[currID];
					} else {
						tmpHash[currID] = 1;
					}
				}
			}
		}

		google::dense_hash_map<T, unsigned> tmpHash2;
		tmpHash2.set_empty_key(0);
		unsigned maxCount = 0;
		for (typename google::dense_hash_map<T, unsigned>::iterator itr =
				tmpHash.begin(); itr != tmpHash.end(); ++itr) {
			for (unsigned i = 0; i < colliIDs[itr->first]->size(); ++i) {
				T id = colliIDs[itr->first]->at(i);
				if (tmpHash2.find(id) != tmpHash2.end()) {
					tmpHash2[id] += itr->second;
					if (maxCount < tmpHash2[id]) {
						maxCount = tmpHash2[id];
					}
				} else {
					tmpHash2[id] = itr->second;
					if (maxCount < tmpHash2[id]) {
						maxCount = tmpHash2[id];
					}
				}
			}
		}
		for (typename google::dense_hash_map<T, unsigned>::iterator itr =
				tmpHash2.begin(); itr != tmpHash2.end(); ++itr) {
			if (maxCount == itr->second) {
				results.push_back(itr->first);
			}
		}
		return results;
	}

	inline const vector<vector<unsigned> > &getSeedValues() const {
		return m_ssVal;
	}

	inline unsigned getKmerSize() {
		return m_kmerSize;
	}

	inline unsigned getHashNum() {
		return m_hashNum;
	}

	/*
	 * Return FPR based on popcount
	 */
	inline double getFPR() const {
		return pow(double(getPop()) / double(m_bv.size()),
				double(m_hashNum));
	}

	/*
	 * Return FPR based on popcount and minimum number of matches for a hit
	 */
	inline double getFPR(unsigned allowedMiss) const {
		assert(allowedMiss <= m_hashNum);
		double cumulativeProb = 0;
		double popCount = getPop();
		double p = popCount / double(m_bv.size());
		cerr << p << endl;
		for (unsigned i = m_hashNum - allowedMiss; i <= m_hashNum;
				++i) {
			cumulativeProb += double(nChoosek(m_hashNum, i)) * pow(p, i)
					* pow(1.0 - p, (m_hashNum - i));
		}
		return (cumulativeProb);
	}

	/*
	 * Return FPR based on number of inserted elements
	 */
	inline double getFPR_numEle() const {
		assert(m_nEntry > 0);
		return calcFPR_numInserted(m_nEntry);
	}

	inline size_t getPop() const {
		size_t index = m_bv.size() - 1;
		while (m_bv[index] == 0) {
			--index;
		}
		return m_rankSupport(index - 1) + 1;
	}

	inline size_t getUniqueEntries() const {
		return m_nEntry;
	}

	~MIBloomFilter() {
		delete[] m_data;
	}

private:

	/*
	 * Helper function for header storage
	 */
	inline void writeHeader(ofstream &out) const {
		FileHeader header;
		strncpy(header.magic, MAGIC, 8);
		char magic[9];
		strncpy(magic, header.magic, 8);
		magic[8] = '\0';

		header.hlen = sizeof(struct FileHeader) + m_kmerSize * m_sseeds.size();
		header.kmer = m_kmerSize;
		header.size = m_dSize;
		header.nhash = m_hashNum;
		header.dFPR = m_dFPR;
		header.nEntry = m_nEntry;
		header.tEntry = m_tEntry;

		cerr << "Writing header... magic: " << magic << " hlen: " << header.hlen
				<< " nhash: " << header.nhash << " size: " << header.size
				<< " dFPR: " << header.dFPR << " nEntry: " << header.nEntry
				<< " tEntry: " << header.tEntry << endl;

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
	inline static unsigned calcOptiHashNum(double fpr) {
		return unsigned(-log(fpr) / log(2));
	}

	/*
	 * Calculate FPR based on hash functions, size and number of entries
	 * see http://en.wikipedia.org/wiki/Bloom_filter
	 */
	inline double calcFPR_numInserted(size_t numEntr) const {
		return pow(
				1.0
						- pow(1.0 - 1.0 / double(m_bv.size()),
								double(numEntr) * double(m_hashNum)),
				double(m_hashNum));
	}

	/*
	 * Calculates the optimal FPR to use based on hash functions
	 */
	inline double calcFPR_hashNum(int hashFunctNum) const {
		return pow(2.0, -hashFunctNum);
	}

	/*
	 * Lock free cas value setting for larger element
	 */
	inline void setIfGreater(T *val, T newVal) {
		T oldValue;
		do {
			oldValue = *val;
			if (oldValue >= newVal) {
				break;
			}
		} while (!__sync_bool_compare_and_swap(val, oldValue, newVal));
	}

	inline void setVal(T *val, T newVal,
			boost::numeric::ublas::matrix<unsigned> &mat) {
		T oldValue;
		do {
			oldValue = *val;
			if (oldValue == newVal)
				break;
		} while (!__sync_bool_compare_and_swap(val, oldValue, newVal));
		if (oldValue != 0 && oldValue != newVal) {
#pragma omp atomic
			++mat(oldValue - 1, newVal - 1);
		}
	}

	/*
	 * Lock free cas value setting for larger element
	 * check if it is a collision, setting it accordingly
	 */
	inline void setVal(T *val, T newVal,
			const vector<boost::shared_ptr<google::dense_hash_map<T, T> > > &colliIDs) {
		T oldValue;
		T insertValue;
		do {
			oldValue = *val;
			insertValue = newVal;
			if (oldValue != 0) {
				//NO NET CHANGE
				if (oldValue == insertValue
						|| oldValue == numeric_limits<T>::max()) {
					break;
				}
				//check if oldValue and new value have a collision ID together
				typename google::dense_hash_map<T, T>::iterator colliPtr =
						colliIDs[oldValue]->find(insertValue);
				if (colliPtr != colliIDs[oldValue]->end()) {
					insertValue = colliPtr->second;
					//NO NET CHANGE
					if (oldValue == insertValue) {
						break;
					}
				}
				//IF UNRESOLVED COLLISION
				else {
					insertValue = numeric_limits<T>::max();
				}
			}
		} while (!__sync_bool_compare_and_swap(val, oldValue, insertValue));
	}

	inline unsigned nChoosek(unsigned n, unsigned k) const {
		if (k > n)
			return 0;
		if (k * 2 > n)
			k = n - k;
		if (k == 0)
			return 1;

		int result = n;
		for (unsigned i = 2; i <= k; ++i) {
			result *= (n - i + 1);
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
	unsigned m_hashNum;
	unsigned m_kmerSize;

	typedef vector<vector<unsigned> > SeedVal;
	vector<string> m_sseeds;
	SeedVal m_ssVal;
	const char* MAGIC = "MIBLOOMF";
};

#endif /* MIBLOOMFILTER_HPP_ */
