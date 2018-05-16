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
#include <sdsl/bit_vector_il.hpp>
#include <sdsl/rank_support.hpp>
#include "Options.h"
#include "Common/Options.h"
#include <omp.h>
#include <algorithm>    // std::random_shuffle
#include <boost/math/distributions/binomial.hpp>
#include <google/dense_hash_set>
#include <google/dense_hash_map>

using namespace std;

template<typename T>
class MIBloomFilter {
public:
	static const T s_mask = 1 << (sizeof(T) * 8 - 1);
	static const T s_antiMask = (T)~s_mask;

	static const unsigned BLOCKSIZE = 512;
	//static methods
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

	//helper methods
	//calculates the per frame probability of a random match for single value
	static double calcProbSingleFrame(double occupancy, unsigned hashNum,
			double freq, unsigned allowedMisses = 0) {
		double probTotal = 0.0;
		for (unsigned i = hashNum - allowedMisses; i <= hashNum; i++) {
			double prob = nChoosek(hashNum, i);
			prob *= pow(occupancy, i);
			prob *= pow(1.0 - occupancy, hashNum - i);
			prob *= (1.0 - pow(1.0 - freq, i));
			probTotal += prob;
		}
		return probTotal;
	}

	//calculates the per frame probability of a random multi match given a significant result
	static double calcProbMultiMatchSingleFrame(double occupancy, unsigned hashNum,
			double freq) {
		double prob = 1.0
				- pow(1.0 - freq, hashNum * (1 + occupancy / log(1 - occupancy)));
		return prob;
	}

	/*
	 * Preconditions:
	 * 	frameProbs but be equal in size to multiMatchProbs
	 * 	frameProbs must be preallocated to correct size (number of ids + 1)
	 * Max value is the largest value seen in your set of possible values
	 */
	static void calcFrameProbs(MIBloomFilter<T> &miBF,
			vector<double> &frameProbs, vector<double> &multiMatchProbs) {
		double occupancy = double(miBF.getPop()) / double(miBF.size());
		unsigned hashNum = miBF.getHashNum();
		vector<size_t> countTable = vector<size_t>(frameProbs.size(), 0);
		miBF.getIDCounts(countTable);
		size_t sum = 0;
		for (vector<size_t>::const_iterator itr = countTable.begin();
				itr != countTable.end(); ++itr) {
			sum += *itr;
		}
		for (size_t i = 0; i < countTable.size(); ++i) {
			frameProbs[i] = calcProbSingleFrame(occupancy, hashNum,
					double(countTable[i]) / double(sum));
			multiMatchProbs[i] = calcProbMultiMatchSingleFrame(occupancy,
					hashNum, double(countTable[i]) / double(sum));
		}
	}

	/*
	 * Utility function for querying MiBF for a sequence in a hash itr object
	 * Computes probability of a match for each possible value queried.
	 * In this scheme saturated regions are ignored.
	 *
	 * Parameters:
	 * perFrameProb - per frame probability of each possible value
	 * alpha - significance threshold
	 * itr - hash iterator of type H (ntHash)
	 * maxPos - max number of positions to move on the sequence
	 *
	 * Returns significant values (smaller than alpha threshold)
	 */
	//TODO return saturated frame counts?
	//TODO return pVals?
	template<typename H>
	vector<pair<ID, double>> query(H &itr, const vector<double> &perFrameProb,
			const vector<double> &perMultiMatchFrameProb, double alpha = 0.0001,
			double multimapAlpha = 0.001, unsigned maxMiss = 0, size_t maxPos =
					numeric_limits<size_t>::max()) {
		unsigned evaluatedSeeds = 0;
		unsigned totalCount = 0;
		unsigned saturatedCount = 0;

		google::dense_hash_map<T, unsigned> counts;
		counts.set_empty_key(0);
		while (itr != itr.end() && itr.pos() < maxPos) {
			bool saturated = true;
			vector<T> results = at(*itr, saturated, maxMiss);
			//to determine if already added for this frame
			google::dense_hash_set<T> tempIDs;
			tempIDs.set_empty_key(0);

			if (!saturated) {
				for (typename vector<T>::const_iterator j = results.begin();
						j != results.end(); j++) {
					if (*j != 0) {
						if (tempIDs.find(*j) == tempIDs.end()) {
							typename google::dense_hash_map<T, unsigned>::iterator tempItr =
									counts.find(*j);
							assert(*j > 0);
							if (tempItr == counts.end()) {
								counts[*j] = 1;
							} else {
								++(tempItr->second);
							}
							tempIDs.insert(*j);
						}
					}
				}
				++evaluatedSeeds;
			}
			else{
				++saturatedCount;
			}
			++totalCount;
			++itr;
		}

		//potential signifResults
		vector<pair<T, double>> potSignifResults;
		vector<pair<T, double>> signifResults;

		double adjustedPValThreshold = 1.0
				- pow(1.0 - alpha, 1.0 / double(perFrameProb.size() - 1));
		T bestSignifVal = counts.begin()->first;
		for (typename google::dense_hash_map<T, unsigned>::const_iterator itr =
				counts.begin(); itr != counts.end(); itr++) {
			//TODO use complement cdf? so I don't have to subtract?
			boost::math::binomial bin(evaluatedSeeds, 1.0 - perFrameProb.at(itr->first));
			double cumProb = cdf(bin, evaluatedSeeds - itr->second);
			if (adjustedPValThreshold > cumProb) {
				if (counts[bestSignifVal] < counts[itr->first]) {
					bestSignifVal = itr->first;
				}
				potSignifResults.push_back(pair<T, double>(itr->first, cumProb));
			}
		}

//		adjustedPValThreshold = 1.0
//				- pow(1.0 - multimapAlpha, 1.0 / double(perFrameProb.size() - 1));
		assert(perMultiMatchFrameProb.size());
		assert(multimapAlpha);
		//TODO: generalized because this assumes a = 0, fix me?
//		for (typename vector<pair<T, double>>::const_iterator itr = potSignifResults.begin();
//				itr != potSignifResults.end(); ++itr) {
//			//compute single frame prob
//			boost::math::binomial bin(counts[bestSignifVal], 1.0 - perMultiMatchFrameProb.at(itr->first));
//			double cumProb = cdf(bin, counts[bestSignifVal] - counts[itr->first]);
//			if (adjustedPValThreshold > cumProb) {
//				signifResults.push_back(*itr);
//			}
//		}

		for (typename vector<pair<T, double>>::const_iterator itr = potSignifResults.begin();
				itr != potSignifResults.end(); ++itr) {
			if (counts[bestSignifVal] <= counts[itr->first]) {
				signifResults.push_back(*itr);
			}
		}

		//TODO detect repeat sequence
//
//		if (signifResults.size() == 0 && saturatedCount) {
//			boost::math::binomial bin(totalCount, 1.0 - m_probSaturated);
//			double cumProb = cdf(bin, totalCount - saturatedCount);
//			if (adjustedPValThreshold > cumProb) {
//				//0 is empty
//				signifResults.push_back(pair<T, double>(0, saturatedCount));
//			}
//		}
		sort(signifResults.begin(), signifResults.end(), sortbysec);
		//Best hit considered the class with the most hits and lowest pValue
		return signifResults;
	}

	//TODO refine extraFrameLimit
	template<typename H>
	vector<pair<T, double>> query(H &itr, const vector<unsigned> &minCount,
			const vector<double> &perFrameProb, unsigned extraFrameLimit = 50,
			unsigned extraCount = 5) {
		vector<pair<T, double>> signifResults;
		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned frameCount = 0;
		unsigned secondBestCount = 0;

		google::dense_hash_map<T, unsigned> counts;
		counts.set_empty_key(0);
		google::dense_hash_set<T> candidateMatch;
		candidateMatch.set_empty_key(0);
		bool candidateFound = false;

		while (itr != itr.end() && !candidateFound) {
			unsigned count = 0;
			vector<size_t> rankPos = atPos(*itr, count);
			if (count == m_hashNum) {
				vector<ID> results(m_hashNum);
				for (unsigned i = 0; i < m_hashNum; ++i) {
					results[i] = m_data[rankPos[i]];
				}
				google::dense_hash_set<T> tempIDs;
				tempIDs.set_empty_key(0);
				for (typename vector<T>::const_iterator j = results.begin();
						j != results.end(); j++) {
					if (*j != 0) {
//						bool saturated = true;
						T result = *j;
						//check for saturation
						if (result > s_mask) {
							result = *j & s_antiMask;
						}
						if (tempIDs.find(result) == tempIDs.end()) {
							typename google::dense_hash_map<T, unsigned>::iterator tempItr =
									counts.find(result);
							assert(result > 0);
							if (tempItr == counts.end()) {
								counts[result] = 1;
							} else {
								//check is count is exceeded
								if (minCount[result] <= ++tempItr->second) {
									if (tempItr->second > bestCount) {
										bestCount = counts[result];
									} else if (counts[result] > secondBestCount) {
										secondBestCount = counts[result];
									}
									candidateMatch.insert(result);
								}
							}
							tempIDs.insert(result);
						}
					}
				}
				if (bestCount <= secondBestCount + extraCount) {
					extraFrame = 0;
				}
				if (bestCount && bestCount > secondBestCount) {
					if (extraFrameLimit < extraFrame++) {
						candidateFound = true;
					}
				}
			}
			++frameCount;
			++itr;
		}
		for (typename google::dense_hash_set<T>::const_iterator candidates =
				candidateMatch.begin(); candidates != candidateMatch.end();
				candidates++) {
			if (bestCount <= counts[*candidates] + extraCount) {
				//todo record not counts but pVals?
				signifResults.push_back(
						pair<T, double>(*candidates,
								double(counts[*candidates])
										+ perFrameProb[*candidates]));
			}
		}
		sort(signifResults.begin(), signifResults.end(), sortbysec);
		return signifResults;
	}

	template<typename H>
	vector<pair<ID, double>> query(H &itr1, H &itr2,
			const vector<unsigned> &minCount,
			const vector<double> &perFrameProb, unsigned extraFrameLimit = 10, unsigned extraCount = 1) {

		vector<pair<T, double>> signifResults;
		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned frameCount = 0;
		unsigned secondBestCount = 0;

		google::dense_hash_map<T, unsigned> counts;
		counts.set_empty_key(0);
		google::dense_hash_set<T> candidateMatch;
		candidateMatch.set_empty_key(0);
		bool candidateFound = false;

		while ((itr1 != itr1.end() && itr2 != itr2.end()) && !candidateFound) {
			H &itr = frameCount % 2 == 0 && itr1 != itr1.end() ? itr1 :
						frameCount % 2 == 1 && itr2 != itr2.end() ? itr2 : itr1;
			unsigned count = 0;
			vector<size_t> rankPos = atPos(*itr, count);
			if (count == m_hashNum) {
				vector<ID> results(m_hashNum);
				for (unsigned i = 0; i < m_hashNum; ++i) {
					results[i] = m_data[rankPos[i]];
				}
				google::dense_hash_set<T> tempIDs;
				tempIDs.set_empty_key(0);
				for (typename vector<T>::const_iterator j = results.begin();
						j != results.end(); j++) {
					if (*j != 0) {
//						bool saturated = true;
						T result = *j;
						//check for saturation
						if (result > s_mask) {
							result = *j & s_antiMask;
						}
						if (tempIDs.find(result) == tempIDs.end()) {
							typename google::dense_hash_map<T, unsigned>::iterator tempItr =
									counts.find(result);
							assert(result > 0);
							if (tempItr == counts.end()) {
								counts[result] = 1;
							} else {
								//check is count is exceeded
								if (minCount[result] <= ++tempItr->second) {
									if (tempItr->second > bestCount) {
										bestCount = counts[result];
									} else if (counts[result] > secondBestCount) {
										secondBestCount = counts[result];
									}
									candidateMatch.insert(result);
								}
							}
							tempIDs.insert(result);
						}
					}
				}
				if (bestCount <= secondBestCount + extraCount) {
					extraFrame = 0;
				}
				if (bestCount && bestCount > secondBestCount) {
					if (extraFrameLimit < extraFrame++) {
						candidateFound = true;
					}
				}
			}
			++frameCount;
			++itr;
		}
		for (typename google::dense_hash_set<T>::const_iterator candidates =
				candidateMatch.begin(); candidates != candidateMatch.end(); candidates++) {
			if (bestCount <= counts[*candidates] + extraCount) {
				//todo record not counts but pVals?
				signifResults.push_back(pair<T, double>(*candidates, double(counts[*candidates]) + perFrameProb[*candidates]));
			}
		}
		sort(signifResults.begin(), signifResults.end(), sortbysec);
		return signifResults;
	}

	/*
	 * Returns an a filter size large enough to maintain an occupancy specified
	 */
	static size_t calcOptimalSize(size_t entries, unsigned hashNum,
			double occupancy) {
		size_t non64ApproxVal = size_t(
				-double(entries) * double(hashNum) / log(1.0 - occupancy));
		return non64ApproxVal + (64 - non64ApproxVal % 64);
	}

	/*
	 * Inserts a set of hash values into an sdsl bitvector and returns the number of collisions
	 * Thread safe on the bv, though return values will not be the same run to run
	 */
	static unsigned insert(sdsl::bit_vector &bv, uint64_t * hashValues,
			unsigned hashNum) {
		unsigned colliCount = 0;
		for (size_t i = 0; i < hashNum; ++i) {
			size_t pos = hashValues[i] % bv.size();
			uint64_t *dataIndex = bv.data() + (pos >> 6);
			uint64_t bitMaskValue = (uint64_t) 1 << (pos & 0x3F);
			colliCount += __sync_fetch_and_or(dataIndex, bitMaskValue)
					>> (pos & 0x3F) & 1;
		}
		return colliCount;
	}

	//TODO: include allowed miss in header
#pragma pack(1) //to maintain consistent values across platforms
	struct FileHeader {
		char magic[8];
		uint32_t hlen;	//header length (including spaced seeds)
		uint64_t size;
		uint32_t nhash;
		uint32_t kmer;
//		uint8_t allowedMiss;
	};

	/*
	 * Constructor using a prebuilt bitvector
	 */
	MIBloomFilter<T>(unsigned hashNum, unsigned kmerSize, sdsl::bit_vector &bv,
			const vector<string> seeds = vector<string>(0)) :
			m_dSize(0), m_hashNum(hashNum), m_kmerSize(kmerSize), m_sseeds(
					seeds), m_probSaturated(0) {
		m_bv = sdsl::bit_vector_il<BLOCKSIZE>(bv);
		bv = sdsl::bit_vector();
		if (!seeds.empty()) {
			m_ssVal = parseSeedString(m_sseeds);
			assert(m_sseeds[0].size() == kmerSize);
			for (vector<string>::const_iterator itr = m_sseeds.begin();
					itr != m_sseeds.end(); ++itr) {
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
						<< header.nhash << " kmer: " << header.kmer << endl;

				m_hashNum = header.nhash;
				m_kmerSize = header.kmer;
				m_dSize = header.size;
				m_data = new T[m_dSize]();

				if (header.hlen > sizeof(struct FileHeader)) {
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
		//TODO make more streamlined
		m_probSaturated = pow(double(getPopSaturated())/double(getPop()),m_hashNum);
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
	 * Insert using an existing iterator (returning array of size_t)
	 */
	template<typename ITR>
	inline unsigned insert(ITR &itr, unsigned max, T value) {
		unsigned failedInsert = 0;
		while (itr != itr.end()) {
			//Last iteration check if value was obliterated
			if (!insert(*itr, value, max)) {
				++failedInsert;
			}
			++itr;
		}
		return failedInsert;
	}

	/*
	 * Returns false if unable to insert hashes values
	 * Inserts hash functions in random order
	 */
	inline bool insert(const size_t *hashes, T value, unsigned max) {
		unsigned count = 0;
		std::vector<unsigned> hashOrder;
		bool saturated = true;
		//for random number generator seed
		size_t randValue = value;

		//check values and if value set
		for (unsigned i = 0; i < m_hashNum; ++i) {
			//check if values are already set
			size_t pos = m_rankSupport(hashes[i] % m_bv.size());
			//check for saturation
			T oldVal = m_data[pos];
			if (oldVal > s_mask) {
				oldVal = oldVal & s_antiMask;
			} else {
				saturated = false;
			}
			if (oldVal == value) {
				++count;
			}
			else{
				hashOrder.push_back(i);
			}
			if (count >= max) {
				return true;
			}
			randValue ^= hashes[i];
		}
		std::minstd_rand g(randValue);
		std::shuffle(hashOrder.begin(), hashOrder.end(), g);

		//insert seeds in random order
		for (std::vector<unsigned>::iterator itr = hashOrder.begin();
				itr != hashOrder.end(); ++itr) {
			size_t pos = m_rankSupport(hashes[*itr] % m_bv.size());
			//check for saturation
			T oldVal = setVal(&m_data[pos], value);
			if (oldVal > s_mask) {
				oldVal = oldVal & s_antiMask;
			} else {
				saturated = false;
			}
			if (oldVal == 0) {
				++count;
			}
			if (count >= max) {
				return true;
			}
		}
		if (count == 0) {
			if (!saturated) {
				assert(max == 1); //if this triggers then spaced seed is probably not symmetric
				saturate(hashes);
			}
			return false;
		}
		return true;
	}

	inline void saturate(const size_t *hashes) {
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t pos = m_rankSupport(hashes[i] % m_bv.size());
			__sync_or_and_fetch(&m_data[pos], s_mask);
		}
	}

	/*
	 * Return the position of a hash values
	 */
	inline vector<size_t> atPos(const size_t *hashes, unsigned &finished, unsigned maxMiss = 0){
		vector<size_t> rankPos(m_hashNum);
		unsigned misses = 0;
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_bv.size();
			if (m_bv[pos]) {
				rankPos[i] = m_rankSupport(pos);
			} else {
				++misses;
				if (misses > maxMiss) {
					finished = i;
					return rankPos;
				}
			}
		}
		finished = m_hashNum;
		return rankPos;
	}

	/*
	 * Returns if IDs are saturated
	 * ~2 cache misses
	 */
	inline bool isSaturated(const size_t *hashes) {
		unsigned misses = 0;
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_bv.size();
			if(m_data[m_rankSupport(pos)] < s_mask){
				return false;
			}
		}
		return true;
	}

	/*
	 * No saturation masking
	 */
	inline vector<T> at(const size_t *hashes, unsigned maxMiss = 0) {
		vector<T> results(m_hashNum);
		vector<size_t> rankPos(m_hashNum);
		unsigned misses = 0;
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_bv.size();
			if (m_bv[pos] == 0) {
				++misses;
				if (misses > maxMiss) {
					return vector<T>();
				}
			} else {
				rankPos[i] = m_rankSupport(pos);
			}
		}
		for (unsigned i = 0; i < m_hashNum; ++i){
			results[i] = m_data[rankPos[i]];
		}
		return results;
	}

	/*
	 * No saturation masking
	 */
	inline vector<size_t> getRankPos(const size_t *hashes) const{
		vector<size_t> rankPos(m_hashNum);
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_bv.size();
			if (m_bv[pos] == 0) {
			} else {
				rankPos[i] = m_rankSupport(pos);
			}
		}
		return rankPos;
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
	 * computes id frequency based on datavector
	 */
	inline void getIDCounts(vector<size_t> &counts) {
		for (size_t i = 0; i < m_dSize; ++i) {
			++counts[m_data[i] & s_antiMask];
		}
	}

	/*
	 * Returns max ID seen in filter
	 * Expensive operation.
	 */
	inline T getMaxID() const {
		T max = 0;
		for (size_t i = 0; i < m_dSize; ++i) {
			T tempVal = m_data[i] & s_antiMask;
			if(max < tempVal){
				max = tempVal;
			}
		}
		return max;
	}

	inline size_t getPop() const {
		size_t index = m_bv.size() - 1;
		while (m_bv[index] == 0) {
			--index;
		}
		return m_rankSupport(index) + 1;
	}

	inline size_t getPopNonZero() const {
		size_t count = 0;
		for (size_t i = 0; i < m_dSize; ++i) {
			if (m_data[i] != 0) {
				++count;
			}
		}
		return count;
	}

	/*
	 * Checks data array for abnormal IDs
	 * (i.e. values greater than what is possible)
	 * Returns first abnormal ID or value of maxVal if no abnormal IDs are found
	 * For debugging
	 */
	inline ID checkValues(T maxVal) const {
		for (size_t i = 0; i < m_dSize; ++i) {
			if ((m_data[i] & s_antiMask) > maxVal) {
				return m_data[i];
			}
		}
		return maxVal;
	}

	inline size_t getPopSaturated() const {
		size_t count = 0;
		for (size_t i = 0; i < m_dSize; ++i) {
			if (m_data[i] > s_mask) {
				++count;
			}
		}
		return count;
	}

	inline size_t size() {
		return m_bv.size();
	}

	//overwrites existing value
	inline void setData(size_t pos, T id){
		m_data[pos] = id;
	}

	inline vector<T> getData(const vector<size_t> &rankPos) const{
		vector<T> results(rankPos.size());
		for (unsigned i = 0; i < m_hashNum; ++i) {
			results[i] = m_data[rankPos[i]];
		}
		return results;
	}

	~MIBloomFilter() {
		delete[] m_data;
	}

private:
	// Driver function to sort the vector elements
	// by second element of pairs
	static bool sortbysec(const pair<int,int> &a,
	              const pair<int,int> &b)
	{
	    return (a.second < b.second);
	}

	/*
	 * Helper function for header storage
	 */
	inline void writeHeader(ofstream &out) const {
		FileHeader header;
		char magic[9];
		strncpy(magic, MAGICSTR, 8);
		magic[8] = '\0';
		strncpy(header.magic, magic, 8);

		header.hlen = sizeof(struct FileHeader) + m_kmerSize * m_sseeds.size();
		header.kmer = m_kmerSize;
		header.size = m_dSize;
		header.nhash = m_hashNum;

		cerr << "Writing header... magic: " << magic << " hlen: " << header.hlen
				<< " nhash: " << header.nhash << " size: " << header.size
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
	 * Returns old value that was inside
	 */
	inline T setVal(T *val, T newVal) {
		T oldValue;
		do {
			oldValue = *val;
			if (oldValue != 0)
				break;
		} while (!__sync_bool_compare_and_swap(val, oldValue, newVal));
		return oldValue;
	}

	static inline unsigned nChoosek(unsigned n, unsigned k) {
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

	unsigned m_hashNum;
	unsigned m_kmerSize;

	typedef vector<vector<unsigned> > SeedVal;
	vector<string> m_sseeds;

	double m_probSaturated;
	SeedVal m_ssVal;
	const char* MAGICSTR = "MIBLOOMF";
};

#endif /* MIBLOOMFILTER_HPP_ */
