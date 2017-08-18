/*
 * SeqEval.h
 *
 * Algorithms for bloomfilter based evaluation on sequences
 *
 *  Created on: Mar 10, 2015
 *      Author: cjustin
 *
 * TODO: See if intermediate hash value storage helps speed thing up at all
 */

#ifndef SEQEVAL_H_
#define SEQEVAL_H_ 1

#include <string>
#include <cmath>
#include <cassert>
#include "Common/Options.h"
#include "btl_bloomfilter/BloomFilter.hpp"
#include "btl_bloomfilter/ntHashIterator.hpp"

using namespace std;

namespace SeqEval {

enum EvalMode {
	EVAL_STANDARD, EVAL_MIN_MATCH_LEN
};

extern EvalMode evalMode;

inline double denormalizeScore(double score, unsigned kmerSize, size_t seqLen) {
	assert(score >= 0 && score <= 1);
	return score * (seqLen - kmerSize + 1);
}

inline double normalizeScore(double score, unsigned kmerSize, size_t seqLen) {
	return score / (seqLen - kmerSize + 1);
}

//TODO experiment with hashvalue storage functions

inline bool evalSingle(const string &rec, const BloomFilter &filter,
		double threshold, const BloomFilter *subtract =
		NULL) {

	const double thres = denormalizeScore(threshold, filter.getKmerSize(),
			rec.length());
	const double antiThres = floor(
			denormalizeScore(1.0 - threshold, filter.getKmerSize(),
					rec.length()));

	double score = 0;
	unsigned antiScore = 0;
	unsigned streak = 0;
	ntHashIterator itr(rec, filter.getKmerSize(), filter.getKmerSize());
	unsigned prevPos = itr.pos();
	while (itr != itr.end()) {
		//check if k-mer has deviated/started again
		//TODO try to terminate before itr has to re-init after skipping
		if (itr.pos() != prevPos + 1) {
			antiScore += itr.pos() - prevPos - 1;
			if (antiThres <= antiScore) {
				return false;
			}
			streak = 0;
		}
		if (filter.contains(*itr)) {
			if (streak == 0) {
				if (subtract == NULL || !subtract->contains(*itr))
					score += 0.5;
				if (thres <= score) {
					return true;
				}
			} else {
				if (subtract == NULL || !subtract->contains(*itr))
					++score;
			}
			if (thres <= score) {
				return true;
			}
			prevPos = itr.pos();
			itr.next();
			++streak;
		} else {
			if (streak < opt::streakThreshold) {
				if (antiThres <= ++antiScore)
					return false;
				prevPos = itr.pos();
				itr.next();
			} else {
				//TODO: Determine if faster to force re-init
				unsigned skipEnd = itr.pos() + filter.getKmerSize();
				//skip lookups
				while (itr.pos() < skipEnd) {
					if (antiThres <= ++antiScore)
						return false;
					prevPos = itr.pos();
					itr.next();
				}
			}
			streak = 0;
		}
	}
	return false;
}

/*
 * Evaluation algorithm based on minimum number of contiguous matching bases.
 */
inline bool evalMinMatchLen(const string &rec, const BloomFilter &filter,
		unsigned minMatchLen, const BloomFilter *subtract = NULL) {
	// number of contiguous k-mers matched
	unsigned matchLen = 0;
	size_t l = rec.length();

	ntHashIterator itr(rec, filter.getHashNum(), filter.getKmerSize());
	unsigned prevPos = itr.pos();
	while (itr != itr.end()) {
		// quit early if there is no hope
		if (l - itr.pos() + matchLen < minMatchLen)
			return false;

		//check if k-mer has deviated/started again
		//TODO try to terminate before itr has to re-init after skipping
		if (itr.pos() != prevPos + 1) {
			matchLen = 0;
		}
		if (filter.contains(*itr)) {
			if (subtract == NULL || !subtract->contains(*itr)) {
				if (matchLen == 0)
					matchLen = filter.getKmerSize();
				else
					++matchLen;
			}
		} else {
			matchLen = 0;
		}
		// if min match length reached
		if (matchLen >= minMatchLen)
			return true;
	}
	return false;
}

inline bool evalRead(const string &rec, const BloomFilter &filter,
		double threshold, const BloomFilter *subtract = NULL) {
	switch (evalMode) {
	case EVAL_MIN_MATCH_LEN:
		return evalMinMatchLen(rec, filter, (unsigned) round(threshold),
				subtract);
	case EVAL_STANDARD:
	default:
		return evalSingle(rec, filter, threshold, subtract);
	}
	assert(false);
}

inline bool evalRead(const string &rec, const BloomFilter &filter,
		double threshold, const BloomFilter &subtract) {
	return evalRead(rec, filter, threshold, &subtract);
}

/*
 * Evaluation algorithm with no hashValue storage (optimize speed for single queries)
 * Returns score instead if just a match
 */
inline double evalSingleScore(const string &rec, const BloomFilter &filter) {

	double score = 0;
	unsigned antiScore = 0;
	unsigned streak = 0;
	ntHashIterator itr(rec, filter.getHashNum(), filter.getKmerSize());
	unsigned prevPos = itr.pos();
	while (itr != itr.end()) {
		//check if k-mer has deviated/started again
		//TODO try to terminate before itr has to re-init after skipping
		if (itr.pos() != prevPos + 1) {
			antiScore += itr.pos() - prevPos - 1;
			streak = 0;
		}
		if (filter.contains(*itr)) {
			if (streak == 0) {
				score += 0.5;
			} else {
				++score;
			}
			prevPos = itr.pos();
			itr.next();
			++streak;
		} else {
			if (streak < opt::streakThreshold) {
				streak = 0;
				prevPos = itr.pos();
				itr.next();
			} else {
				unsigned skipEnd = itr.pos() + filter.getKmerSize();
				//skip lookups
				while (itr.pos() < skipEnd) {
					prevPos = itr.pos();
					itr.next();
				}
			}
			streak = 0;
		}
	}
	return normalizeScore(score, filter.getKmerSize(), rec.length());
}

/*
 * Evaluation algorithm with no hashValue storage (optimize speed for single queries)
 * Returns score instead if just a match
 */
inline double evalSingleScore(const string &rec, const BloomFilter &filter,
		double threshold, const BloomFilter *subtract =
		NULL) {

	const double antiThres = floor(
			denormalizeScore(1.0 - threshold, filter.getKmerSize(),
					rec.length()));

	double score = 0;
	unsigned antiScore = 0;
	unsigned streak = 0;
	ntHashIterator itr(rec, filter.getKmerSize(), filter.getKmerSize());
	unsigned prevPos = itr.pos();
	while (itr != itr.end()) {
		//check if k-mer has deviated/started again
		//TODO try to terminate before itr has to re-init after skipping
		if (itr.pos() != prevPos + 1) {
			antiScore += itr.pos() - prevPos - 1;
			if (antiThres <= antiScore) {
				return score;
			}
			streak = 0;
		}
		if (filter.contains(*itr)) {
			if (streak == 0) {
				if (subtract == NULL || !subtract->contains(*itr))
					score += 0.5;
			} else {
				if (subtract == NULL || !subtract->contains(*itr))
					++score;
			}
			prevPos = itr.pos();
			itr.next();
			++streak;
		} else {
			if (streak < opt::streakThreshold) {
				if (antiThres <= ++antiScore)
					return score;
				prevPos = itr.pos();
				itr.next();
			} else {
				//TODO: Determine if faster to force re-init
				unsigned skipEnd = itr.pos() + filter.getKmerSize();
				//skip lookups
				while (itr.pos() < skipEnd) {
					if (antiThres <= ++antiScore)
						return score;
					prevPos = itr.pos();
					itr.next();
				}
			}
			streak = 0;
		}
	}
	return normalizeScore(score, filter.getKmerSize(), rec.length());
}

///*
// * Core evaluation algorithm, with ability start evaluating sequence midway
// * Evaluation algorithm with hashValue storage (minimize redundant work)
// * Also stores if position has already been visited to minimize work
// * Takes in last position visited and score and updates them accordingly
// */
//inline bool eval(const string &rec, unsigned kmerSize,
//		const BloomFilter &filter, double threshold, double antiThreshold,
//		vector<bool> &visited, unsigned &currentLoc, double &score) {
//	threshold = denormalizeScore(threshold, kmerSize, rec.length());
//	antiThreshold = denormalizeScore(antiThreshold, kmerSize, rec.length());
//	score = denormalizeScore(score, kmerSize, rec.length());
//
//	unsigned antiScore = 0;
//	unsigned streak = 0;
//	bool hit = false;
//
//	while (rec.length() >= currentLoc + kmerSize) {
//
//		//prepare hash values for filter
//
//		//check if hash value is already generated
//		if (hashValues[currentLoc].size() == 0) {
//			if (!visited[currentLoc]) {
//				const unsigned char* currentSeq = proc.prepSeq(rec, currentLoc);
//				if (currentSeq != NULL) {
//					hashValues[currentLoc] = multiHash(currentSeq,
//							filter.getHashNum(), kmerSize);
//				}
//				visited[currentLoc] = true;
//			}
//		}
//
//		if (streak == 0) {
//			if (hashValues[currentLoc].size() > 0) {
//				if (filter.contains(hashValues[currentLoc])) {
//					score += 0.5;
//					++streak;
//					if (threshold <= score) {
//						++currentLoc;
//						hit = true;
//						break;
//					}
//				} else if (antiThreshold <= ++antiScore) {
//					++currentLoc;
//					hit = false;
//					break;
//				}
//				++currentLoc;
//			} else {
//				if (currentLoc > kmerSize) {
//					currentLoc += kmerSize + 1;
//					antiScore += kmerSize + 1;
//				} else {
//					++antiScore;
//					++currentLoc;
//				}
//				if (antiThreshold <= antiScore) {
//					hit = false;
//					break;
//				}
//			}
//		} else {
//			if (hashValues[currentLoc].size() > 0) {
//				if (filter.contains(hashValues[currentLoc])) {
//					++streak;
//					++score;
//					++currentLoc;
//
//					if (threshold <= score) {
//						hit = true;
//						break;
//					}
//					continue;
//				} else if (antiThreshold <= ++antiScore) {
//					++currentLoc;
//					hit = false;
//					break;
//				}
//			} else {
//				//if has non atcg character
//				currentLoc += kmerSize + 1;
//				antiScore += kmerSize + 1;
//			}
//			if (streak < opt::streakThreshold) {
//				++currentLoc;
//			} else {
//				currentLoc += kmerSize;
//				antiScore += kmerSize;
//			}
//			if (antiThreshold <= antiScore) {
//				hit = false;
//				break;
//			}
//			streak = 0;
//		}
//	}
//	score = normalizeScore(score, kmerSize, rec.length());
//	return hit;
//}
}
;

#endif /* SEQEVAL_H_ */
