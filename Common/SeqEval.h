/*
 * SeqEval.h
 *
 * Algorithms for bloomfilter based evaluation on sequences
 *
 *  Created on: Mar 10, 2015
 *      Author: cjustin
 *
 * Todo: try to expand and transfer methods from BioBloomClassifier
 */

#ifndef SEQEVAL_H_
#define SEQEVAL_H_

#include <string>
#include "boost/unordered/unordered_map.hpp"
#include "DataLayer/FastaReader.h"
#include "Common/Options.h"

using namespace std;
using namespace boost;

namespace SeqEval {
/*
 * Evaluation algorithm with no hashValue storage
 */
inline bool evalSingle(FastqRecord rec, unsigned kmerSize, BloomFilter filter,
		double threshold, unsigned antiThreshold)
{
	ReadsProcessor proc(kmerSize);
	size_t currentLoc = 0;
	double score = 0;
	unsigned antiScore = 0;
	unsigned streak = 0;
	while (rec.seq.length() >= currentLoc + kmerSize) {
		const unsigned char* currentKmer = proc.prepSeq(rec.seq, currentLoc);
		if (streak == 0) {
			if (currentKmer != NULL) {
				if (filter.contains(currentKmer)) {
					score += 0.5;
					++streak;
				}
				else if (antiThreshold <= ++antiScore) {
					return false;
				}
				++currentLoc;
			} else {
				if (currentLoc > kmerSize) {
					currentLoc += kmerSize + 1;
					antiScore += kmerSize + 1;
				} else {
					++antiScore;
					++currentLoc;
				}
				if (antiThreshold <= antiScore) {
					return false;
				}
			}
		} else {
			if (currentKmer != NULL) {
				if (filter.contains(currentKmer)) {
					++streak;
					score += 1 - 1 / (2 * streak);
					++currentLoc;

					if (threshold <= score) {
						return true;
					}
					continue;
				}
				else if (antiThreshold <= ++antiScore) {
					return false;
				}
			} else {
				currentLoc += kmerSize + 1;
				antiScore += kmerSize + 1;
			}
			if (streak < opt::streakThreshold) {
				++currentLoc;
			} else {
				currentLoc += kmerSize;
				antiScore += kmerSize;
			}
			if (antiThreshold <= antiScore) {
				return false;
			}
			streak = 0;
		}
	}
	return false;
}

///*
// * Evaluation algorithm with hashValue storage
// */
//inline bool evalSingle(FastqRecord rec, unsigned kmerSize, BloomFilter filter,
//		double threshold, double antiThreshold, unsigned hashNum,
//		vector<vector<size_t> > &hashValues)
//{
//	ReadsProcessor proc(kmerSize);
//	size_t currentLoc = 0;
//	double score = 0;
//	unsigned antiScore = 0;
//	unsigned streak = 0;
//	while (rec.seq.length() >= currentLoc + kmerSize) {
//		const unsigned char* currentKmer = proc.prepSeq(rec.seq, currentLoc);
//		if (streak == 0) {
//			if (currentKmer != NULL) {
//				if (filter.contains(currentKmer)) {
//					score += 0.5;
//					++streak;
//				}
//				else if (antiThreshold <= ++antiScore) {
//					return false;
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
//					return false;
//				}
//			}
//		} else {
//			if (currentKmer != NULL) {
//				if (filter.contains(currentKmer)) {
//					++streak;
//					score += 1 - 1 / (2 * streak);
//					++currentLoc;
//
//					if (threshold <= score) {
//						return true;
//					}
//					continue;
//				}
//				else if (antiThreshold <= ++antiScore) {
//					return false;
//				}
//			} else {
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
//				return false;
//			}
//			streak = 0;
//		}
//}
}
;

#endif /* SEQEVAL_H_ */
