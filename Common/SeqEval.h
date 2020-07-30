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
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"
#include <boost/math/distributions/binomial.hpp>
#include "Common/SDust.hpp"

using namespace std;

namespace SeqEval {

inline double denormalizeScore(double score, unsigned kmerSize, size_t seqLen) {
	assert(score >= 0 && score <= 1);
	return score * (seqLen - kmerSize + 1);
}

inline double normalizeScore(double score, unsigned kmerSize, size_t seqLen) {
	return score / (seqLen - kmerSize + 1);
}

inline bool evalSimple(const string &rec, const BloomFilter &filter,
		double threshold, const BloomFilter *subtract =
		NULL, SDust *sduster = NULL) {

	const double thres = denormalizeScore(threshold, filter.getKmerSize(),
			rec.length());
	const double antiThres = floor(
			denormalizeScore(1.0 - threshold, filter.getKmerSize(),
					rec.length()));

	double score = 0;
	unsigned antiScore = 0;
	unsigned streak = 0;
	ntHashIterator itr(rec, filter.getKmerSize(), filter.getKmerSize());
	unsigned prevPos = 0;
	if (itr != itr.end()) {
		if (!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr)) {
			if (subtract == NULL || !subtract->contains(*itr))
				score += 0.5;
			if (thres <= score) {
				return true;
			}
			++streak;
		} else {
			if (antiThres <= ++antiScore)
				return false;
		}
		prevPos = itr.pos();
		++itr;
	}
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
		if (!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr)) {
			if (streak == 0) {
				if (subtract == NULL || !subtract->contains(*itr))
					score += 0.5;
			} else {
				if (subtract == NULL || !subtract->contains(*itr))
					++score;
			}
			if (thres <= score) {
				return true;
			}
			prevPos = itr.pos();
			++itr;
			++streak;
		} else {
			if (streak < opt::streakThreshold) {
				if (antiThres <= ++antiScore)
					return false;
				prevPos = itr.pos();
				++itr;
			} else {
				//TODO: Determine if faster to force re-init
				unsigned skipEnd = itr.pos() + filter.getKmerSize();
				//skip lookups
				while (itr.pos() < skipEnd) {
					if (antiThres <= ++antiScore)
						return false;
					prevPos = itr.pos();
					++itr;
				}
			}
			streak = 0;
		}
	}
	return false;
}

inline bool evalHarmonic(const string &rec, const BloomFilter &filter,
		double threshold, const BloomFilter *subtract =
		NULL, SDust *sduster = NULL) {

	const double thres = denormalizeScore(threshold, filter.getKmerSize(),
			rec.length());
	const double antiThres = floor(
			denormalizeScore(1.0 - threshold, filter.getKmerSize(),
					rec.length()));

	double score = 0;
	unsigned antiScore = 0;
	unsigned streak = 0;
	ntHashIterator itr(rec, filter.getKmerSize(), filter.getKmerSize());
	unsigned prevPos = 0;
	if (itr != itr.end()) {
		if (!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr)) {
			if (subtract == NULL || !subtract->contains(*itr))
				score += 0.5;
			if (thres <= score) {
				return true;
			}
			++streak;
		} else {
			if (antiThres <= ++antiScore)
				return false;
		}
		prevPos = itr.pos();
		++itr;
	}
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
		if (!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr)) {
			if (streak == 0) {
				if (subtract == NULL || !subtract->contains(*itr))
					score += 0.5;
			} else {
				if (subtract == NULL || !subtract->contains(*itr))
					score += 1.0 - 1.0 / (1.0 + double(streak));
			}
			if (thres <= score) {
				return true;
			}
			prevPos = itr.pos();
			++itr;
			++streak;
		} else {
			if (streak < opt::streakThreshold) {
				if (antiThres <= ++antiScore)
					return false;
				prevPos = itr.pos();
				++itr;
			} else {
				//TODO: Determine if faster to force re-init
				unsigned skipEnd = itr.pos() + filter.getKmerSize();
				//skip lookups
				while (itr.pos() < skipEnd) {
					if (antiThres <= ++antiScore)
						return false;
					prevPos = itr.pos();
					++itr;
				}
			}
			streak = 0;
		}
	}
	return false;
}

/*
 * Calculate minimum number of matches to classify sequences, accounting
 * for false positive rate
 */
inline size_t calcMinCount(size_t frameLen, double bfFPR, double minFPR) {
	using namespace boost::math::policies;
	using namespace boost::math;
	typedef boost::math::binomial_distribution<double,
			policy<discrete_quantile<integer_round_up> > > binom_round_up;
	binom_round_up bin(frameLen, bfFPR);
	unsigned i = quantile(complement(bin, minFPR));
	return i > 1 ? i : 1;
}

inline double calcProbMatches(size_t frameLen, double bfFPR, size_t matches) {
	using namespace boost::math::policies;
	using namespace boost::math;
	typedef boost::math::binomial_distribution<double,
			policy<discrete_quantile<integer_round_up> > > binom_round_up;
	binom_round_up bin(frameLen, bfFPR);
	return cdf(complement(bin, matches));
}

inline bool evalBinomial(const string &rec, const BloomFilter &filter,
		double threshold, const BloomFilter *subtract =
		NULL, SDust *sduster = NULL) {
	if (rec.size() < filter.getKmerSize()) {
		return false;
	}
	const unsigned frameLen = rec.size() - filter.getKmerSize() + 1;
	const unsigned thres = calcMinCount(frameLen, filter.getFPRPrecompute(), threshold);
//	cerr << filter.getFilterSize() << " " << threshold << " " << rec.size() << " " << thres << " " << frameLen << endl;
	const unsigned antiThres = frameLen - thres;

	unsigned score = 0;
	unsigned antiScore = 0;
	unsigned streak = 0;
	ntHashIterator itr(rec, filter.getKmerSize(), filter.getKmerSize());
	unsigned prevPos = 0;
	if (itr != itr.end()) {
		if (!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr)) {
			if (subtract == NULL || !subtract->contains(*itr))
				score++;
			if (thres <= score) {
				return true;
			}
			++streak;
		} else {
			if (antiThres <= ++antiScore)
				return false;
		}
		prevPos = itr.pos();
		++itr;
	}
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
		if (!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr)) {
			if (subtract == NULL || !subtract->contains(*itr))
					++score;
			if (thres <= score) {
				return true;
			}
			prevPos = itr.pos();
			++itr;
			++streak;
		} else {
			if (streak < opt::streakThreshold) {
				if (antiThres <= ++antiScore)
					return false;
				prevPos = itr.pos();
				++itr;
			} else {
				unsigned skipEnd = itr.pos() + filter.getKmerSize();
				//skip lookups
				while (itr.pos() < skipEnd) {
					if (antiThres <= ++antiScore)
						return false;
					prevPos = itr.pos();
					++itr;
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
		unsigned minMatchLen, const BloomFilter *subtract = NULL,
		SDust *sduster = NULL) {
	// number of contiguous k-mers matched
	unsigned matchLen = 0;
	size_t l = rec.length();

	ntHashIterator itr(rec, filter.getHashNum(), filter.getKmerSize());
	unsigned prevPos = 0;
	while (itr != itr.end()) {
		// quit early if there is no hope
		if (l - itr.pos() + matchLen < minMatchLen)
			return false;

		//check if k-mer has deviated/started again
		//TODO try to terminate before itr has to re-init after skipping
		if (itr.pos() != prevPos + 1) {
			matchLen = 0;
		}
		if (!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr)) {
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
		prevPos = itr.pos();
		++itr;
	}
	return false;
}

inline double evalSimpleScore(const string &rec, const BloomFilter &filter,
		const BloomFilter *subtract = NULL, SDust *sduster = NULL) {

	double score = 0;
	unsigned streak = 0;
	ntHashIterator itr(rec, filter.getKmerSize(), filter.getKmerSize());
	unsigned prevPos = 0;
	while (itr != itr.end()) {
		//check if k-mer has deviated/started again
		//TODO try to terminate before itr has to re-init after skipping
		if (itr.pos() != prevPos + 1) {
			streak = 0;
		}
		if (!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr)) {
			if (streak == 0) {
				if (subtract == NULL || !subtract->contains(*itr))
					score += 0.5;
			} else {
				if (subtract == NULL || !subtract->contains(*itr))
					++score;
			}
			prevPos = itr.pos();
			++itr;
			++streak;
		} else {
			if (streak < opt::streakThreshold) {
				prevPos = itr.pos();
				++itr;
			} else {
				//TODO: Determine if faster to force re-init
				unsigned skipEnd = itr.pos() + filter.getKmerSize();
				//skip lookups
				while (itr.pos() < skipEnd) {
					prevPos = itr.pos();
					++itr;
				}
			}
			streak = 0;
		}
	}
	return normalizeScore(score, filter.getKmerSize(), rec.length());
}

inline double evalHarmonicScore(const string &rec, const BloomFilter &filter,
		const BloomFilter *subtract = NULL, SDust *sduster = NULL) {
	double score = 0;
	unsigned streak = 0;
	ntHashIterator itr(rec, filter.getKmerSize(), filter.getKmerSize());
	unsigned prevPos = 0;
	while (itr != itr.end()) {
		//check if k-mer has deviated/started again
		//TODO try to terminate before itr has to re-init after skipping
		if (itr.pos() != prevPos + 1) {
			streak = 0;
		}
		if (!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr)) {
			if (streak == 0) {
				if (subtract == NULL || !subtract->contains(*itr))
					score += 0.5;
			} else {
				if (subtract == NULL || !subtract->contains(*itr))
					score += 1.0 - 1.0 / (1.0 + double(streak));
			}
			prevPos = itr.pos();
			++itr;
			++streak;
		} else {
			if (streak < opt::streakThreshold) {
				prevPos = itr.pos();
				++itr;
			} else {
				//TODO: Determine if faster to force re-init
				unsigned skipEnd = itr.pos() + filter.getKmerSize();
				//skip lookups
				while (itr.pos() < skipEnd) {
					prevPos = itr.pos();
					++itr;
				}
			}
			streak = 0;
		}
	}
	return normalizeScore(score, filter.getKmerSize(), rec.length());
}

inline unsigned evalMinMatchLenScore(const string &rec,
		const BloomFilter &filter, const BloomFilter *subtract = NULL,
		SDust *sduster = NULL) {
	unsigned matchLen = 0;
	ntHashIterator itr(rec, filter.getHashNum(), filter.getKmerSize());
	unsigned prevPos = 0;
	while (itr != itr.end()) {
		//check if k-mer has deviated/started again
		//TODO try to terminate before itr has to re-init after skipping
		if (itr.pos() != prevPos + 1) {
			matchLen = 0;
		}
		if ((!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr))) {
			if (subtract == NULL || !subtract->contains(*itr)) {
				if (matchLen == 0)
					matchLen = filter.getKmerSize();
				else
					++matchLen;
			}
		} else {
			matchLen = 0;
		}
		prevPos = itr.pos();
		++itr;
	}
	return matchLen;
}

inline double evalBinomialScore(const string &rec, const BloomFilter &filter,
		const BloomFilter *subtract = NULL, SDust *sduster = NULL) {
	if (rec.size() < filter.getKmerSize()) {
		return 1.0;
	}
	const unsigned frameLen = rec.size() - filter.getKmerSize() + 1;
	unsigned score = 0;
	unsigned streak = 0;
	ntHashIterator itr(rec, filter.getKmerSize(), filter.getKmerSize());
	unsigned prevPos = 0;

	while (itr != itr.end()) {
		//check if k-mer has deviated/started again
		//TODO try to terminate before itr has to re-init after skipping
		if (itr.pos() != prevPos + 1) {
			streak = 0;
		}
		if (!(sduster != NULL && sduster->isLowComp(itr.pos()))
				&& filter.contains(*itr)) {
			if (subtract == NULL || !subtract->contains(*itr))
				++score;
			prevPos = itr.pos();
			++itr;
			++streak;
		} else {
			if (streak < opt::streakThreshold) {
				prevPos = itr.pos();
				++itr;
			} else {
				unsigned skipEnd = itr.pos() + filter.getKmerSize();
				//skip lookups
				//TODO Add skip functionality to nthash (faster skip or reconstruct)
				while (itr.pos() < skipEnd) {
					prevPos = itr.pos();
					++itr;
				}
			}
			streak = 0;
		}
	}
	return calcProbMatches(frameLen, filter.getFPRPrecompute(), score);
}

inline bool evalRead(const string &rec, const BloomFilter &filter,
		double threshold, const BloomFilter *subtract = NULL) {

	if(opt::dust){
		SDust sduster(rec);
		switch (opt::scoringMethod) {
		case opt::LENGTH:
			return evalMinMatchLen(rec, filter, (unsigned) round(threshold), subtract, &sduster);
		case opt::HARMONIC:
			return evalHarmonic(rec, filter, threshold, subtract, &sduster);
		case opt::BINOMIAL:
			return evalBinomial(rec, filter, threshold, subtract, &sduster);
		case opt::SIMPLE:
		default:
			return evalSimple(rec, filter, threshold, subtract, &sduster);
		}
	}
	else{
		switch (opt::scoringMethod) {
		case opt::LENGTH:
			return evalMinMatchLen(rec, filter, (unsigned) round(threshold),
					subtract);
		case opt::HARMONIC:
			return evalHarmonic(rec, filter, threshold, subtract);
		case opt::BINOMIAL:
			return evalBinomial(rec, filter, threshold, subtract);
		case opt::SIMPLE:
		default:
			return evalSimple(rec, filter, threshold, subtract);
		}
	}
}

inline bool evalRead(const string &rec, const BloomFilter &filter,
		double threshold, const BloomFilter &subtract) {
	return evalRead(rec, filter, threshold, &subtract);
}

/*
 * Computes exhaustively (no ending early to save speed)
 */
inline double evalScore(const string &rec, const BloomFilter &filter, const BloomFilter *subtract =
		NULL) {

	if(opt::dust){
		SDust sduster(rec);
		switch (opt::scoringMethod) {
		case opt::LENGTH:
			return evalMinMatchLenScore(rec, filter, subtract, &sduster);
		case opt::HARMONIC:
			return evalHarmonicScore(rec, filter, subtract, &sduster);
		case opt::BINOMIAL:
			return log10(evalBinomialScore(rec, filter, subtract, &sduster)) * -10;
		case opt::SIMPLE:
		default:
			return evalSimpleScore(rec, filter, subtract, &sduster);
		}
	}
	else{
		switch (opt::scoringMethod) {
		case opt::LENGTH:
			return evalMinMatchLenScore(rec, filter, subtract);
		case opt::HARMONIC:
			return evalHarmonicScore(rec, filter, subtract);
		case opt::BINOMIAL:
			return log10(evalBinomialScore(rec, filter, subtract)) * -10;
		case opt::SIMPLE:
		default:
			return evalSimpleScore(rec, filter, subtract);
		}
	}
}
}
;

#endif /* SEQEVAL_H_ */
