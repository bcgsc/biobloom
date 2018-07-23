/*
 * MIBFQuerySupport.hpp
 *
 * Purpose: To provide support for complex classification
 *
 * Functions for most accurate classification and faster heuristic classification
 * are split into sections.
 *
 * Contains support objects intended to be private per thread (copied)
 *
 *
 *  Created on: Jun 6, 2018
 *      Author: justin
 */

#ifndef MIBFQUERYSUPPORT_HPP_
#define MIBFQUERYSUPPORT_HPP_

#include "MIBloomFilter.hpp"
//#include <set>
#include <boost/math/distributions/binomial.hpp>
#include "btl_bloomfilter/stHashIterator.hpp"
#include "btl_bloomfilter/ntHashIterator.hpp"

using namespace std;
using boost::math::binomial;

//T = ID type, H = rolling hash itr
template<typename T>
class MIBFQuerySupport {
public:
	MIBFQuerySupport(const MIBloomFilter<T> &miBF,
			const vector<double> &perFrameProb, unsigned extraCount,
			unsigned extraFrameLimit, unsigned maxMiss, unsigned minCount,
			bool bestHitAgree) :
			m_miBF(miBF), m_perFrameProb(perFrameProb), m_extraCount(
					extraCount), m_extraFrameLimit(extraFrameLimit), m_maxMiss(
					maxMiss), m_minCount(minCount), m_bestHitAgree(bestHitAgree), m_satCount(
					0), m_evalCount(0), m_rankPos(miBF.getHashNum()), m_hits(
					miBF.getHashNum(), true), m_counts(
					vector<CountResult>(perFrameProb.size(),
							{ 0, 0, 0, 0, 0, 0 })) {
		//this should always be a small array
		m_seenSet.reserve(miBF.getHashNum());
	}

	struct QueryResult {
		T id;
		uint16_t count;
		uint16_t nonSatCount;
		uint16_t totalCount;
		uint16_t totalNonSatCount;
		uint16_t nonSatFrameCount;
		uint16_t solidCount;
		double frameProb;
	};

	struct CountResult {
		uint16_t count;
		uint16_t nonSatCount;
		uint16_t totalCount;
		uint16_t totalNonSatCount;
		uint16_t nonSatFrameCount;
		uint16_t solidCount;
	};

	/*
	 * totalTrials = number of possible trials that can be checked
	 */
	const vector<QueryResult> &query(stHashIterator &itr,
			const vector<unsigned> &minCount) {
		//reset reusable values
		m_candidateMatches.clear();
		std::fill(m_counts.begin(), m_counts.end(), CountResult());
		m_signifResults.clear();
		m_satCount = 0;
		m_evalCount = 0;

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		while (itr != itr.end() && !candidateFound) {
			candidateFound = updateCountsSeeds(itr, minCount, bestCount,
					secondBestCount, extraFrame);
			++itr;
		}
		summarizeCandiates(bestCount);

		return m_signifResults;
	}

	/*
	 * Normal query using k-mers
	 */
	const vector<QueryResult> &query(ntHashIterator &itr,
			const vector<unsigned> &minCount) {
		//reset reusable values
		m_candidateMatches.clear();
		std::fill(m_counts.begin(), m_counts.end(), CountResult());
		m_signifResults.clear();
		m_satCount = 0;
		m_evalCount = 0;

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		while (itr != itr.end() && !candidateFound) {
			candidateFound = updateCountsKmer(itr, minCount, bestCount,
					secondBestCount, extraFrame);
			++itr;
		}
		summarizeCandiates(bestCount);
		return m_signifResults;
	}

	const vector<QueryResult> &query(stHashIterator &itr1, stHashIterator &itr2,
			const vector<unsigned> &minCount) {
		m_candidateMatches.clear();
		std::fill(m_counts.begin(), m_counts.end(), CountResult());
		m_signifResults.clear();
		m_satCount = 0;
		m_evalCount = 0;

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned frameCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		while ((itr1 != itr1.end() && itr2 != itr2.end()) && !candidateFound) {
			stHashIterator &itr = frameCount % 2 == 0 && itr1 != itr1.end() ? itr1 :
						frameCount % 2 == 1 && itr2 != itr2.end() ? itr2 : itr1;
			candidateFound = updateCountsSeeds(itr, minCount, bestCount,
					secondBestCount, extraFrame);
			++itr;
		}
		summarizeCandiates(bestCount);
		return m_signifResults;
	}

	const vector<QueryResult> &query(ntHashIterator &itr1, ntHashIterator &itr2,
			const vector<unsigned> &minCount) {
		m_candidateMatches.clear();
		std::fill(m_counts.begin(), m_counts.end(), CountResult());
		m_signifResults.clear();
		m_satCount = 0;
		m_evalCount = 0;

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned frameCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		while ((itr1 != itr1.end() && itr2 != itr2.end()) && !candidateFound) {
			ntHashIterator &itr = frameCount % 2 == 0 && itr1 != itr1.end() ? itr1 :
						frameCount % 2 == 1 && itr2 != itr2.end() ? itr2 : itr1;
			candidateFound = updateCountsKmer(itr, minCount, bestCount,
					secondBestCount, extraFrame);
			++itr;
		}
		summarizeCandiates(bestCount);
		return m_signifResults;
	}

	unsigned getSatCount() const {
		return m_satCount;
	}

	unsigned getEvalCount() const {
		return m_evalCount;
	}

	//debugging functions:

	void printAllCounts(const vector<string> &ids){
		for(size_t i = 0; i< m_counts.size(); ++i){
			if(m_counts[i].totalCount > 0)
			{
			cout << i << "\t" << ids[i] << "\t" << m_counts[i].nonSatFrameCount <<
					"\t" << m_counts[i].count <<
					"\t" << m_counts[i].solidCount <<
					"\t" << m_counts[i].nonSatCount <<
					"\t" << m_counts[i].totalNonSatCount <<
					"\t" << m_counts[i].totalCount <<
					"\n";
			}
		}
	}

	/*
	 * Debugging
	 * Computes criteria used for judging a read consisting of:
	 * Position of matches
	 * Number of actually evaluated k-mers
	 * Return count of matching k-mers to set
	 */
	//TODO saturation not handle correctly
	inline vector<unsigned> getMatchSignature(const string &seq,
			unsigned &evaluatedSeeds, vector<vector<pair<ID, bool> > > &hitsPattern) {
		vector<unsigned> matchPos;
		matchPos.reserve(seq.size() - m_miBF.getKmerSize());

		if (m_miBF.getSeedValues().size() > 0) {
			stHashIterator itr(seq, m_miBF.getSeedValues(), m_miBF.getHashNum(),
					m_miBF.getKmerSize());
			while (itr != itr.end()) {
				if (m_miBF.atRank(*itr, m_rankPos, m_hits, m_maxMiss)) {
					vector<ID> results = m_miBF.getData(m_rankPos);
					vector<pair<ID, bool>> processedResults(results.size(),
							pair<ID, bool>(0, false));
					for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
						if (m_hits[i]) {
							ID tempResult = results[i];
							if (tempResult > MIBloomFilter<ID>::s_mask) {
								processedResults[i] = pair<ID, bool>(
										tempResult
												& MIBloomFilter<ID>::s_antiMask,
										true);
							} else {
								processedResults[i] = pair<ID, bool>(
										tempResult
												& MIBloomFilter<ID>::s_antiMask,
										false);
							}
						}
					}
					matchPos.push_back(itr.pos());
					hitsPattern.push_back(processedResults);
				}
				++itr;
				++evaluatedSeeds;
			}
		} else {
			ntHashIterator itr(seq, m_miBF.getHashNum(),
					m_miBF.getKmerSize());
			while (itr != itr.end()) {
				if (m_miBF.atRank(*itr, m_rankPos)) {
					vector<ID> results = m_miBF.getData(m_rankPos);
					vector<pair<ID, bool>> processedResults(results.size(),
							pair<ID, bool>(0, false));
					if (results.size() > 0) {
						for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
							ID tempResult = results[i];
							if (tempResult > MIBloomFilter<ID>::s_mask) {
								processedResults[i] = pair<ID, bool>(
										tempResult
												& MIBloomFilter<ID>::s_antiMask,
										true);
							} else {
								processedResults[i] = pair<ID, bool>(
										tempResult
												& MIBloomFilter<ID>::s_antiMask,
										false);
							}
						}
						matchPos.push_back(itr.pos());
						hitsPattern.push_back(processedResults);
					}
				}
				++itr;
				++evaluatedSeeds;
			}
		}
		return matchPos;
	}

private:

	/*
	 * Sort in order of
	 * nonSatFrameCount
	 * count
	 * solidCount
	 * nonSatCount
	 * totalNonSatCount
	 * totalCount
	 * frameProb
	 */
	static inline bool sortCandidates(const QueryResult &a,
			const QueryResult &b) {
		return (b.nonSatFrameCount == a.nonSatFrameCount ?
				(b.count == a.count ?
				(b.solidCount == a.solidCount ?
				(b.nonSatCount == a.nonSatCount ?
				(b.totalNonSatCount == a.totalNonSatCount ?
				(b.totalCount == a.totalCount ?
				(a.frameProb > b.frameProb) :
					a.totalCount > b.totalCount) :
					a.totalNonSatCount > b.totalNonSatCount) :
					a.nonSatCount > b.nonSatCount) :
					a.solidCount > b.solidCount) :
					a.count > b.count) :
					a.nonSatFrameCount > b.nonSatFrameCount);
	}

	bool checkCountAgreement(QueryResult b, QueryResult a){
		return (
				b.nonSatFrameCount >= a.nonSatFrameCount
				&& b.count >= a.count
				&& b.solidCount >= a.solidCount
				&& b.nonSatCount >= a.nonSatCount
				&& b.totalNonSatCount >= a.totalNonSatCount
				&& b.totalCount >= a.totalCount
				);
	}

	//contains reference to parent
	const MIBloomFilter<T> &m_miBF;
	const vector<double> &m_perFrameProb;

	//not references, but shared other objects or static variables
	const unsigned m_extraCount;
	const unsigned m_extraFrameLimit;
	const unsigned m_maxMiss;
	const unsigned m_minCount;
	const bool m_bestHitAgree;
//	const double m_rateSaturated;

	//resusable variables
	unsigned m_satCount;
	unsigned m_evalCount;

	//resusable objects
	vector<size_t> m_rankPos;
	vector<bool> m_hits;
	vector<QueryResult> m_signifResults;
	vector<CountResult> m_counts;
	vector<T> m_candidateMatches;
	vector<T> m_seenSet;


	bool updateCountsSeeds(const stHashIterator &itr,
			const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame) {
		bool candidateFound = false;
		unsigned misses = m_miBF.atRank(*itr, m_rankPos, m_hits, m_maxMiss);
		if (misses <= m_maxMiss) {
			candidateFound = updatesCounts(minCount, bestCount, secondBestCount,
					extraFrame, misses);
		}
		return candidateFound;
	}

	bool updateCountsKmer(const ntHashIterator &itr,
			const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame) {
		bool candidateFound = false;
		if (m_miBF.atRank(*itr, m_rankPos)) {
			candidateFound = updatesCounts(minCount, bestCount, secondBestCount,
					extraFrame);
		}
		++m_evalCount;
		return candidateFound;
	}

	bool updatesCounts(const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame, unsigned misses = 0){
		m_seenSet.clear();
		unsigned satCount = 0;
		for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
			if (m_hits[i]) {
				T resultRaw = m_miBF.getData(m_rankPos[i]);
//				assert(resultRaw);
//				assert(minCount.size());
//				assert(bestCount == 0);
//				assert(secondBestCount == 0);
//				assert(extraFrame == 0);
//				assert(misses < m_miBF.getHashNum());
//				assert(satCount == 0);
				++m_evalCount;
				bool saturated = false;
				T result = resultRaw;
				//check for saturation
				if (result > m_miBF.s_mask) {
					result &= m_miBF.s_antiMask;
					saturated = true;
					satCount++;
				} else {
					++m_counts[result].totalNonSatCount;
				}
				++m_counts[result].totalCount;
				if (find(m_seenSet.begin(), m_seenSet.end(), resultRaw)
						== m_seenSet.end()) {
					//check for saturation
					if (saturated) {
						//if the non-saturated version has not been seen before
						if (find(m_seenSet.begin(), m_seenSet.end(), result)
								== m_seenSet.end()) {
							//check is count is exceeded
							++m_counts[result].count;
						}
					} else {
						++m_counts[result].nonSatCount;
						//check is count is exceeded
						++m_counts[result].count;
					}
					m_seenSet.push_back(resultRaw);
				}
			}
		}
		if (satCount == 0) {
			for (typename vector<T>::iterator itr = m_seenSet.begin();
					itr != m_seenSet.end(); ++itr) {
				++m_counts[*itr].nonSatFrameCount;
				if (misses == 0) {
					++m_counts[*itr].solidCount;
				}
			}
		} else {
			++m_satCount;
		}
		for (typename vector<T>::iterator itr = m_seenSet.begin();
				itr != m_seenSet.end(); ++itr) {
			T result = *itr &= m_miBF.s_antiMask;
			if (result > m_miBF.s_mask) {
				//if not saturated version already exists
				if (find(m_seenSet.begin(), m_seenSet.end(), result)
						== m_seenSet.end()) {
					continue;
				}
				result &= m_miBF.s_antiMask;
			}
			if (m_counts[result].count >= minCount[result]) {
				if (find(m_candidateMatches.begin(), m_candidateMatches.end(),
						result) == m_candidateMatches.end()) {
					m_candidateMatches.push_back(*itr);
				}
				if (m_counts[result].nonSatFrameCount > bestCount) {
					bestCount = m_counts[result].nonSatFrameCount;
				} else if (m_counts[result].nonSatFrameCount > secondBestCount) {
					secondBestCount = m_counts[result].nonSatFrameCount;
				}
			}
		}
		if (bestCount <= secondBestCount + m_extraCount) {
			extraFrame = 0;
		}
		if (bestCount > secondBestCount) {
			if (m_extraFrameLimit < extraFrame++) {
				return true;
			}
		}
		return false;
	}

	double calcSat(unsigned evaluatedValues,
			double singleEventProbSaturted, unsigned saturatedCount) {
		double probSaturated = 0;
		if (saturatedCount) {
			binomial bin(evaluatedValues, singleEventProbSaturted);
			probSaturated = cdf(bin, saturatedCount - 1);
		}
		return probSaturated;
	}

	void summarizeCandiates(unsigned bestCount) {
		if (m_candidateMatches.size() && m_minCount <= bestCount) {
			for (typename vector<T>::const_iterator candidate =
					m_candidateMatches.begin();
					candidate != m_candidateMatches.end(); candidate++) {
				CountResult resultCount = m_counts[*candidate];
				if (bestCount <= resultCount.nonSatFrameCount + m_extraCount) {
					QueryResult result;
					result.id = *candidate;
					result.count = resultCount.count;
					result.nonSatCount = resultCount.nonSatCount;
					result.totalCount = resultCount.totalCount;
					result.totalNonSatCount = resultCount.totalNonSatCount;
					result.nonSatFrameCount = resultCount.nonSatFrameCount;
					result.solidCount = resultCount.solidCount;
					result.frameProb = m_perFrameProb.at(*candidate);
					m_signifResults.push_back(result);
				}
			}
			sort(m_signifResults.begin(), m_signifResults.end(),
					sortCandidates);
			if (m_bestHitAgree && m_signifResults.size() >= 2
					&& !checkCountAgreement(m_signifResults[0],
							m_signifResults[1])) {
				m_signifResults.clear();
			}
		}
	}
};

#endif /* MIBFQUERYSUPPORT_HPP_ */
