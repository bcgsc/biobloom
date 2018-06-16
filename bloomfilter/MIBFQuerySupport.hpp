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

#ifndef CHROMIUMMAP_MIBFQUERYSUPPORT_HPP_
#define CHROMIUMMAP_MIBFQUERYSUPPORT_HPP_

#include "MIBloomFilter.hpp"
#include <set>
#include <boost/math/distributions/binomial.hpp>
#include "btl_bloomfilter/stHashIterator.hpp"
#include "btl_bloomfilter/ntHashIterator.hpp"

using namespace std;
using boost::math::binomial;

//T = ID type, H = rolling hash itr
template<typename T>
class MIBFQuerySupport {
public:
	struct QueryResult {
		T id;
		unsigned count;
		bool strand;
	};

	MIBFQuerySupport(const MIBloomFilter<T> &miBF,
			const vector<double> &perFrameProb, unsigned extraCount,
			unsigned extraFrameLimit, unsigned maxMiss) :
			m_miBF(miBF), m_perFrameProb(perFrameProb), m_extraCount(
					extraCount), m_extraFrameLimit(extraFrameLimit), m_maxMiss(
					maxMiss), m_rankPos(miBF.getHashNum()), m_hits(
					miBF.getHashNum()) {

		m_signifResults.reserve(m_perFrameProb.size());
		m_strandCounts.set_empty_key(0);
		m_counts.set_empty_key(0);
	}

	//default destructor should be fine
//	virtual ~MIBFQuerySupport();

	/*
	 * Strand & region aware query
	 * Takes a read, region hash table and computes most likely hit
	 * TODO: If junction exists, return position of junction to position read
	 * TODO: use different perFrameProb (generalized to shared frames to increase sensitivity)
	 */
	const vector<QueryResult> &queryStrand(stHashIterator &itr, const vector<unsigned> &minCount,
			double rateSaturated, double &probSaturated) {
		//reset reusable values
		m_candidateMatches.clear();
		m_strandCounts.clear();
		m_signifResults.clear();

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		unsigned saturatedCount = 0;
		unsigned evaluatedValues = 0;

		while (itr != itr.end() && !candidateFound) {
			candidateFound = updateCountsSeedsStrand(itr, minCount, bestCount,
					secondBestCount, extraFrame, evaluatedValues,
					saturatedCount);
			++itr;
		}
		probSaturated = summarizeCandiates(evaluatedValues, rateSaturated,
				saturatedCount, bestCount);
		return m_signifResults;
	}

	/*
	 * Strand & region aware query
	 * Takes a read, region hash table and computes most likely hit
	 * TODO: If junction exists, return position of junction to position read
	 * TODO: use different perFrameProb (generalized to shared frames to increase sensitivity)
	 */
	const vector<QueryResult> &queryStrand(ntHashIterator &itr, const vector<unsigned> &minCount,
			double rateSaturated, double &probSaturated) {
		//reset reusable values
		m_candidateMatches.clear();
		m_strandCounts.clear();
		m_signifResults.clear();

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		unsigned saturatedCount = 0;
		unsigned evaluatedValues = 0;

		while (itr != itr.end() && !candidateFound) {
			candidateFound = updateCountsKmerStrand(itr, minCount, bestCount,
					secondBestCount, extraFrame, evaluatedValues,
					saturatedCount);
			++itr;
		}
		probSaturated = summarizeCandiates(evaluatedValues, pow(rateSaturated,m_miBF.getHashNum()),
				saturatedCount, bestCount);
		return m_signifResults;
	}

	/*
	 * totalTrials = number of possible trials that can be checked
	 */
	const vector<QueryResult> &query(stHashIterator &itr, const vector<unsigned> &minCount,
			double rateSaturated, double &probSaturated) {
		//reset reusable values
		m_candidateMatches.clear();
		m_counts.clear();
		m_signifResults.clear();

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		unsigned saturatedCount = 0;
		unsigned evaluatedValues = 0;

		while (itr != itr.end() && !candidateFound) {
			candidateFound = updateCountsSeeds(itr, minCount, bestCount,
					secondBestCount, extraFrame, evaluatedValues,
					saturatedCount);
			++itr;
		}
		probSaturated = summarizeCandiates(evaluatedValues, rateSaturated,
				saturatedCount, bestCount);

		return m_signifResults;
	}

	/*
	 * Normal query using k-mers
	 */
	const vector<QueryResult> &query(ntHashIterator &itr, const vector<unsigned> &minCount,
			double rateSaturated, double &probSaturated) {
		//reset reusable values
		m_candidateMatches.clear();
		m_counts.clear();
		m_signifResults.clear();

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		unsigned saturatedCount = 0;
		unsigned evaluatedValues = 0;

		while (itr != itr.end() && !candidateFound) {
			candidateFound = updateCountsKmer(itr, minCount, bestCount,
					secondBestCount, extraFrame, evaluatedValues,
					saturatedCount);
			++itr;
		}
		probSaturated = summarizeCandiates(evaluatedValues, pow(rateSaturated,m_miBF.getHashNum()),
				saturatedCount, bestCount);
		return m_signifResults;
	}

	const vector<QueryResult> &query(stHashIterator &itr1, stHashIterator &itr2,
			const vector<unsigned> &minCount, double rateSaturated,
			double &probSaturated) {
		m_candidateMatches.clear();
		m_counts.clear();
		m_signifResults.clear();

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned frameCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		unsigned saturatedCount = 0;
		unsigned evaluatedValues = 0;


		while ((itr1 != itr1.end() && itr2 != itr2.end()) && !candidateFound) {
			stHashIterator &itr = frameCount % 2 == 0 && itr1 != itr1.end() ? itr1 :
						frameCount % 2 == 1 && itr2 != itr2.end() ? itr2 : itr1;
			candidateFound = updateCountsSeeds(itr, minCount, bestCount,
					secondBestCount, extraFrame, evaluatedValues,
					saturatedCount);
			++itr;
		}
		probSaturated = summarizeCandiates(evaluatedValues, rateSaturated,
				saturatedCount, bestCount);
		return m_signifResults;
	}

	const vector<QueryResult> &query(ntHashIterator &itr1, ntHashIterator &itr2,
			const vector<unsigned> &minCount, double rateSaturated,
			double &probSaturated) {
		m_candidateMatches.clear();
		m_counts.clear();
		m_signifResults.clear();

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned frameCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		unsigned saturatedCount = 0;
		unsigned evaluatedValues = 0;


		while ((itr1 != itr1.end() && itr2 != itr2.end()) && !candidateFound) {
			ntHashIterator &itr = frameCount % 2 == 0 && itr1 != itr1.end() ? itr1 :
						frameCount % 2 == 1 && itr2 != itr2.end() ? itr2 : itr1;
			candidateFound = updateCountsKmer(itr, minCount, bestCount,
					secondBestCount, extraFrame, evaluatedValues,
					saturatedCount);
			++itr;
		}
		probSaturated = summarizeCandiates(evaluatedValues, pow(rateSaturated,m_miBF.getHashNum()),
				saturatedCount, bestCount);
		return m_signifResults;
	}

private:

	//TODO: sort by counts and then p-value, needs refinement!
	static bool sortCandidates(const QueryResult &a, const QueryResult &b) {
		return (b.count > a.count);
	}

	//contains reference to MIBF object
	const MIBloomFilter<T> &m_miBF;
	const vector<double> &m_perFrameProb;

	const unsigned m_extraCount;
	const unsigned m_extraFrameLimit;
	const unsigned m_maxMiss;

	//resusable objects
	vector<size_t> m_rankPos;
	vector<bool> m_hits;

	vector<QueryResult> m_signifResults;

	google::dense_hash_map<T, pair<unsigned, unsigned>> m_strandCounts;
	google::dense_hash_map<T, unsigned> m_counts;
	set<T> m_candidateMatches; //should be a small set
	set<T> m_seenSet; //always a small set

	inline bool updateCountsSeedsStrand(const stHashIterator &itr,
			const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame,
			unsigned &evaluatedValues, unsigned &saturatedCount) {
		bool candidateFound = false;
		unsigned misses = m_miBF.atRank(*itr, m_rankPos, m_hits,
				m_maxMiss);
		if (misses <= m_maxMiss) {
			m_seenSet.clear();
			for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
				if (m_hits[i]) {
					T result = m_miBF.getData(m_rankPos[i]);
					//check for saturation
					++evaluatedValues;
					if (result > m_miBF.s_mask) {
						++saturatedCount;
						result &= m_miBF.s_antiMask;
					}
					//derive strand
					bool curStrand = result > m_miBF.s_strand;
					bool strandAgree = curStrand == itr.strandArray()[i];
					if (curStrand) {
						result &= m_miBF.s_antiStrand;
					}
					if (m_seenSet.find(result) == m_seenSet.end()) {
						typename google::dense_hash_map<T,
								pair<unsigned, unsigned>>::iterator tempItr =
								m_strandCounts.find(result);
						if (tempItr == m_strandCounts.end()) {
							m_strandCounts[result] =
									strandAgree ?
											pair<unsigned, unsigned>(1, 0) :
											pair<unsigned, unsigned>(0, 1);
						} else {
							unsigned tempCount =
									strandAgree ?
											++tempItr->second.first :
											++tempItr->second.second;
							//check is count is exceeded
							if (minCount[result] <= tempCount) {
								if (tempCount > bestCount) {
									bestCount = tempCount;
								} else if (tempCount > secondBestCount) {
									secondBestCount = tempCount;
								}
								m_candidateMatches.insert(result);
							}
						}
						m_seenSet.insert(result);
					}
				}
			}
			if (bestCount <= secondBestCount + m_extraCount) {
				extraFrame = 0;
			}
			if (bestCount && bestCount > secondBestCount) {
				if (m_extraFrameLimit < extraFrame++) {
					candidateFound = true;
				}
			}
		}
		return candidateFound;
	}

	inline bool updateCountsSeeds(const stHashIterator &itr,
			const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame,
			unsigned &evaluatedValues, unsigned &saturatedCount) {
		bool candidateFound = false;
		unsigned misses = m_miBF.atRank(*itr, m_rankPos, m_hits,
				m_maxMiss);
		if (misses <= m_maxMiss) {
			m_seenSet.clear();
			for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
				if (m_hits[i]) {
					T result = m_miBF.getData(m_rankPos[i]);
					++evaluatedValues;
					//check for saturation
					if (result > m_miBF.s_mask) {
						result &= m_miBF.s_antiMask;
						++saturatedCount;
					}
					if (m_seenSet.find(result) == m_seenSet.end()) {
						typename google::dense_hash_map<T, unsigned>::iterator tempItr =
								m_counts.find(result);
						if (tempItr == m_counts.end()) {
							m_counts[result] = 1;
						} else {
							unsigned tempCount = ++tempItr->second;
							//check is count is exceeded
							if (minCount[result] <= tempCount) {
								if (tempCount > bestCount) {
									bestCount = tempCount;
								} else if (tempCount > secondBestCount) {
									secondBestCount = tempCount;
								}
								m_candidateMatches.insert(result);
							}
						}
						m_seenSet.insert(result);
					}
				}
			}
			if (bestCount <= secondBestCount + m_extraCount) {
				extraFrame = 0;
			}
			if (bestCount && bestCount > secondBestCount) {
				if (m_extraFrameLimit < extraFrame++) {
					//TODO check if saturation not resolved
					candidateFound = true;
				}
			}
		}
		return candidateFound;
	}

	inline bool updateCountsKmerStrand(const ntHashIterator &itr,
			const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame,
			unsigned &evaluatedValues, unsigned &saturatedCount) {
		bool candidateFound = false;
		if (m_miBF.atRank(*itr, m_rankPos)) {
			unsigned tempSaturatedCount = 0;
			m_seenSet.clear();
			for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
				T result = m_miBF.getData(m_rankPos[i]);
				//check for saturation
				if (result > m_miBF.s_mask) {
					result &= m_miBF.s_antiMask;
				}
				//derive strand
				bool curStrand = result > m_miBF.s_strand;
				bool strandAgree = curStrand == itr.strandArray()[i];
				if (curStrand) {
					result &= m_miBF.s_antiStrand;
				}
				if (m_seenSet.find(result) == m_seenSet.end()) {
					typename google::dense_hash_map<T, pair<unsigned, unsigned>>::iterator tempItr =
							m_strandCounts.find(result);
					if (tempItr == m_strandCounts.end()) {
						m_strandCounts[result] =
								strandAgree ?
										pair<unsigned, unsigned>(1, 0) :
										pair<unsigned, unsigned>(0, 1);
					} else {
						unsigned tempCount =
								strandAgree ?
										++tempItr->second.first :
										++tempItr->second.second;
						//check is count is exceeded
						if (minCount[result] <= tempCount) {
							if (tempCount > bestCount) {
								bestCount = tempCount;
							} else if (tempCount > secondBestCount) {
								secondBestCount = tempCount;
							}
							m_candidateMatches.insert(result);
						}
					}
					m_seenSet.insert(result);
				}
			}
			if (bestCount <= secondBestCount + m_extraCount) {
				extraFrame = 0;
			}
			if (bestCount && bestCount > secondBestCount) {
				if (m_extraFrameLimit < extraFrame++) {
					candidateFound = true;
				}
			}
			if (tempSaturatedCount == m_miBF.getHashNum()) {
				++saturatedCount;
			}
		}
		++evaluatedValues;
		return candidateFound;
	}


	inline bool updateCountsKmer(const ntHashIterator &itr,
			const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame,
			unsigned &evaluatedValues, unsigned &saturatedCount) {
		bool candidateFound = false;
		if (m_miBF.atRank(*itr, m_rankPos)) {
			unsigned tempSaturatedCount = 0;
			m_seenSet.clear();
			for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
				T result = m_miBF.getData(m_rankPos[i]);
				//check for saturation
				if (result > m_miBF.s_mask) {
					result &= m_miBF.s_antiMask;
					++tempSaturatedCount;
				}
				if (m_seenSet.find(result) == m_seenSet.end()) {
					typename google::dense_hash_map<T, unsigned>::iterator tempItr =
							m_counts.find(result);
					if (tempItr == m_counts.end()) {
						m_counts[result] = 1;
					} else {
						unsigned tempCount = ++tempItr->second;
						//check is count is exceeded
						if (minCount[result] <= tempCount) {
							if (tempCount > bestCount) {
								bestCount = tempCount;
							} else if (tempCount > secondBestCount) {
								secondBestCount = tempCount;
							}
							m_candidateMatches.insert(result);
						}
					}
					m_seenSet.insert(result);
				}
			}
			if (bestCount <= secondBestCount + m_extraCount) {
				extraFrame = 0;
			}
			if (bestCount && bestCount > secondBestCount) {
				if (m_extraFrameLimit < extraFrame++) {
					candidateFound = true;
				}
			}
			if (tempSaturatedCount == m_miBF.getHashNum()) {
				++saturatedCount;
			}
		}
		++evaluatedValues;
		return candidateFound;
	}

	inline double summarizeCandiates(unsigned evaluatedValues,
			double singleEventProbSaturted, unsigned saturatedCount, unsigned bestCount) {
		//do a statistical test if saturation rate occurs at a rate higher than random chance
		//TODO learn how do use complement cdf?
		binomial bin(evaluatedValues, 1.0 - singleEventProbSaturted);
		double probSaturated = cdf(bin, evaluatedValues - saturatedCount);

		if (m_candidateMatches.size()) {
			for (typename set<T>::const_iterator candidate =
					m_candidateMatches.begin();
					candidate != m_candidateMatches.end(); candidate++) {
				unsigned tempCount = m_counts[*candidate];
				if (bestCount <= tempCount + m_extraCount) {
					QueryResult result;
					result.id = *candidate;
					result.count = tempCount;
					m_signifResults.push_back(result);
				}
			}
			sort(m_signifResults.begin(), m_signifResults.end(),
					sortCandidates);
		}
		return probSaturated;
	}
};

#endif /* CHROMIUMMAP_MIBFQUERYSUPPORT_HPP_ */
