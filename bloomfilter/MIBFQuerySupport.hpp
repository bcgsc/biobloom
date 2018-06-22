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
//#include <set>
#include <boost/math/distributions/binomial.hpp>
#include "stHashIterator.hpp"
#include "ntHashIterator.hpp"

using namespace std;
using boost::math::binomial;

//T = ID type, H = rolling hash itr
template<typename T>
class MIBFQuerySupport {
public:
	struct QueryResult {
		T id;
		unsigned count;
//		unsigned nonSatCount;
//		unsigned nonSatFrameCounts;
//		unsigned totalCount;
//		unsigned totalnonSatCount;
		bool strand;
	};

	MIBFQuerySupport(const MIBloomFilter<T> &miBF,
			const vector<double> &perFrameProb, unsigned extraCount,
			unsigned extraFrameLimit, unsigned maxMiss, double rateSaturated) :
			m_miBF(miBF), m_perFrameProb(perFrameProb), m_extraCount(
					extraCount), m_extraFrameLimit(extraFrameLimit), m_maxMiss(
					maxMiss), m_rateSaturated(rateSaturated), m_satCount(0), m_evalCount(
					0), m_rankPos(miBF.getHashNum()), m_hits(miBF.getHashNum()) {

		//this should be a very small array most of the time
		m_signifResults.reserve(m_perFrameProb.size());
		//this should always be a small array
		m_seenSet.reserve(miBF.getHashNum());

		m_strandCounts.resize(m_perFrameProb.size());
		m_strandCounts.min_load_factor(0.0);
		m_strandCounts.set_empty_key(0);
		m_counts.resize(m_perFrameProb.size());
		m_counts.min_load_factor(0.0);
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
	const vector<QueryResult> &queryStrandJunction(stHashIterator &itr,
			const vector<unsigned> &minCount) {
		//reset reusable values
		m_candidateMatches.clear();
		m_strandCounts.clear();
		m_signifResults.clear();
		m_satCount = 0;
		m_evalCount = 0;

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		while (itr != itr.end() && !candidateFound) {
			candidateFound = updateCountsSeedsStrand(itr, minCount, bestCount,
					secondBestCount, extraFrame);
			++itr;
		}
//		probSaturated = calcSat(totalEvaluated, m_rateSaturated,
//				saturatedCount);
		summarizeCandiatesStrand(bestCount);
		return m_signifResults;
	}

	/*
	 * Strand & region aware query
	 * Takes a read, region hash table and computes most likely hit
	 * TODO: If junction exists, return position of junction to position read
	 * TODO: use different perFrameProb (generalized to shared frames to increase sensitivity)
	 */
	const vector<QueryResult> &queryStrand(stHashIterator &itr,
			const vector<unsigned> &minCount) {
		//reset reusable values
		m_candidateMatches.clear();
		m_strandCounts.clear();
		m_signifResults.clear();
		m_satCount = 0;
		m_evalCount = 0;

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		while (itr != itr.end() && !candidateFound) {
			candidateFound = updateCountsSeedsStrand(itr, minCount, bestCount,
					secondBestCount, extraFrame);
			++itr;
		}
//		probSaturated = calcSat(totalEvaluated, m_rateSaturated,
//				saturatedCount);
		summarizeCandiatesStrand(bestCount);
		return m_signifResults;
	}

	/*
	 * Strand & region aware query
	 * Takes a read, region hash table and computes most likely hit
	 * TODO: If junction exists, return position of junction to position read
	 * TODO: use different perFrameProb (generalized to shared frames to increase sensitivity)
	 */
	const vector<QueryResult> &queryStrand(ntHashIterator &itr,
			const vector<unsigned> &minCount) {
		//reset reusable values
		m_candidateMatches.clear();
		m_strandCounts.clear();
		m_signifResults.clear();
		m_satCount = 0;
		m_evalCount = 0;

		unsigned extraFrame = 0;
		unsigned bestCount = 0;
		unsigned secondBestCount = 0;
		bool candidateFound = false;

		while (itr != itr.end() && !candidateFound) {
			candidateFound = updateCountsKmerStrand(itr, minCount, bestCount,
					secondBestCount, extraFrame);
			++itr;
		}
//		probSaturated = calcSat(evaluatedValues,
//				pow(m_rateSaturated, m_miBF.getHashNum()), saturatedCount);
		summarizeCandiatesStrand(bestCount);
		return m_signifResults;
	}

	/*
	 * totalTrials = number of possible trials that can be checked
	 */
	const vector<QueryResult> &query(stHashIterator &itr,
			const vector<unsigned> &minCount) {
		//reset reusable values
		m_candidateMatches.clear();
		m_counts.clear();
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
//		probSaturated = calcSat(evaluatedValues, m_rateSaturated,
//				saturatedCount);
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
		m_counts.clear();
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
//		probSaturated = calcSat(evaluatedValues,
//				pow(m_rateSaturated, m_miBF.getHashNum()), saturatedCount);
		summarizeCandiates(bestCount);
		return m_signifResults;
	}

	const vector<QueryResult> &query(stHashIterator &itr1, stHashIterator &itr2,
			const vector<unsigned> &minCount) {
		m_candidateMatches.clear();
		m_counts.clear();
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
//		probSaturated = calcSat(evaluatedValues, m_rateSaturated,
//				saturatedCount);
		summarizeCandiates(bestCount);
		return m_signifResults;
	}

	const vector<QueryResult> &query(ntHashIterator &itr1, ntHashIterator &itr2,
			const vector<unsigned> &minCount) {
		m_candidateMatches.clear();
		m_counts.clear();
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
//		probSaturated = calcSat(evaluatedValues,
//				pow(m_rateSaturated, m_miBF.getHashNum()), saturatedCount);
		summarizeCandiates(bestCount);
		return m_signifResults;
	}

	unsigned getSatCount() const {
		return m_satCount;
	}

	unsigned getEvalCount() const {
		return m_evalCount;
	}

private:

	//TODO: sort by counts and then p-value, needs refinement!
	static bool sortCandidates(const QueryResult &a, const QueryResult &b) {
		return (b.count > a.count);
	}

	//contains reference to parent
	const MIBloomFilter<T> &m_miBF;
	const vector<double> &m_perFrameProb;

	//not references, but shared other objects or static variables
	const unsigned m_extraCount;
	const unsigned m_extraFrameLimit;
	const unsigned m_maxMiss;
	const double m_rateSaturated;

	//resusable variables
	unsigned m_satCount;
	unsigned m_evalCount;

	//resusable objects
	vector<size_t> m_rankPos;
	vector<bool> m_hits;
	vector<QueryResult> m_signifResults;
	google::dense_hash_map<T, pair<unsigned, unsigned>> m_strandCounts;
	google::dense_hash_map<T, unsigned> m_counts;
//	google::dense_hash_map<T, unsigned> m_nonSatCounts;
//	google::dense_hash_map<T, unsigned> m_nonSatFrameCounts;
//	google::dense_hash_map<T, unsigned> m_totalCounts;
//	google::dense_hash_map<T, unsigned> m_totalnonSatCounts;
//	set<T> m_candidateMatches; //should be a small set
//	set<T> m_seenSet; //always a small set
	vector<T> m_candidateMatches;
	vector<T> m_seenSet;

	bool updateCountsSeedsStrand(const stHashIterator &itr,
			const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame) {
		bool candidateFound = false;
		unsigned misses = m_miBF.atRank(*itr, m_rankPos, m_hits,
				m_maxMiss);
		if (misses <= m_maxMiss) {
			m_seenSet.clear();
			for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
				if (m_hits[i]) {
					T result = m_miBF.getData(m_rankPos[i]);
					//check for saturation
					++m_evalCount;
					if (result > m_miBF.s_mask) {
						++m_satCount;
						result &= m_miBF.s_antiMask;
					}
					//derive strand
					bool curStrand = result > m_miBF.s_strand;
					bool strandAgree = curStrand == itr.strandArray()[i];
					if (curStrand) {
						result &= m_miBF.s_antiStrand;
					}
					if (find(m_seenSet.begin(), m_seenSet.end(), result) == m_seenSet.end()) {
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
								if(find(m_candidateMatches.begin(), m_candidateMatches.end(), result) == m_candidateMatches.end()){
									m_candidateMatches.push_back(result);
								}
							}
						}
						m_seenSet.push_back(result);
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

	bool updateCountsSeeds(const stHashIterator &itr,
			const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame) {
		bool candidateFound = false;
		unsigned misses = m_miBF.atRank(*itr, m_rankPos, m_hits,
				m_maxMiss);
		if (misses <= m_maxMiss) {
			m_seenSet.clear();
			for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
				if (m_hits[i]) {
					T result = m_miBF.getData(m_rankPos[i]);
					++m_evalCount;
					//check for saturation
					if (result > m_miBF.s_mask) {
						result &= m_miBF.s_antiMask;
						++m_satCount;
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

	bool updateCountsKmerStrand(const ntHashIterator &itr,
			const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame) {
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
				//derive strand
				bool curStrand = result > m_miBF.s_strand;
				bool strandAgree = curStrand == itr.strandArray()[i];
				if (curStrand) {
					result &= m_miBF.s_antiStrand;
				}
				if (find(m_seenSet.begin(), m_seenSet.end(), result)
						== m_seenSet.end()) {
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
							if (find(m_candidateMatches.begin(),
									m_candidateMatches.end(), result)
									== m_candidateMatches.end()) {
								m_candidateMatches.push_back(result);
							}
						}
					}
					m_seenSet.push_back(result);
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
				++m_satCount;
			}
		}
		++m_evalCount;
		return candidateFound;
	}


	bool updateCountsKmer(const ntHashIterator &itr,
			const vector<unsigned> &minCount, unsigned &bestCount,
			unsigned &secondBestCount, unsigned &extraFrame) {
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
				++m_satCount;
			}
		}
		++m_evalCount;
		return candidateFound;
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
	}

	void summarizeCandiatesStrand(unsigned bestCount) {
		if (m_candidateMatches.size()) {
			for (typename vector<T>::const_iterator candidate =
					m_candidateMatches.begin();
					candidate != m_candidateMatches.end(); candidate++) {
				//true = rv, false = fw
				bool strand = m_strandCounts[*candidate].first
						> m_strandCounts[*candidate].second;
				unsigned tempCount =
						strand ?
								m_strandCounts[*candidate].first :
								m_strandCounts[*candidate].second;
				if (bestCount <= tempCount + m_extraCount) {
					QueryResult result;
					result.id = *candidate;
					result.strand = strand;
					result.count = tempCount;
					m_signifResults.push_back(result);
				}
			}
			sort(m_signifResults.begin(), m_signifResults.end(),
					sortCandidates);
		}
	}
};

#endif /* CHROMIUMMAP_MIBFQUERYSUPPORT_HPP_ */
