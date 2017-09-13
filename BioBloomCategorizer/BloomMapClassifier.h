/*
 * BloomMapClassifier.h
 *
 *  Created on: Mar 23, 2016
 *      Author: cjustin
 */

#ifndef BLOOMMAPCLASSIFIER_H_
#define BLOOMMAPCLASSIFIER_H_

#include <string>
#include "boost/shared_ptr.hpp"
#include "Common/Options.h"
#include "Options.h"
#include <google/dense_hash_map>
#include <google/dense_hash_set>
#include "Common/Dynamicofstream.h"
#include <vector>
#include "BioBloomClassifier.h"
#include "bloomfilter/MIBloomFilter.hpp"
#include "bloomfilter/RollingHashIterator.h"
#include "btl_bloomfilter/ntHashIterator.hpp"

using namespace std;

class BloomMapClassifier {
public:
	explicit BloomMapClassifier(const string &filterFile);
	void filter(const vector<string> &inputFiles);

	void filterPair(const string &file1, const string &file2);

	virtual ~BloomMapClassifier();
private:
	MIBloomFilter<ID> m_filter;
	vector<string> m_fullIDs;
	vector<boost::shared_ptr<vector<ID> > > m_colliIDs;

	bool fexists(const string &filename) const {
		ifstream ifile(filename.c_str());
		return ifile.good();
	}

	/*
	 * Returns the number of hits
	 */
	inline unsigned evaluateRead(const string &seq,
			google::dense_hash_map<ID, unsigned> &hitCounts) {
		unsigned nonZeroCount = 0;
		if (m_filter.getSeedValues().empty()) {
			if (m_filter.getType() == MIBloomFilter<ID>::MIBFMVAL) {
				for (ntHashIterator itr(seq, m_filter.getHashNum(),
						m_filter.getKmerSize()); itr != itr.end(); ++itr) {
					unsigned misses = 0;
					ID id = m_filter.at(*itr, opt::allowMisses, misses);
					if (id != 0) {
						if (id != opt::COLLI) {
							if (hitCounts.find(id) != hitCounts.end()) {
								++hitCounts[id];
							} else {
								hitCounts[id] = 1;
							}
						}
						++nonZeroCount;
					}
				}
			} else {
				for (ntHashIterator itr(seq, m_filter.getHashNum(),
						m_filter.getKmerSize()); itr != itr.end(); ++itr) {
					unsigned misses = 0;
					vector<ID> ids = m_filter.at(*itr, m_colliIDs,
							opt::allowMisses, misses);
					nonZeroCount += ids.size() > 0;
					for (unsigned i = 0; i < ids.size(); ++i) {
						ID id = ids[i];
						if (hitCounts.find(id) == hitCounts.end()) {
							hitCounts[id] = 0;
						}
						hitCounts[id] += ids.size() > 0;
					}
				}
			}
		} else {
			RollingHashIterator itr(seq, m_filter.getKmerSize(),
					m_filter.getSeedValues());
			if (m_filter.getType() == MIBloomFilter<ID>::MIBFMVAL) {
				while (itr != itr.end()) {
					unsigned misses = 0;
					ID id = m_filter.at(*itr, opt::allowMisses, misses);
					if (id != 0) {
						if (id != opt::COLLI) {
							if (hitCounts.find(id) == hitCounts.end()) {
								hitCounts[id] = 0;
							}
							++hitCounts[id];
							hitCounts[id] += misses == 0;
						}
						nonZeroCount += misses == 0;
						++nonZeroCount;
					}
					++itr;
				}
				nonZeroCount /= 2;
			} else {
				while (itr != itr.end()) {
					unsigned misses = 0;
					vector<ID> ids = m_filter.at(*itr, m_colliIDs,
							opt::allowMisses, misses);
					nonZeroCount += misses == 0;
					nonZeroCount += ids.size() > 0;
					for (unsigned i = 0; i < ids.size(); ++i) {
						ID id = ids[i];
						if (hitCounts.find(id) == hitCounts.end()) {
							hitCounts[id] = 0;
						}
						++hitCounts[id];
						hitCounts[id] += misses == 0;
					}
					++itr;
				}
				nonZeroCount /= 2;
			}
		}
		return nonZeroCount;
	}

	/*
	 * Returns a vector of hits to a specific ID
	 */
	void convertToHits(const google::dense_hash_map<ID, unsigned> &hitCounts,
			vector<ID> &hits) {
		unsigned bestHit = 0;
		for (google::dense_hash_map<ID, unsigned>::const_iterator i =
				hitCounts.begin(); i != hitCounts.end(); ++i) {
			if (bestHit < i->second) {
				bestHit = i->second;
			}
		}
		for (google::dense_hash_map<ID, unsigned>::const_iterator i =
				hitCounts.begin(); i != hitCounts.end(); ++i) {
			if (bestHit <= i->second + opt::delta) {
				hits.push_back(i->first);
			}
		}
	}

	/*
	 * Returns a vector of best hits to a specific read
	 * Both reads must match
	 */
	unsigned convertToHitsBoth(
			const google::dense_hash_map<ID, unsigned> &hitCounts1,
			const google::dense_hash_map<ID, unsigned> &hitCounts2,
			vector<ID> &hits) {
		unsigned bestHit = 0;
		if (m_filter.getType() == MIBloomFilter<ID>::MIBFMVAL) {
			ID best = opt::EMPTY;
			//values shared between both reads
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				const google::dense_hash_map<ID, unsigned>::const_iterator &j =
						hitCounts2.find(i->first);
				if (j != hitCounts2.end()
						&& bestHit <= (i->second + j->second)) {
					if (bestHit == (i->second + j->second) && best < i->first) {
						continue;
					}
					bestHit = i->second + j->second;
					best = i->first;
				}
			}
			//if ensure pair makes any sense
			if (best != 0) {
				hits.push_back(best);
			}
		} else {
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				unsigned hitCount = i->second;
				google::dense_hash_map<ID, unsigned>::const_iterator itr =
						hitCounts2.find(i->first);
				if (itr != hitCounts2.end()) {
					hitCount += itr->second;
				} else {
					continue;
				}
				if (bestHit < hitCount) {
					bestHit = hitCount;
				}
			}
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				unsigned hitCount = i->second;
				google::dense_hash_map<ID, unsigned>::const_iterator itr =
						hitCounts2.find(i->first);
				if (itr != hitCounts2.end()) {
					hitCount += itr->second;
				} else {
					continue;
				}
				if (bestHit == hitCount) {
					hits.push_back(i->first);
				}
			}
		}
		return bestHit;
	}

	/*
	 * Returns a vector of best hits to a specific read
	 * Only one read needs to match
	 * outputs additional score based on results
	 */
	//TODO: possible to optimize more!
	unsigned convertToHitsOnlyOne(
			const google::dense_hash_map<ID, unsigned> &hitCounts1,
			const google::dense_hash_map<ID, unsigned> &hitCounts2,
			vector<ID> &hits) {
		unsigned bestHit = 0;
		if (m_filter.getType() == MIBloomFilter<ID>::MIBFMVAL) {
			ID best = opt::EMPTY;
			//values shared between both reads
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts2.begin(); i != hitCounts2.end(); ++i) {
				if (bestHit <= i->second) {
					if (bestHit == i->second && best < i->first) {
						continue;
					}
					bestHit = i->second;
					best = i->first;
				}
			}
			//values found in only 1 read try to pick stronger one
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				const google::dense_hash_map<ID, unsigned>::const_iterator &j =
						hitCounts2.find(i->first);
				if (j != hitCounts2.end()) {
					if (bestHit <= (i->second + j->second)) {
						if (bestHit == (i->second + j->second)
								&& best < i->first) {
							continue;
						}
						bestHit = i->second + j->second;
						best = i->first;
					}
				} else {
					if (bestHit <= i->second) {
						if (bestHit == i->second && best < i->first) {
							continue;
						}
						bestHit = i->second;
						best = i->first;
					}
				}
			}
			hits.push_back(best);
		} else {

			//for tracking already found element in first set
			google::dense_hash_set<ID> tempSet;
			tempSet.set_empty_key(opt::EMPTY);

			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				unsigned hitCount = i->second;
				google::dense_hash_map<ID, unsigned>::const_iterator itr =
						hitCounts2.find(i->first);
				if (itr != hitCounts2.end()) {
					tempSet.insert(i->first);
					hitCount += itr->second;
				}
				if (bestHit < hitCount) {
					bestHit = hitCount;
				}
			}
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts2.begin(); i != hitCounts2.end(); ++i) {
				if (bestHit < i->second) {
					bestHit = i->second;
				}
			}

			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts1.begin(); i != hitCounts1.end(); ++i) {
				unsigned hitCount = i->second;
				google::dense_hash_map<ID, unsigned>::const_iterator itr =
						hitCounts2.find(i->first);
				if (itr != hitCounts2.end()) {
					hitCount += itr->second;
				}
				if (bestHit <= hitCount + opt::delta) {
					hits.push_back(i->first);
				}
			}
			for (google::dense_hash_map<ID, unsigned>::const_iterator i =
					hitCounts2.begin(); i != hitCounts2.end(); ++i) {
				if (tempSet.find(i->first) == tempSet.end()) {
					if (bestHit <= i->second + opt::delta) {
						hits.push_back(i->first);
					}
				}
			}
		}
		return bestHit;
	}
};

#endif /* BLOOMMAPCLASSIFIER_H_ */
