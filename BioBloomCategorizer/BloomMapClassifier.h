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

//	/*
//	 * Returns the number of hits
//	 */
//	inline unsigned evaluateRead(const string &seq,
//			google::dense_hash_map<ID, unsigned> &hitCounts,
//			google::dense_hash_map<ID, unsigned> &solidHitCounts,
//			unsigned &solidCounts) {
//		RollingHashIterator itr(seq, m_filter.getKmerSize(),
//				m_filter.getSeedValues());
//		unsigned nonZeroCount = 0;
//
//		while (itr != itr.end()) {
//			unsigned missCount = 0;
//			vector<ID> ids = m_filter.at(*itr, m_colliIDs, missCount);
//			if (missCount == 0) {
//				for (unsigned i = 0; i < ids.size(); ++i) {
//					ID id = ids[i];
//					if (solidHitCounts.find(id) != solidHitCounts.end()) {
//						++solidHitCounts[id];
//					} else {
//						solidHitCounts[id] = 1;
//					}
//					if (hitCounts.find(id) != hitCounts.end()) {
//						++hitCounts[id];
//					} else {
//						hitCounts[id] = 1;
//					}
//				}
//				++solidCounts;
//				++nonZeroCount;
//			} else if (missCount <= opt::allowMisses) {
//				for (unsigned i = 0; i < ids.size(); ++i) {
//					ID id = ids[i];
//					if (hitCounts.find(id) != hitCounts.end()) {
//						++hitCounts[id];
//					} else {
//						hitCounts[id] = 1;
//					}
//				}
//				++nonZeroCount;
//			}
//			++itr;
//		}
//		return nonZeroCount;
//	}

	/*
	 * Returns the number of hits
	 */
	inline unsigned evaluateRead(const string &seq,
			google::dense_hash_map<ID, unsigned> &hitCounts) {
		RollingHashIterator itr(seq, m_filter.getKmerSize(),
				m_filter.getSeedValues());
		unsigned nonZeroCount = 0;
		while (itr != itr.end()) {
			unsigned misses = 0;
			vector<ID> ids = m_filter.at(*itr, m_colliIDs, opt::allowMisses, misses);
			nonZeroCount += ids.size() > 0;
			nonZeroCount += misses == 0 * opt::allowMisses;
			for (unsigned i = 0; i < ids.size(); ++i) {
				ID id = ids[i];
				if (hitCounts.find(id) == hitCounts.end()) {
					hitCounts[id] = 0;
				}
				hitCounts[id] += ids.size() > 0;
				hitCounts[id] += misses == 0 * opt::allowMisses;
			}
			++itr;
		}
		return nonZeroCount;
	}

	/*
	 * Returns a vector of hits to a specific ID, favoring smaller IDs
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
	void convertToHitsBoth(
			const google::dense_hash_map<ID, unsigned> &hitCounts1,
			const google::dense_hash_map<ID, unsigned> &hitCounts2,
			vector<ID> &hits) {
		unsigned bestHit = 0;
		for (google::dense_hash_map<ID, unsigned>::const_iterator i =
				hitCounts1.begin(); i != hitCounts1.end(); ++i) {
			unsigned hitCount = i->second;
			google::dense_hash_map<ID, unsigned>::const_iterator itr =
					hitCounts2.find(i->first);
			if (itr != hitCounts2.end()) {
				hitCount += itr->second;
			}
			else{
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
			}
			else {
				continue;
			}
			if (bestHit == hitCount) {
				hits.push_back(i->first);
			}
		}
	}

	/*
	 * Returns a vector of best hits to a specific read
	 * Only one read needs to match
	 * outputs additional score based on results
	 */
	//TODO: possible to optimize more!
	void convertToHitsOnlyOne(
			const google::dense_hash_map<ID, unsigned> &hitCounts1,
			const google::dense_hash_map<ID, unsigned> &hitCounts2,
			vector<ID> &hits) {
		unsigned bestHit = 0;

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
};

#endif /* BLOOMMAPCLASSIFIER_H_ */
