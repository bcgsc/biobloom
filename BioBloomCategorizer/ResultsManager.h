/*
 * ResultsManager.h
 *
 *  Created on: Jun 25, 2013
 *      Author: cjustin
 */

#ifndef RESULTSMANAGER_H_
#define RESULTSMANAGER_H_

#include <vector>
#include <string>
#include "boost/unordered/unordered_map.hpp"
#include "boost/shared_ptr.hpp"
#include "Common/BloomFilterInfo.h"
#include "MultiFilter.h"

using namespace std;
using namespace boost;

class ResultsManager {
public:
	explicit ResultsManager(const vector<string> &hashSigsRef,
			const unordered_map<string, shared_ptr<MultiFilter> > &filtersRef,
			const unordered_map<string, vector<shared_ptr<BloomFilterInfo> > > &infoFilesRef,
			float minHit, size_t maxHitValue,
			uint8_t tileModifier);

	const string updateSummaryData(size_t seqLen,
			unordered_map<string, float> &hits);
	const string updateSummaryData(size_t seqLen1, size_t seqLen2,
			unordered_map<string, float> &hits1,
			unordered_map<string, float> &hits2);

	const string getResultsSummary(size_t readCount) const;
	const string getCountSummary(size_t readCount) const;
	virtual ~ResultsManager();
private:
	//Variables copied from biobloomcategorizer
	const unordered_map<string, vector<shared_ptr<BloomFilterInfo> > > infoFiles;
	const unordered_map<string, shared_ptr<MultiFilter> > filters;
	const vector<string> hashSigs;
	const float minHit;
	const size_t maxHitValue;
	const uint8_t tileModifier;

	unordered_map<string, size_t> aboveThreshold;
	unordered_map<string, size_t> belowThreshold;
	unordered_map<string, vector<size_t> > rawHits;
};

#endif /* RESULTSMANAGER_H_ */
