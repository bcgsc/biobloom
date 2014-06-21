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
			size_t minHit, double percHit, size_t maxHitValue, uint8_t tileModifier);

	const string updateSummaryData(size_t seqLen,
			unordered_map<string, size_t> &hits);
	const string updateSummaryData(size_t seqLen1, size_t seqLen2,
			unordered_map<string, size_t> &hits1,
			unordered_map<string, size_t> &hits2);

	const string getResultsSummary(size_t readCount) const;
	const string getCountSummary(size_t readCount) const;
	virtual ~ResultsManager();
private:
	//Variables copied from biobloomcategorizer
	const unordered_map<string, vector<shared_ptr<BloomFilterInfo> > > infoFiles;
	const unordered_map<string, shared_ptr<MultiFilter> > filters;
	const vector<string> hashSigs;
	const size_t minHit;
	const double percentMinHit;
	const size_t maxHitValue;
	const uint8_t tileModifier;

	unordered_map<string, size_t> aboveThreshold;
	unordered_map<string, size_t> unique;
	unordered_map<string, vector<size_t> > rawHits;
	size_t multiMatch;
	size_t noMatch;
};

#endif /* RESULTSMANAGER_H_ */
