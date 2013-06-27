/*
 * ResultsManager.h
 *
 *  Created on: Jun 25, 2013
 *      Author: cjustin
 */

#ifndef RESULTSMANAGER_H_
#define RESULTSMANAGER_H_

#include <string>
#include <vector>
#include "boost/unordered/unordered_map.hpp"
#include "boost/shared_ptr.hpp"
#include "MultiFilter.h"

using namespace std;
using namespace boost;

class ResultsManager {
public:
	explicit ResultsManager(const vector<string> &hashSigs,
			const unordered_map<string, shared_ptr<MultiFilter> > &filtersRef,
			const unordered_map<string, vector<shared_ptr<BloomFilterInfo> > > &infoFiles,
			size_t minHit, double percHit, size_t maxHitValue);

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
	unordered_map<string, vector<shared_ptr<BloomFilterInfo> > > &infoFiles;
	unordered_map<string, shared_ptr<MultiFilter> > &filters;
	vector<string> &hashSigs;
	size_t minHit;
	double percentMinHit;
	size_t maxHitValue;

	unordered_map<string, size_t> aboveThreshold;
	unordered_map<string, size_t> belowThreshold;
	unordered_map<string, vector<const size_t> > rawHits;
};

#endif /* RESULTSMANAGER_H_ */
