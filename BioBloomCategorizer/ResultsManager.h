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
#include "Common/Options.h"

using namespace std;
using namespace boost;

class ResultsManager {
public:
	explicit ResultsManager(const vector<string> &m_filterOrder,
			bool inclusive);

	//TODO: replace usage with vector<ID> version
	const string updateSummaryData(const unordered_map<string, bool> &hits);

	const string updateSummaryData(const vector<ID> &hits);

	//TODO: replace usage with vector<ID> version
	const string updateSummaryData(const unordered_map<string, bool> &hits1,
			const unordered_map<string, bool> &hits2);

	const string getResultsSummary(size_t readCount) const;
	virtual ~ResultsManager();
private:
	//Variables copied from biobloomcategorizer
	const vector<string> &m_filterOrder;

	size_t m_multiMatch;
	size_t m_noMatch;
	bool m_inclusive;

	//TODO: change to IDs vector type (for speed)
	vector<size_t> m_aboveThreshold;
	vector<size_t> m_unique;
};

#endif /* RESULTSMANAGER_H_ */
