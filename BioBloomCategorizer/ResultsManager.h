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
#include "Common/BloomFilterInfo.h"
#include "Common/Options.h"

using namespace std;

class ResultsManager {
public:
	explicit ResultsManager(const vector<string> &m_filterOrder,
			bool inclusive);

	unsigned updateSummaryData(const vector<unsigned> &hits);
	unsigned updateSummaryData(const vector<unsigned> &hits1,
			const vector<unsigned> &hits2);

	const string getResultsSummary(size_t readCount) const;
	virtual ~ResultsManager();
private:
	//Variables owned by biobloomcategorizer
	const vector<string> &m_filterOrder;

	const unsigned m_noMatchIndex;
	const unsigned m_multiMatchIndex;

	vector<size_t> m_aboveThreshold;
	vector<size_t> m_unique;
	size_t m_multiMatch;
	size_t m_noMatch;
	bool m_inclusive;
};

#endif /* RESULTSMANAGER_H_ */
