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
			const unordered_map<string, boost::shared_ptr<MultiFilter> > &filtersRef,
			const unordered_map<string, vector<boost::shared_ptr<BloomFilterInfo> > > &infoFilesRef,
			bool inclusive);

	const string updateSummaryData(const unordered_map<string, bool> &hits);
	const string updateSummaryData(const unordered_map<string, bool> &hits1,
			const unordered_map<string, bool> &hits2);

	const string getResultsSummary(size_t readCount) const;
	virtual ~ResultsManager();
private:
	//Variables copied from biobloomcategorizer
	const vector<string> m_hashSigs;
	const unordered_map<string, boost::shared_ptr<MultiFilter> > m_filters;
	const unordered_map<string, vector<boost::shared_ptr<BloomFilterInfo> > > m_infoFiles;

	bool m_inclusive;

	unordered_map<string, size_t> m_aboveThreshold;
	unordered_map<string, size_t> m_unique;
	size_t m_multiMatch;
	size_t m_noMatch;
};

#endif /* RESULTSMANAGER_H_ */
