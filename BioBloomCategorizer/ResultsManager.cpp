/*
 * ResultsManager.cpp
 *
 *  Created on: Jun 25, 2013
 *      Author: cjustin
 */

#include "ResultsManager.h"
#include <sstream>
#include <iostream>
#include <BioBloomClassifier.h>
#if _OPENMP
# include <omp.h>
#endif

ResultsManager::ResultsManager(const vector<string> &filterOrderRef,
		bool inclusive) :
		m_filterOrder(filterOrderRef), m_multiMatch(0), m_noMatch(0), m_unknown(
				0), m_inclusive(inclusive), m_aboveThreshold(
				vector<size_t>(filterOrderRef.size(), 0)), m_unique(
				vector<size_t>(filterOrderRef.size(), 0)) {
}

/*
 * Records data for read summary based on thresholds
 * Returns filter ID that this read equals
 */
const string ResultsManager::updateSummaryData(
		const unordered_map<string, bool> &hits)
{
	string filterID;
	unsigned filterIdx = 0;
	bool noMatchFlag = true;
	bool multiMatchFlag = false;

	for (unsigned i = 0; i < m_filterOrder.size(); ++i) {
		if (hits.at(m_filterOrder.at(i))) {
#pragma omp atomic
			++m_aboveThreshold[i];
			if (noMatchFlag) {
				noMatchFlag = false;
				filterID = m_filterOrder.at(i);
				filterIdx = i;
			} else {
				multiMatchFlag = true;
			}
		}
	}
	if (noMatchFlag) {
		filterID = NO_MATCH;
#pragma omp atomic
		++m_noMatch;
	} else {
		if (multiMatchFlag) {
			filterID = MULTI_MATCH;
#pragma omp atomic
			++m_multiMatch;
		} else {
#pragma omp atomic
			++m_unique[filterIdx];
		}
	}
	return filterID;
}

/*
 * Records data for read summary based on thresholds
 * Returns filter ID index (index -1) that this read equals
 */
ID ResultsManager::updateSummaryData(const vector<ID> &hits,
		bool aboveThreshold) {
	ID filterID = 0;

	for (vector<ID>::const_iterator i = hits.begin(); i != hits.end(); ++i) {
#pragma omp atomic
		++m_aboveThreshold[*i];
		filterID = *i;
	}
	if (hits.size() == 0) {
		if (aboveThreshold) {
			++m_unknown;
			return opt::COLLI;
		} else {
#pragma omp atomic
			++m_noMatch;
			return opt::EMPTY;
		}
	} else {
		if (hits.size() > 1) {
#pragma omp atomic
			++m_multiMatch;
			//TODO return collisionID?
			return opt::COLLI;
		} else {
#pragma omp atomic
			++m_unique[filterID];
		}
	}
	return filterID;
}

/*
 * Records data for read summary based on thresholds
 * Returns filter ID that this read pair equals
 */
const string ResultsManager::updateSummaryData(
		const unordered_map<string, bool> &hits1,
		const unordered_map<string, bool> &hits2)
{
	string filterID;
	unsigned filterIdx = 0;
	bool noMatchFlag = true;
	bool multiMatchFlag = false;

	for (unsigned i = 0; i < m_filterOrder.size(); ++i) {
		if (m_inclusive) {
			if (hits1.at(m_filterOrder.at(i)) || hits2.at(m_filterOrder.at(i))) {
#pragma omp atomic
				++m_aboveThreshold[i];
				if (noMatchFlag) {
					noMatchFlag = false;
					filterID = m_filterOrder.at(i);
					filterIdx = i;
				} else {
					multiMatchFlag = true;
				}
			}
		} else {
			if (hits1.at(m_filterOrder.at(i)) && hits2.at(m_filterOrder.at(i))) {
#pragma omp atomic
				++m_aboveThreshold[i];
				if (noMatchFlag) {
					noMatchFlag = false;
					filterID = m_filterOrder.at(i);
					filterIdx = i;
				} else {
					multiMatchFlag = true;
				}
			}
		}
	}
	if (noMatchFlag) {
		filterID = NO_MATCH;
#pragma omp atomic
		++m_noMatch;
	} else {
		if (multiMatchFlag) {
			filterID = MULTI_MATCH;
#pragma omp atomic
			++m_multiMatch;
		} else {
#pragma omp atomic
			++m_unique[filterIdx];
		}
	}
	return filterID;
}

const string ResultsManager::getResultsSummary(size_t readCount) const
{

	stringstream summaryOutput;

	//print header
	summaryOutput
			<< "filter_id\thits\tmisses\tshared\trate_hit\trate_miss\trate_shared\n";

	for (unsigned i = 1; i < m_filterOrder.size(); ++i) {
		summaryOutput << m_filterOrder[i];
		summaryOutput << "\t" << m_aboveThreshold.at(i);
		summaryOutput << "\t" << readCount - m_aboveThreshold.at(i);
		summaryOutput << "\t"
				<< (m_aboveThreshold.at(i)
						- m_unique.at(i));
		summaryOutput << "\t"
				<< double(m_aboveThreshold.at(i))
						/ double(readCount);
		summaryOutput << "\t"
				<< double(readCount - m_aboveThreshold.at(i))
						/ double(readCount);
		summaryOutput << "\t"
				<< double(m_aboveThreshold.at(i) - m_unique.at(i))
						/ double(readCount);
		summaryOutput << "\n";
	}

	summaryOutput << MULTI_MATCH;
	summaryOutput << "\t" << m_multiMatch;
	summaryOutput << "\t" << readCount - m_multiMatch;
	summaryOutput << "\t" << 0;
	summaryOutput << "\t" << double(m_multiMatch) / double(readCount);
	summaryOutput << "\t"
			<< double(readCount - m_multiMatch) / double(readCount);
	summaryOutput << "\t" << 0.0;
	summaryOutput << "\n";

	summaryOutput << UNKNOWN;
	summaryOutput << "\t" << m_unknown;
	summaryOutput << "\t" << readCount - m_unknown;
	summaryOutput << "\t" << 0;
	summaryOutput << "\t" << double(m_unknown) / double(readCount);
	summaryOutput << "\t"
			<< double(readCount - m_unknown) / double(readCount);
	summaryOutput << "\t" << 0.0;
	summaryOutput << "\n";

	summaryOutput << NO_MATCH;
	summaryOutput << "\t" << m_noMatch;
	summaryOutput << "\t" << readCount - m_noMatch;
	summaryOutput << "\t" << 0;
	summaryOutput << "\t" << double(m_noMatch) / double(readCount);
	summaryOutput << "\t" << double(readCount - m_noMatch) / double(readCount);
	summaryOutput << "\t" << 0.0;
	summaryOutput << "\n";

	return summaryOutput.str();
}

ResultsManager::~ResultsManager()
{
}

