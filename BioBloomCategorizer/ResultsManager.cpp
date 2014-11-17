/*
 * ResultsManager.cpp
 *
 *  Created on: Jun 25, 2013
 *      Author: cjustin
 */

#include "ResultsManager.h"
#include <sstream>
#include <iostream>
#if _OPENMP
# include <omp.h>
#endif

ResultsManager::ResultsManager(const vector<string> &hashSigsRef,
		const unordered_map<string, shared_ptr<MultiFilter> > &filtersRef,
		const unordered_map<string, vector<shared_ptr<BloomFilterInfo> > > &infoFilesRef,
		double scoreThreshold) :
		m_hashSigs(hashSigsRef), m_filters(filtersRef), m_infoFiles(
				infoFilesRef), m_multiMatch(0), m_noMatch(0)
{
	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = m_hashSigs.begin();
			j != m_hashSigs.end(); ++j)
	{
		const shared_ptr<MultiFilter> &temp = m_filters.at(*j);
		const vector<string> &idsInFilter = temp->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			m_aboveThreshold[*i] = 0;
			m_unique[*i] = 0;
		}
	}
}

/*
 * Records data for read status file based on thresholds
 * Returns qualifying read IDs that meet threshold
 */
const string ResultsManager::updateSummaryData(
		const unordered_map<string, bool> &hits)
{
	string filterID;
	bool noMatchFlag = true;
	bool multiMatchFlag = false;

	for (vector<string>::const_iterator j = m_hashSigs.begin();
			j != m_hashSigs.end(); ++j)
	{
		//update summary
		const shared_ptr<MultiFilter> &temp = m_filters.at(*j);
		const vector<string> &idsInFilter = temp->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			if (hits.at(*i)) {
#pragma omp atomic
				++m_aboveThreshold[*i];
				if (noMatchFlag) {
					noMatchFlag = false;
					filterID = *i;
				} else {
					multiMatchFlag = true;
				}
			}
		}
	}
	if (noMatchFlag) {
		filterID = "m_noMatch";
#pragma omp atomic
		++m_noMatch;
	} else {
		if (multiMatchFlag) {
			filterID = "m_multiMatch";
#pragma omp atomic
			++m_multiMatch;
		} else {
#pragma omp atomic
			++m_unique[filterID];
		}
	}
	return filterID;
}

/*
 * Records data for read status file based on thresholds
 * Returns qualifying read IDs that meet threshold
 * both reads must qualify
 */
const string ResultsManager::updateSummaryData(
		const unordered_map<string, bool> &hits1,
		const unordered_map<string, bool> &hits2)
{
	string filterID;
	bool noMatchFlag = true;
	bool multiMatchFlag = false;

	for (vector<string>::const_iterator j = m_hashSigs.begin();
			j != m_hashSigs.end(); ++j)
	{
		//update summary
		const shared_ptr<MultiFilter> &temp = m_filters.at(*j);
		const vector<string> &idsInFilter = temp->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			if (hits1.at(*i) && hits2.at(*i)) {
#pragma omp atomic
				++m_aboveThreshold[*i];
				if (noMatchFlag) {
					noMatchFlag = false;
					filterID = *i;
				} else {
					multiMatchFlag = true;
				}
			}
		}
	}
	if (noMatchFlag) {
		filterID = "m_noMatch";
#pragma omp atomic
		++m_noMatch;
	} else {
		if (multiMatchFlag) {
			filterID = "m_multiMatch";
#pragma omp atomic
			++m_multiMatch;
		} else {
#pragma omp atomic
			++m_unique[filterID];
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

	for (vector<string>::const_iterator j = m_hashSigs.begin();
			j != m_hashSigs.end(); ++j)
	{
		const shared_ptr<MultiFilter> &temp = m_filters.at(*j);
		const vector<string> &idsInFilter = temp->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			const vector<shared_ptr<BloomFilterInfo> > &tempVect =
					m_infoFiles.at(*j);
			summaryOutput << *i << "_" << tempVect.front()->getKmerSize();
			summaryOutput << "\t" << m_aboveThreshold.at(*i);
			summaryOutput << "\t" << readCount - m_aboveThreshold.at(*i);
			summaryOutput << "\t"
					<< (m_aboveThreshold.at(*i) - m_unique.at(*i));
			summaryOutput << "\t"
					<< double(m_aboveThreshold.at(*i)) / double(readCount);
			summaryOutput << "\t"
					<< double(readCount - m_aboveThreshold.at(*i))
							/ double(readCount);
			summaryOutput << "\t"
					<< double(m_aboveThreshold.at(*i) - m_unique.at(*i))
							/ double(readCount);
			summaryOutput << "\n";
		}
	}

	summaryOutput << "m_multiMatch";
	summaryOutput << "\t" << m_multiMatch;
	summaryOutput << "\t" << readCount - m_multiMatch;
	summaryOutput << "\t" << 0;
	summaryOutput << "\t" << double(m_multiMatch) / double(readCount);
	summaryOutput << "\t"
			<< double(readCount - m_multiMatch) / double(readCount);
	summaryOutput << "\t" << 0.0;
	summaryOutput << "\n";

	summaryOutput << "m_noMatch";
	summaryOutput << "\t" << m_noMatch;
	summaryOutput << "\t" << readCount - m_noMatch;
	summaryOutput << "\t" << 0;
	summaryOutput << "\t" << double(m_noMatch) / double(readCount);
	summaryOutput << "\t" << double(readCount - m_noMatch) / double(readCount);
	summaryOutput << "\t" << 0.0;
	summaryOutput << "\n";

	cerr << summaryOutput.str() << endl;
	return summaryOutput.str();
}

ResultsManager::~ResultsManager()
{
}

