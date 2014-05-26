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
		hashSigs(hashSigsRef), filters(filtersRef), infoFiles(infoFilesRef), scoreThreshold(
				scoreThreshold), multiMatch(0), noMatch(0)
{
	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const shared_ptr<MultiFilter> &temp = filters.at(*j);
		const vector<string> &idsInFilter = temp->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			aboveThreshold[*i] = 0;
			unique[*i] = 0;
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

	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		//update summary
		const shared_ptr<MultiFilter> &temp = filters.at(*j);
		const vector<string> &idsInFilter = temp->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			if (hits.at(*i)) {
#pragma omp atomic
				++aboveThreshold[*i];
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
		filterID = "noMatch";
#pragma omp atomic
		++noMatch;
	} else {
		if (multiMatchFlag) {
			filterID = "multiMatch";
#pragma omp atomic
			++multiMatch;
		} else {
#pragma omp atomic
			++unique[filterID];
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

	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		//update summary
		const shared_ptr<MultiFilter> &temp = filters.at(*j);
		const vector<string> &idsInFilter = temp->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			if (hits1.at(*i) && hits2.at(*i)) {
#pragma omp atomic
				++aboveThreshold[*i];
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
		filterID = "noMatch";
#pragma omp atomic
		++noMatch;
	} else {
		if (multiMatchFlag) {
			filterID = "multiMatch";
#pragma omp atomic
			++multiMatch;
		} else {
#pragma omp atomic
			++unique[filterID];
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

	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const shared_ptr<MultiFilter> &temp = filters.at(*j);
		const vector<string> &idsInFilter = temp->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			const vector<shared_ptr<BloomFilterInfo> > &tempVect = infoFiles.at(
					*j);
			summaryOutput << *i << "_" << tempVect.front()->getKmerSize();
			summaryOutput << "\t" << aboveThreshold.at(*i);
			summaryOutput << "\t" << readCount - aboveThreshold.at(*i);
			summaryOutput << "\t" << (aboveThreshold.at(*i) - unique.at(*i));
			summaryOutput << "\t"
					<< double(aboveThreshold.at(*i)) / double(readCount);
			summaryOutput << "\t"
					<< double(readCount - aboveThreshold.at(*i))
							/ double(readCount);
			summaryOutput << "\t"
					<< double(aboveThreshold.at(*i) - unique.at(*i))
							/ double(readCount);
			summaryOutput << "\n";
		}
	}

	summaryOutput << "multiMatch";
	summaryOutput << "\t" << multiMatch;
	summaryOutput << "\t" << readCount - multiMatch;
	summaryOutput << "\t" << 0;
	summaryOutput << "\t" << double(multiMatch) / double(readCount);
	summaryOutput << "\t" << double(readCount - multiMatch) / double(readCount);
	summaryOutput << "\t" << 0.0;
	summaryOutput << "\n";

	summaryOutput << "noMatch";
	summaryOutput << "\t" << noMatch;
	summaryOutput << "\t" << readCount - noMatch;
	summaryOutput << "\t" << 0;
	summaryOutput << "\t" << double(noMatch) / double(readCount);
	summaryOutput << "\t" << double(readCount - noMatch) / double(readCount);
	summaryOutput << "\t" << 0.0;
	summaryOutput << "\n";

	cout << summaryOutput.str() << endl;
	return summaryOutput.str();
}

ResultsManager::~ResultsManager()
{
}

