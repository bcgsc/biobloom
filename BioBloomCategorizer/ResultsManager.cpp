/*
 * ResultsManager.cpp
 *
 *  Created on: Jun 25, 2013
 *      Author: cjustin
 */

#include "ResultsManager.h"
#include <sstream>

ResultsManager::ResultsManager(const vector<string> &hashSigs,
		const unordered_map<string, shared_ptr<MultiFilter> > &filtersRef,
		const unordered_map<string, vector<shared_ptr<BloomFilterInfo> > > &infoFiles,
		size_t minHit, double percHit, size_t maxHitValue) :
		hashSigs(hashSigs), filters(filtersRef), infoFiles(infoFiles), minHit(
				minHit), percentMinHit(percHit), maxHitValue(maxHitValue)
{
	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> &idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			aboveThreshold[*i] = 0;
			belowThreshold[*i] = 0;

			//initialize rawHits
			vector<size_t> temp(maxHitValue);
			uint16_t counter = 0;
			rawHits[*i] = temp;
			while (counter < maxHitValue) {
				(*rawHits[*i])[counter] = 0;
				counter++;
			}
		}
	}
}

/*
 * Records data for read status file based on thresholds
 * Returns qualifying read IDs that meet threshold
 */
const string ResultsManager::updateSummaryData(size_t seqLen,
		unordered_map<string, size_t> &hits)
{
	string filterID = "noMatch";

	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		//update summary
		const vector<string> &idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			//pick threshold, by percent or by absolute value
			uint16_t kmerSize = (*(infoFiles[*j].front())).getKmerSize();
			size_t threshold = size_t(percentMinHit * (seqLen / kmerSize));
			if (minHit > threshold) {
				threshold = minHit;
			}

			if (hits[*i] >= threshold) {
				++aboveThreshold[*i];
				if (filterID == "noMatch") {
					filterID = *i;
				} else {
					filterID = "multiMatch";
				}
			} else if (hits[*i] != 0) {
				++belowThreshold[*i];
			}

			//modify total reads
			if (rawHits[*i].size() > hits[*i]) {
				(*rawHits[*i])[hits[*i]]++;
			}
		}
	}
	return filterID;
}

/*
 * Records data for read status file based on thresholds
 * Returns qualifying read IDs that meet threshold
 * both reads must qualify
 */
const string ResultsManager::updateSummaryData(size_t seqLen1,
		size_t seqLen2, unordered_map<string, size_t> &hits1,
		unordered_map<string, size_t> &hits2)
{
	string filterID = "noMatch";

	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		//update summary
		const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			//pick threshold, by percent or by absolute value
			uint16_t kmerSize = (*(infoFiles[*j].front())).getKmerSize();
			size_t threshold1 = size_t(percentMinHit * (seqLen1 / kmerSize));
			size_t threshold2 = size_t(percentMinHit * (seqLen2 / kmerSize));
			if (minHit > threshold1) {
				threshold1 = minHit;
			}
			if (minHit > threshold2) {
				threshold2 = minHit;
			}

			if (hits1[*i] >= threshold1 && hits2[*i] >= threshold2) {
				++aboveThreshold[*i];
				if (filterID == "noMatch") {
					filterID = *i;
				} else {
					filterID = "multiMatch";
				}
			} else if (hits1[*i] != 0 && hits2[*i] != 0) {
				++belowThreshold[*i];
			}
			//modify total reads
			if (rawHits[*i].size() > hits1[*i] + hits2[*i]) {
				rawHits[*i][hits1[*i] + hits2[*i]]++;
			}
		}
	}
	return filterID;
}

const string ResultsManager::getResultsSummary(size_t readCount) const
{

	stringstream summaryOutput;
	summaryOutput << "type";
	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> &idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << *i << "_"
					<< (*(infoFiles[*j].front())).getKmerSize();
		}
	}
	summaryOutput << "\n";

	//print summary information and close filehandles
	summaryOutput << "\nHits";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(aboveThreshold[*i]) / double(readCount);
		}
	}
	summaryOutput << "\nMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(belowThreshold[*i]) / double(readCount);
		}
	}
	summaryOutput << "\nConfidentMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(
							readCount - belowThreshold[*i] - aboveThreshold[*i])
							/ double(readCount);
		}
	}

	summaryOutput << "\n\nHits";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << aboveThreshold[*i];
		}
	}
	summaryOutput << "\nMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << belowThreshold[*i];
		}
	}
	summaryOutput << "\nConfidentMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< readCount - belowThreshold[*i] - aboveThreshold[*i];
		}
	}
	return summaryOutput.str();
}
const string ResultsManager::getCountSummary(size_t readCount) const
{
	stringstream summaryOutput;
	summaryOutput << "type";
	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> &idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << *i << "_"
					<< (*(infoFiles[*j].front())).getKmerSize();
		}
	}
	summaryOutput << "\n";

	uint16_t currentHitVal = 0;
	size_t runningTotal = 0;

	while (currentHitVal < maxHitValue) {
		summaryOutput << currentHitVal;
		for (vector<string>::const_iterator j = hashSigs.begin();
				j != hashSigs.end(); ++j)
		{
			const vector<string> &idsInFilter = filters[*j]->getFilterIds();
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				const vector<size_t> &temp = rawHits[*i];
				if ( temp.size() < currentHitVal) {
					summaryOutput << "\t0";
				} else {
					runningTotal += (*rawHits[*i])[currentHitVal];
					summaryOutput << "\t" << (*rawHits[*i])[currentHitVal];
				}
			}
		}
		summaryOutput << "\n";
		currentHitVal++;
	}
	summaryOutput << ">" << maxHitValue << "\t" << readCount - runningTotal
			<< "\n";
	return summaryOutput.str();
}

ResultsManager::~ResultsManager()
{
}

