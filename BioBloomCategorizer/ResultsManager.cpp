/*
 * ResultsManager.cpp
 *
 *  Created on: Jun 25, 2013
 *      Author: cjustin
 */

#include "ResultsManager.h"
#include <sstream>

ResultsManager::ResultsManager(const vector<string> &hashSigsRef,
		const unordered_map<string, shared_ptr<MultiFilter> > &filtersRef,
		const unordered_map<string, vector<shared_ptr<BloomFilterInfo> > > &infoFilesRef,
		size_t minHit, double percHit, size_t maxHitValue, uint8_t tileModifier) :
		hashSigs(hashSigsRef), filters(filtersRef), infoFiles(infoFilesRef), minHit(
				minHit), percentMinHit(percHit), maxHitValue(maxHitValue), tileModifier(
				tileModifier), multiMatch(0), noMatch(0)
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

			//initialize rawHits
			vector<size_t> temp(maxHitValue);
			uint16_t counter = 0;
			rawHits[*i] = temp;
			while (counter < maxHitValue) {
				rawHits[*i][counter] = 0;
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
			//pick threshold, by percent or by absolute value
			const vector<shared_ptr<BloomFilterInfo> > &tempVect = infoFiles.at(
					*j);
			uint16_t kmerSize = tempVect.front()->getKmerSize();
			size_t threshold = size_t(
					percentMinHit * (seqLen / (kmerSize + tileModifier)));
			if (minHit > threshold) {
				threshold = minHit;
			}

			if (hits[*i] >= threshold) {
				++aboveThreshold[*i];
				if (noMatchFlag) {
					noMatchFlag = false;
					filterID = *i;
				} else {
					multiMatchFlag = true;
				}
			}

			//modify total reads
			if (rawHits[*i].size() > hits[*i]) {
				rawHits[*i][hits[*i]]++;
			}
		}
	}
	if (noMatchFlag) {
		filterID = "noMatch";
		++noMatch;
	} else {
		if (multiMatchFlag) {
			filterID = "multiMatch";
			++multiMatch;
		} else {
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
const string ResultsManager::updateSummaryData(size_t seqLen1, size_t seqLen2,
		unordered_map<string, size_t> &hits1,
		unordered_map<string, size_t> &hits2)
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
			//pick threshold, by percent or by absolute value
			const vector<shared_ptr<BloomFilterInfo> > &tempVect = infoFiles.at(
					*j);
			uint16_t kmerSize = (*(tempVect.front())).getKmerSize();
			size_t threshold1 = size_t(
					percentMinHit * (seqLen1 / (kmerSize + tileModifier)));
			size_t threshold2 = size_t(
					percentMinHit * (seqLen2 / (kmerSize + tileModifier)));
			if (minHit > threshold1) {
				threshold1 = minHit;
			}
			if (minHit > threshold2) {
				threshold2 = minHit;
			}

			if (hits1[*i] >= threshold1 && hits2[*i] >= threshold2) {
				++aboveThreshold[*i];
				if (noMatchFlag) {
					noMatchFlag = false;
					filterID = *i;
				} else {
					multiMatchFlag = true;
				}
			}
			//modify total reads
			if (rawHits[*i].size() > hits1[*i] + hits2[*i]) {
				rawHits[*i][hits1[*i] + hits2[*i]]++;
			}
		}
	}
	if (noMatchFlag) {
		filterID = "noMatch";
		++noMatch;
	} else {
		if (multiMatchFlag) {
			filterID = "multiMatch";
			++multiMatch;
		} else {
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
		const shared_ptr<MultiFilter> &temp = filters.at(*j);
		const vector<string> &idsInFilter = temp->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			const vector<shared_ptr<BloomFilterInfo> > &tempVect = infoFiles.at(
					*j);
			summaryOutput << "\t" << *i << "_"
					<< tempVect.front()->getKmerSize();
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
			const shared_ptr<MultiFilter> &temp = filters.at(*j);
			const vector<string> &idsInFilter = temp->getFilterIds();
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				const vector<size_t> &temp = rawHits.at(*i);
				if (temp.size() < currentHitVal) {
					summaryOutput << "\t0";
				} else {
					runningTotal += temp[currentHitVal];
					summaryOutput << "\t" << temp[currentHitVal];
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

