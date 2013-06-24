/*
 * BioBloomClassifier.h
 *
 *  Created on: Oct 17, 2012
 *      Author: cjustin
 */

#ifndef BIOBLOOMCLASSIFIER_H_
#define BIOBLOOMCLASSIFIER_H_
#include <vector>
#include <string>
#include "boost/unordered/unordered_map.hpp"
#include "boost/shared_ptr.hpp"
#include "Common/BloomFilterInfo.h"
#include "MultiFilter.h"
#include "DataLayer/FastaReader.h"
#include "Common/ReadsProcessor.h"
#include "Common/Uncompress.h"

using namespace std;
using namespace boost;

class BioBloomClassifier {
public:
	explicit BioBloomClassifier(const vector<string> &filterFilePaths, size_t minHit,
			double percentMinHit, size_t maxHitValue,
			const string &outputPrefix, const string &outputPostFix, uint8_t tileModifier);
	void filter(const vector<string> &inputFiles);
	void filterPrint(const vector<string> &inputFiles);
	void filterPair(const string &file1, const string &file2);
	void filterPairPrint(const string &file1, const string &file2);
	void filterPairBAM(const string &file);
	void filterPairBAMPrint(const string &file);

	virtual ~BioBloomClassifier();
private:
	void loadFilters(const vector<string> &filterFilePaths);
	//Todo: should be const qualified, need to refactor function
	void printSummary(const string &outputPrefix,
			unordered_map<string, size_t> &aboveThreshold,
			unordered_map<string, size_t> &belowThreshold, size_t totalReads);
	void printCountSummary(const string &outputPrefix,
			unordered_map<string, vector<size_t> > &rawHits, size_t total);
	bool fexists(const string &filename) const;
	void evaluateRead(const FastqRecord &rec, const string &hashSig,
			unordered_map<string, size_t> &hits, uint8_t tileModifier);
	const string initSummaryVars(vector<string> &hashSigs,
			unordered_map<string, size_t> &aboveThreshold,
			unordered_map<string, size_t> &belowThreshold,
			unordered_map<string, vector<size_t> > &rawHits);
	void initHits(unordered_map<string, size_t> &hits);
	const string getReadStatStr(string const &readID, size_t readLength,
			unordered_map<string, size_t> &hits);
	const string getReadStatStrPair(string const &readID, size_t readLength1,
			size_t readLength2, unordered_map<string, size_t> &hits1,
			unordered_map<string, size_t> &hits2);
	const string updateSummaryData(size_t seqLen,
			unordered_map<string, size_t> &hits,
			unordered_map<string, size_t> &aboveThreshold,
			unordered_map<string, size_t> &belowThreshold,
			unordered_map<string, vector<size_t> > &rawHits);
	const string updateSummaryDataPair(size_t seqLen1, size_t seqLen2,
			unordered_map<string, size_t> &hits1,
			unordered_map<string, size_t> &hits2,
			unordered_map<string, size_t> &aboveThreshold,
			unordered_map<string, size_t> &belowThreshold,
			unordered_map<string, vector<size_t> > &rawHits);

	//group filters with same hash signature
	unordered_map<string, vector<shared_ptr<BloomFilterInfo> > > infoFiles;
	unordered_map<string, shared_ptr<MultiFilter> > filters;
	vector<string> hashSigs;
	size_t minHit;
	double percentMinHit;
	uint8_t filterNum;
	size_t maxHitValue;
	const string &postfix;
	const string &prefix;
	const uint8_t tileModifier;

	//Todo: is this really better than hard-coding them in the class?
	const string noMatch;
	const string multiMatch;
};

#endif /* BIOBLOOMCLASSIFIER_H_ */
