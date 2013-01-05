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
#include "Common/HashManager.h"

using namespace std;

class BioBloomClassifier {
public:
	BioBloomClassifier(const vector<string> &filterFilePaths, int16_t minHit,
			double percentMinHit, int16_t maxHitValue);
	void filter(const vector<string> &inputFiles, const string &outputPrefix);
	void filterPrintReads(const vector<string> &inputFiles,
			const string &outputPrefix);
	void filterPair(const string &file1, const string &file2,
			const string &outputPrefix);
	void filterPairPrint(const string &file1, const string &file2,
			const string &outputPrefix);

	virtual ~BioBloomClassifier();
private:
	void loadFilters(const vector<string> &filterFilePaths);
	//Todo: should be const qualified, need to refactor function
	void printSummary(const string &outputPrefix,
			boost::unordered_map<string, size_t> &aboveThreshold,
			boost::unordered_map<string, size_t> &belowThreshold,
			size_t totalReads);
	void printCountSummary(const string &outputPrefix,
			boost::unordered_map<string, vector<size_t> > &rawHits,
			size_t total, size_t nonATCG);
	bool fexists(const string &filename) const;
	bool evaluateRead(const FastqRecord &rec, const string &hashSig,
			boost::unordered_map<string, size_t> &hits);

	//group filters with same hash signature
	boost::unordered_map<string, vector<boost::shared_ptr<BloomFilterInfo> > > infoFiles;
	boost::unordered_map<string, boost::shared_ptr<MultiFilter> > filters;
	vector<string> hashSigs;
	int16_t minHit;
	double percentMinHit;
	int16_t filterNum;
	int16_t maxHitValue;
};

#endif /* BIOBLOOMCLASSIFIER_H_ */
