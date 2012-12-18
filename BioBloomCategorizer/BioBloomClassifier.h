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
			double percentMinHit);
	void filter(const vector<string> &inputFiles, const string &outputPrefix);
	void filterPrintReads(const vector<string> &inputFiles,
			const string &outputPrefix);
	void filterPair(const string &file1, const string &file2,
			const string &outputPrefix);

	virtual ~BioBloomClassifier();
private:
	void loadFilters(const vector<string> &filterFilePaths);
	//todo: should be const qualifier, but getting a compilation bug
	void printSummary(const string &outputDir,
			boost::unordered_map<string, size_t> &aboveThreshold,
			boost::unordered_map<string, size_t> &belowThreshold,
			size_t totalReads);
	bool fexists(const string &filename) const;

	//group filters with same hash signature
	boost::unordered_map<string, vector<boost::shared_ptr<BloomFilterInfo> > > infoFiles;
	boost::unordered_map<string, boost::shared_ptr<MultiFilter> > filters;
	vector<string> hashSigs;
	int16_t minHit;
	double percentMinHit;
	int16_t filterNum;
};

#endif /* BIOBLOOMCLASSIFIER_H_ */
