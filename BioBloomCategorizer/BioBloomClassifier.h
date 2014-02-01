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
#include "Common/BloomFilter.h"

using namespace std;
using namespace boost;

class BioBloomClassifier {
public:
	explicit BioBloomClassifier(const vector<string> &filterFilePaths,
			double scoreThreshold, const string &outputPrefix,
			const string &outputPostFix, uint16_t streakThreshold, uint16_t minHit,
			bool minHitOnly);
	void filter(const vector<string> &inputFiles);
	void filterPrint(const vector<string> &inputFiles,
			const string &outputType);
	void filterPair(const string &file1, const string &file2);
	void filterPairPrint(const string &file1, const string &file2,
			const string &outputType);
	void filterPairBAM(const string &file);
	void filterPairBAMPrint(const string &file, const string &outputType);

	virtual ~BioBloomClassifier();

private:
	void loadFilters(const vector<string> &filterFilePaths);
	const bool fexists(const string &filename) const;
	void evaluateReadStd(const FastqRecord &rec, const string &hashSig,
			unordered_map<string, bool> &hits);
	void evaluateRead(const FastqRecord &rec, const string &hashSig,
			unordered_map<string, bool> &hits);

	//group filters with same hash number
	unordered_map<string, vector<shared_ptr<BloomFilterInfo> > > infoFiles;
	unordered_map<string, shared_ptr<MultiFilter> > filters;
	unordered_map<string, shared_ptr<BloomFilter> > filtersSingle;
	vector<string> hashSigs;
	double scoreThreshold;
	uint8_t filterNum;
	const string &prefix;
	const string &postfix;
	const uint16_t streakThreshold;
	const uint16_t minHit;
	const bool minHitOnly;

	//Todo: is this really better than hard-coding them in the class?
	const string noMatch;
	const string multiMatch;
};

#endif /* BIOBLOOMCLASSIFIER_H_ */
