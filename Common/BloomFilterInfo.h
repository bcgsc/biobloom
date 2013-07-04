/*
 * BloomFilterInfo.h
 * Intended to be used to collect/store and compute information about bloom filters
 * Can output information in IEE format text file
 *
 *  Created on: Aug 20, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTERINFO_H_
#define BLOOMFILTERINFO_H_
#include <string>
#include <vector>
#include <boost/unordered/unordered_map.hpp>

using namespace std;

class BloomFilterInfo {
public:
	explicit BloomFilterInfo(string const &filterID, uint16_t kmerSize,
			float desiredFPR, size_t numEntries, const vector<string> &seqSrcs,
			uint16_t hashNum);
	explicit BloomFilterInfo(string const &fileName);
	void addHashFunction(const string &fnName, size_t seed);
	void setReduanacy(size_t redunSeq);

	void printInfoFile(const string &fileName) const;
	virtual ~BloomFilterInfo();

	//getters
	const uint16_t getKmerSize() const;
	const size_t getCalcuatedFilterSize() const;
	const vector<string> &getHashFunctionNames() const;
	const vector<size_t> &getSeedValues() const;
	const string getSeedHashSigniture() const;
	const string &getFilterID() const;
	const double &getRedunancyFPR() const;

private:
	//user specified input
	string filterID;
	uint16_t kmerSize;
	float desiredFPR;
	vector<string> seqSrcs;
	uint16_t hashNum;

	//determined at run time
	struct runtime {
		size_t size;
		size_t numEntries;
		vector<string> hashFunctions;
		vector<size_t> seeds;
		size_t redundantSequences;
		double redundantFPR;
	};

	runtime runInfo;

	const vector<string> convertSeqSrcString(const string &seqSrcStr) const;
	const vector<string> convertHashFuncString(const string &hashFnStr) const;
	const vector<size_t> convertSeedString(const string &seedStr) const;
	const double calcApproxFPR(size_t size, size_t numEntr,
			uint16_t hashFunctNum) const;
	const double calcRedunancyFPR(size_t size, size_t numEntr,
			uint16_t hashFunctNum, size_t redundantSeqs) const;
	const size_t calcOptimalSize(size_t entries, float fpr) const;
	const size_t calcOptimalSize(size_t entries, float fpr,
			uint16_t hashNum) const;
};

#endif /* BLOOMFILTERINFO_H_ */
