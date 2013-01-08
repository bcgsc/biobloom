/*
 * BloomFilterInfo.cpp
 *
 *  Created on: Aug 20, 2012
 *      Author: cjustin
 */
#include "BloomFilterInfo.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

BloomFilterInfo::BloomFilterInfo(string const &filterID, uint16_t kmerSize,
		float desiredFPR, size_t numEntries, const vector<string> &seqSrcs,
		uint16_t hashNum) :
		filterID(filterID), kmerSize(kmerSize), desiredFPR(desiredFPR), seqSrcs(
				seqSrcs), hashNum(hashNum)
{
	runInfo.numEntries = numEntries;
	runInfo.size = calcOptimalSize(numEntries, desiredFPR, hashNum);
}

/*
 * loads bloom filter information from a file
 */
//Todo: convert to having variables stored in property tree for more modularity
BloomFilterInfo::BloomFilterInfo(string const &fileName)
{
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(fileName, pt);
	filterID = pt.get<string>("user_input_options.filter_id");
	kmerSize = pt.get<uint16_t>("user_input_options.kmer_size");
	desiredFPR = pt.get<float>("user_input_options.desired_false_positve_rate");
	string tempSeqSrcs = pt.get<string>("user_input_options.sequence_sources");
	seqSrcs = convertSeqSrcString(tempSeqSrcs);

	//runtime params
	runInfo.size = pt.get<size_t>("runtime_options.size");
	runInfo.numEntries = pt.get<size_t>("runtime_options.num_entries");
	float tempFPR = pt.get<float>(
			"runtime_options.approximate_false_positive_rate");
	string tempHashFunc = pt.get<string>("runtime_options.hash_functions");
	string tempSeeds = pt.get<string>("runtime_options.seeds");
	runInfo.hashFunctions = convertHashFuncString(tempHashFunc);
	runInfo.seeds = convertSeedString(tempSeeds);
}

//todo: change so that class is not mutable after object construction
void BloomFilterInfo::addHashFunction(const string &fnName, size_t seed)
{
	runInfo.hashFunctions.push_back(fnName);
	runInfo.seeds.push_back(seed);
}

/*
 * Prints out INI format file
 */
void BloomFilterInfo::printInfoFile(const string &fileName) const
{

	//calcuate actual FPR
	float maxFPR = calcApproxFPR(runInfo.size, runInfo.numEntries,
			runInfo.seeds.size());

	//cannot output unless runtime values set
	assert(runInfo.size !=0);
	assert(runInfo.numEntries !=0);
	assert(maxFPR !=0);
	assert(runInfo.hashFunctions.size() !=0);
	assert(runInfo.seeds.size() !=0);
	assert(hashNum == runInfo.hashFunctions.size());

	ofstream output(fileName.c_str(), ios::out);
	//user specified
	output << "[user_input_options]\nfilter_id=" << filterID << "\nkmer_size="
			<< kmerSize << "\ndesired_false_positve_rate=" << desiredFPR
			<< "\nsequence_sources=";
	//print out sources as a list
	size_t tempCounter = 0;
	for (vector<string>::const_iterator it = seqSrcs.begin();
			it != seqSrcs.end(); ++it)
	{
		output << *it;
		output << " ";
	}

	//runtime determined options
	output << "\n\n[runtime_options]\nsize=" << runInfo.size << "\nnum_entries="
			<< runInfo.numEntries << "\napproximate_false_positive_rate="
			<< maxFPR << "\n";
	//print out hash functions as a list

	output << getSeedHashSigniture();

	output.close();
}

//getters

const uint16_t BloomFilterInfo::getKmerSize() const
{
	return kmerSize;
}

/*
 * returns a unique string representing the hash function and seeds used
 */
const string BloomFilterInfo::getSeedHashSigniture() const
{
	stringstream ss;
	ss << "hash_functions=";
	//print out hash functions as a list
	uint16_t tempCounter = 0;
	for (vector<string>::const_iterator it = runInfo.hashFunctions.begin();
			it != runInfo.hashFunctions.end(); ++it)
	{
		ss << (*it);
		//print closing quotes
		++tempCounter;
		if (tempCounter == runInfo.hashFunctions.size()) {
			ss << "";
		} else {
			ss << ", ";
		}
	}
	ss << "\nseeds=";
	tempCounter = 0;
	for (vector<size_t>::const_iterator it = runInfo.seeds.begin();
			it != runInfo.seeds.end(); ++it)
	{
		ss << *it;
		++tempCounter;
		if (tempCounter != runInfo.seeds.size()) {
			ss << ", ";
		}
	}
	return ss.str();
}

const size_t BloomFilterInfo::getCalcuatedFilterSize() const
{
	return runInfo.size;
}

const vector<string> &BloomFilterInfo::getHashFunctionNames() const
{
	return runInfo.hashFunctions;

}

const vector<size_t> &BloomFilterInfo::getSeedValues() const
{
	return runInfo.seeds;
}

const string &BloomFilterInfo::getFilterID() const
{
	return filterID;
}

const vector<string> BloomFilterInfo::convertSeqSrcString(
		string const &seqSrcStr) const
{
	vector<string> inputs;
	string currentFileName = "";
	string temp;
	stringstream converter(seqSrcStr);
	while (converter >> temp) {
		inputs.push_back(temp);
	}
	return inputs;
}

/*
 * converts string obtained from hash values in info file and converts to list
 */
const vector<string> BloomFilterInfo::convertHashFuncString(
		const string &hashFnStr) const
{
	vector<string> fnNames;
	string temp;
	stringstream converter(hashFnStr);
	while (converter >> temp) {
		if (temp.at(temp.length() - 1) == ',') {
			temp.resize(temp.length() - 1);
		}
		fnNames.push_back(temp);
	}
	return fnNames;

}

/*
 * converts string obtained from seeds in info file and converts to list
 */
const vector<size_t> BloomFilterInfo::convertSeedString(
		const string &seedStr) const
{
	vector<size_t> seedVal;
	string temp;
	stringstream converter(seedStr);
	while (converter >> temp) {
		if (temp.at(temp.length() - 1) == ',') {
			temp.resize(temp.length() - 1);
		}
		size_t num;
		stringstream converter2(temp);
		converter2 >> num;
		seedVal.push_back(num);
	}
	return seedVal;
}
// functions for calculations regarding bloomfilter

// todo: Tweak calculations as they are approximations and may not be 100% optimal
// see http://en.wikipedia.org/wiki/Bloom_filter

//Private functions
/*
 * Value uses approximated FPR formula
 */
const float BloomFilterInfo::calcApproxFPR(size_t size, size_t numEntr,
		uint16_t hashFunctNum) const
{
	return pow(
			1.0 - exp(-double(hashFunctNum) * double(numEntr) / double(size)),
			hashFunctNum);
}

/*
 * Only returns multiples of 64 for filter building purposes
 * Is an estimated size using approximations of FPR formula
 * assuming optimal # of hash functions used
 * see http://en.wikipedia.org/wiki/Bloom_filter
 */
//NOTE: Not currently used.
const size_t BloomFilterInfo::calcOptimalSize(size_t entries, float fpr) const
{
	size_t non64ApproxVal = size_t(entries * -log(fpr) / pow(log(2), 2));
	return non64ApproxVal + (64 - non64ApproxVal % 64);
}

/*
 * Only returns multiples of 64 for filter building purposes
 * Is an estimated size using approximations of FPR formula
 * given the number of hash functions
 * see http://en.wikipedia.org/wiki/Bloom_filter
 */
const size_t BloomFilterInfo::calcOptimalSize(size_t entries, float fpr,
		uint16_t hashNum) const
{
	size_t non64ApproxVal = size_t(
			-double(entries) * double(hashNum)
					/ log(1.0 - pow(fpr, float(1 / (float(hashNum))))));

	return non64ApproxVal + (64 - non64ApproxVal % 64);
}

BloomFilterInfo::~BloomFilterInfo()
{
}

