/*
 * BioBloom.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: cjustin
 */

#include <sstream>
#include <string>
#include <vector>
#include "BloomFilterGenerator.h"
#include "Common/BloomFilterInfo.h"
#include <boost/unordered/unordered_map.hpp>
#include <getopt.h>
using namespace std;

/*
 * Parses the input string to store each files name with a list of headers
 * in a map full of vectors
 */
//todo: resolve the fact that this code is also found in BloomFilterInfo
boost::unordered_map<string, vector<string> > convertInputString(
		const string &inputString) {
	boost::unordered_map<string, vector<string> > inputs;
	vector<string> currentList;
	string currentFileName = "";
	string temp;
	stringstream converter(inputString);
	while (converter >> temp) {
		if (temp.at(temp.length() - 1) == ',') {
			temp.resize(temp.length() - 1);
			if (currentFileName.empty()) {
				currentFileName = temp;
			} else {
				currentList.push_back(temp);
			}
			inputs[currentFileName] = currentList;
			currentFileName = "";
		} else {
			if (currentFileName.empty()) {
				currentFileName = temp;
			} else {
				currentList.push_back(temp);
			}
		}
	}
	inputs[currentFileName] = currentList;
	return inputs;
}

int main(int argc, char *argv[]) {

	//switch statement variable
	int c;

	//command line variables
	string rawInputFiles = "";
	double fpr = 1;
	string filterPrefix = "";
	string outputDir = "";
	int16_t kmerSize = 0;
	int16_t hashNum = 0;

	//long form arguments
	//each option format { "optionName", necessary option or not, I have no idea, 'symbol'}
	static struct option long_options[] = { { "input_file_statement", 1, NULL,
			'i' }, { "false_positive_rate", 1, NULL, 'f' }, { "file_prefix", 1,
			NULL, 'p' }, { "output_dir", 0, NULL, 'o' }, { "separate", 0, NULL,
			's' }, { "hash_num", 0, NULL, 'g' }, { "kmer_size", 1, NULL, 'k' },
			{ NULL, 0, NULL, 0 } };

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "i:f:p:o:k:n:", long_options,
			&option_index)) != -1) {
		switch (c) {
		case 'i': {
			rawInputFiles = optarg;
			break;
		}
		case 'f': {
			stringstream convert(optarg);
			if (!(convert >> fpr)) {
				cout << "Error - Invalid set of bloom filter parameters! f: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'p': {
			filterPrefix = optarg;
			break;
		}
		case 'o': {
			outputDir = optarg;
			if (outputDir.at(outputDir.length() - 1) != '/') {
				outputDir = outputDir + '/';
			}
			break;
		}
		case 'k': {
			stringstream convert(optarg);
			if (!(convert >> kmerSize)) {
				cout << "Error - Invalid set of bloom filter parameters! k: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'g': {
			stringstream convert(optarg);
			if (!(convert >> hashNum)) {
				cout << "Error - Invalid set of bloom filter parameters! g: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		}
	}

	boost::unordered_map<string, vector<string> > inputFiles =
			convertInputString(rawInputFiles);
	//create filter
	BloomFilterGenerator filterGen(inputFiles, kmerSize);

	//get expected number of entries
	size_t entryNum = filterGen.getExpectedEntries();

	//set number of hash functions used
	if (hashNum == 0) {
		//get optimal number of hash functions
		hashNum = filterGen.calcOptiHashNum(fpr);
	}

	BloomFilterInfo info(filterPrefix, kmerSize, fpr, entryNum, inputFiles,
			hashNum);

	//get calcuated size of Filter
	size_t filterSize = info.getCalcuatedFilterSize();

	//Add seed to bloom filter and get them for info file
	vector<size_t> seeds = filterGen.addHashFuncs(hashNum);
	assert(seeds.size() == hashNum);
	filterGen.setFilterSize(filterSize);

	//output filter
	filterGen.generate(outputDir + filterPrefix + ".bf");
	//set hashfunction info
	vector<string> hashFunctions = filterGen.getHashFuncNames();

	for (int16_t i = 0; i < hashNum; ++i) {
		info.addHashFunction(hashFunctions[i], seeds[i]);
	}

	//output info
	info.printInfoFile(outputDir + filterPrefix + ".txt");
	cout << "Filter Creation Complete" << endl;

	return 0;
}
