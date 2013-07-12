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
#include "config.h"

using namespace std;

#define PROGRAM "biobloommaker"

void printVersion()
{
	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
	"Written by Justin Chu.\n"
	"\n"
	"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog()
{
	static const char dialog[] =
			"Usage: biobloommaker -p [FILTERID] [OPTION]... [FILE]...\n"
					"Creates a bf and txt file from a list of fasta files. The input sequences are\n"
					"cut into a k-mers with a sliding window and their hash signatures are inserted\n"
					"into a bloom filter.\n"
					"\n"
					"  -p, --file_prefix=N    Filter prefix and filter ID. Required option.\n"
					"  -o, --output_dir=N     Output location of the filter and filter info files.\n"
					"  -h, --help             Display this dialog.\n"
					"\nAdvanced options:\n"
					"  -f, --fal_pos_rate=N   Maximum false positive rate to use in filter. [0.02]\n"
					"  -g, --hash_num=N       Set number of hash functions to use in filter instead\n"
					"                         of automatically using calculated optimal number of\n"
					"                         functions.\n"
					"  -k, --kmer_size        K-mer size to use to create filter. [25]\n"
					"\nOption presets:\n"
					"      --default          Run categorizer assuming default presets (ie. no\n"
					"                         advanced options toggled) [default]\n"
					"      --low_mem          Run categorizer assuming low memory presets.\n"
					"      --minimum_fpr      Run categorizer assuming minimized false rate presets.\n"
					"\n"
					"Report bugs to <cjustin@bcgsc.ca>.";
	cerr << dialog << endl;
	exit(0);
}

int main(int argc, char *argv[])
{

	bool die = false;

	//switch statement variable
	int c;

	//command line variables
	double fpr = 0.02;
	string filterPrefix = "";
	string outputDir = "";
	uint16_t kmerSize = 25;
	uint16_t hashNum = 0;

	//preset options
	int defaultSettings = 0;
	int lowMem = 0;
	int minimumFPR = 0;
	string presetType = "default";

	//long form arguments
	static struct option long_options[] = {
			{
					"fal_pos_rate", required_argument, NULL, 'f' }, {
					"file_prefix", required_argument, NULL, 'p' }, {
					"output_dir", required_argument, NULL, 'o' }, {
					"hash_num", required_argument, NULL, 'g' }, {
					"kmer_size", required_argument, NULL, 'k' }, {
					"default", no_argument, &defaultSettings, 0 }, {
					"low_mem", no_argument, &lowMem, 0 }, {
					"minimum_fpr", no_argument, &minimumFPR, 0 }, {
					"help", no_argument, NULL, 'h' }, {
					NULL, 0, NULL, 0 } };

	//check if only one preset was set
	if (defaultSettings || minimumFPR || lowMem) {
		if (!(defaultSettings ^ minimumFPR ^ lowMem)) {
			cerr << "Error: Cannot mix option presets" << endl;
			exit(1);
		}
	}

	//set presets
	if (lowMem) {
		fpr = 0.14;
		kmerSize = 24;
		presetType = "low_mem";
	} else if (minimumFPR) {
		kmerSize = 24;
		presetType = "minimum_fpr";
	}

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:p:o:k:n:g:hv", long_options,
			&option_index)) != -1)
	{
		switch (c) {
		case 'f': {
			stringstream convert(optarg);
			if (!(convert >> fpr)) {
				cerr << "Error - Invalid set of bloom filter parameters! f: "
						<< optarg << endl;
				return 0;
			}
			if (fpr > 1) {
				cerr << "Error -f cannot be greater than 1 " << optarg << endl;
				return 0;
			}
			presetType = "custom";
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
				cerr << "Error - Invalid set of bloom filter parameters! k: "
						<< optarg << endl;
				return 0;
			}
			presetType = "custom";
			break;
		}
		case 'g': {
			stringstream convert(optarg);
			if (!(convert >> hashNum)) {
				cerr << "Error - Invalid set of bloom filter parameters! g: "
						<< optarg << endl;
				return 0;
			}
			presetType = "custom";
			break;
		}
		case 'h': {
			printHelpDialog();
			break;
		}
		case 'v': {
			printVersion();
			break;
		}
		default: {
			die = true;
			break;
		}
		}
	}

	//Stores fasta input file names
	vector<string> inputFiles;

	while (optind < argc) {
		inputFiles.push_back(argv[optind]);
		optind++;
	}

	//Check needed options
	if (inputFiles.size() == 0) {
		cerr << "Need Input File" << endl;
		die = true;
	}
	if (filterPrefix.size() == 0) {
		cerr << "Need Filter Prefix ID" << endl;
		die = true;
	}
	if (filterPrefix.find('/') != string::npos) {
		cerr << "Prefix ID cannot have '/' characters" << endl;
		die = true;
	}
	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

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
			hashNum, presetType);

	//get calculated size of Filter
	size_t filterSize = info.getCalcuatedFilterSize();

	//Add seed to bloom filter and get them for info file
	vector<size_t> seeds = filterGen.addHashFuncs(hashNum);
	assert(seeds.size() == hashNum);
	filterGen.setFilterSize(filterSize);

	//output filter
	size_t redundNum = filterGen.generate(outputDir + filterPrefix + ".bf");
	info.setReduanacy(redundNum);

	//set hash function info
	vector<string> hashFunctions = filterGen.getHashFuncNames();

	for (uint16_t i = 0; i < hashNum; ++i) {
		info.addHashFunction(hashFunctions[i], seeds[i]);
	}

	//code for redundancy checking
	//calcuate redundancy rate
	double redunRate = double(redundNum) / double(entryNum)
			- info.getRedunancyFPR();
	if (redunRate > 0.2) {
		cerr << "Redunancy Rate is approximately: " << redunRate << endl;
		cerr
				<< "Consider checking your files for duplicate sequences and adjusting them accordingly."
				<< endl;
		cerr
				<< "High redundancy will cause overestimation of filter sizes used."
				<< endl;
	}

	//output info
	info.printInfoFile(outputDir + filterPrefix + ".txt");
	cerr << "Filter Creation Complete." << endl;

	return 0;
}
