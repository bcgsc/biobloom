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

void printVersion() {
	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
	"Written by Justin Chu.\n"
	"\n"
	"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog() {
	static const char dialog[] =
			"Usage: biobloommaker -p [FILTERID] [OPTION]... [FILE]...\n"
					"Creates a bf and txt file from a list of fasta files. The input sequences are\n"
					"cut into a k-mers with a sliding window and their hash signatures are inserted\n"
					"into a bloom filter.\n"
					"\n"
					"  -p, --file_prefix=N    Filter prefix and filter ID. Required option.\n"
					"  -o, --output_dir=N     Output location of the filter and filter info files.\n"
					"  -h, --help             Display this dialog.\n"
					"  -v  --version          Display version information.\n"
					"\nAdvanced options:\n"
					"  -f, --fal_pos_rate=N   Maximum false positive rate to use in filter. [0.02]\n"
					"  -g, --hash_num=N       Set number of hash functions to use in filter instead\n"
					"                         of automatically using calculated optimal number of\n"
					"                         functions.\n"
					"  -k, --kmer_size=N      K-mer size to use to create filter. [25]\n"
					"  -s, --subtract=N       Path to filter that you want to uses to prevent the\n"
					"                         addition of k-mers contained into new filter. You may\n"
					"                         only use filters with k-mer sizes <= the one you wish\n"
					"                         to create.\n"
					"  -n, --num_ele=N        Set the number of expected elements. If set to 0 number\n"
					"                         is determined from sequences sizes within files. [0]"
					"\n"
					"Report bugs to <cjustin@bcgsc.ca>.";
	cerr << dialog << endl;
	exit(0);
}

int main(int argc, char *argv[]) {

	bool die = false;

	//switch statement variable
	int c;

	//command line variables
	double fpr = 0.02;
	string filterPrefix = "";
	string outputDir = "";
	uint16_t kmerSize = 25;
	uint16_t hashNum = 0;
	string subtractFilter = "";
	size_t entryNum = 0;

	//long form arguments
	static struct option long_options[] = { { "fal_pos_rate", required_argument,
	NULL, 'f' }, { "file_prefix", required_argument, NULL, 'p' }, {
			"output_dir", required_argument, NULL, 'o' }, { "version",
	no_argument, NULL, 'v' }, { "hash_num", required_argument, NULL, 'g' }, {
			"kmer_size",
			required_argument, NULL, 'k' }, { "subtract",
	required_argument, NULL, 's' }, { "num_ele",
	required_argument, NULL, 'n' }, { "help", no_argument, NULL, 'h' }, {
	NULL, 0, NULL, 0 } };

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:p:o:k:n:g:hvs:n:", long_options,
			&option_index)) != -1) {
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
			break;
		}
		case 'g': {
			stringstream convert(optarg);
			if (!(convert >> hashNum)) {
				cerr << "Error - Invalid set of bloom filter parameters! g: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'n': {
			stringstream convert(optarg);
			if (!(convert >> entryNum)) {
				cerr << "Error - Invalid set of bloom filter parameters! n: "
						<< optarg << endl;
				return 0;
			}
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
		case 's': {
			stringstream convert(optarg);
			if (!(convert >> subtractFilter)) {
				cerr << "Error - Invalid set of bloom filter parameters! s: "
						<< optarg << endl;
				return 0;
			}
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

	//set number of hash functions used
	if (hashNum == 0) {
		//get optimal number of hash functions
		hashNum = calcOptiHashNum(fpr);
	}

	//create filter
	BloomFilterGenerator filterGen(inputFiles, kmerSize, hashNum, entryNum);

	if (entryNum == 0) {
		filterGen = BloomFilterGenerator(inputFiles, kmerSize, hashNum);
		entryNum = filterGen.getExpectedEntries();
	}

	BloomFilterInfo info(filterPrefix, kmerSize, hashNum, fpr, entryNum,
			inputFiles);

	//get calculated size of Filter
	size_t filterSize = info.getCalcuatedFilterSize();
	filterGen.setFilterSize(filterSize);

	size_t redundNum = 0;
	//output filter
	if (subtractFilter == "") {
		redundNum = filterGen.generate(outputDir + filterPrefix + ".bf");
	} else {
		redundNum = filterGen.generate(outputDir + filterPrefix + ".bf",
				subtractFilter);
	}
	info.setTotalNum(filterGen.getTotalEntries());
	info.setRedundancy(redundNum);

	//code for redundancy checking
	//calcuate redundancy rate
	double redunRate = double(redundNum) / double(entryNum)
			- info.getRedundancyFPR();
	if (redunRate > 0.25) {
		cerr << "Redundancy Rate is approximately: " << redunRate << endl;
		cerr
				<< "Consider checking your files for duplicate sequences and adjusting them accordingly.\n"
						"High redundancy will cause filter sizes used overestimated, potentially resulting in a larger than needed filter.\n"
						"Alternatively you can set the size of filter wanted with (-n) and ignore this message."
				<< endl;
	}

	//output info
	info.printInfoFile(outputDir + filterPrefix + ".txt");
	cerr << "Filter Creation Complete." << endl;

	return 0;
}
