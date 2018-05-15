/*
 * BioBloom.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: cjustin
 */

#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <getopt.h>
#include "config.h"
#include "Options.h"
#include "Common/Options.h"
#include <fstream>
#include <omp.h>
#include "MIBFGen.hpp"

using namespace std;

#define PROGRAM "biobloommimaker"

void printVersion() {
	const char VERSION_MESSAGE[] =
	PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	"Written by Justin Chu.\n"
	"\n"
	"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

/*
 * Parses input string into separate strings, returning a vector.
 */
vector<string> convertInputString(const string &inputString)
{
	vector<string> currentString;
	string temp;
	stringstream converter(inputString);
	while (converter >> temp) {
		currentString.push_back(temp);
	}
	assert(currentString.size() > 0);
	return currentString;
}

void printHelpDialog() {
	static const char dialog[] =
		"Usage: biobloommimaker -p [FILTERID] [OPTION]... [FILE]...\n"
		"Usage: biobloommimaker -p [FILTERID] -r 0.2 [FILE]... [FASTQ1] [FASTQ2] \n"
		"Creates a multi-index Bloom Filter from a list of fasta files.\n"
		"\n"
		"  -p, --file_prefix=N    Filter prefix and filter ID. Required option.\n"
		"  -h, --help             Display this dialog.\n"
		"      --version          Display version information.\n"
		"  -v  --verbose          Display verbose output.\n"
		"  -t, --threads=N        The number of threads to use.\n"
		"  -b, --occupancy=N      Occupancy of Bloom filter.[0.5]\n"
		"  -n, --num_ele=N        Set the number of expected distinct k-mer frames.\n"
		"                         If set to 0 number is determined from sequences\n"
		"                         sizes within files. [0]\n"
		"  -S, --seed_str=N       Generate a miBF using multiple spaced seeds. Expects\n"
		"                         list of seed 1s & 0s separated by spaces.\n"
		"  -F, --by_file          Assign IDs by file rather than by fasta header\n"
		"k-mer mode options (disabled when using spaced seeds):\n"
		"  -g, --hash_num=N       Set number of hash functions when using k-mers.\n"
		"  -k, --kmer_size=N      K-mer size to use to create filter. [25]\n"
		"\n"
		"Report bugs to <cjustin@bcgsc.ca>.";
	cerr << dialog << endl;
	exit(0);
}

enum {
	OPT_VERSION
};


int main(int argc, char *argv[]) {

	bool die = false;

	//switch statement variable
	int c;

	//long form arguments
	static struct option long_options[] = {
		{
			"file_prefix", required_argument, NULL, 'p' }, {
			"help", no_argument, NULL, 'h' }, {
			"threads", required_argument, NULL, 't' }, {
			"occupancy", required_argument, NULL, 'b' }, {
			"seed_str", required_argument, NULL, 'S' }, {
			"hash_num", required_argument, NULL, 'g' }, {
			"kmer_size", required_argument, NULL, 'k' }, {
			"num_ele", required_argument, NULL, 'n' }, {
			"by_file", no_argument, NULL, 'F' }, {
			"verbose", no_argument, NULL, 'v' }, {
			"version", no_argument, NULL, OPT_VERSION }, {
			NULL, 0, NULL, 0 } };

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "p:ht:b:S:g:k:n:Fv",
			long_options, &option_index)) != -1) {
		switch (c) {
		case 'b': {
			stringstream convert(optarg);
			if (!(convert >> opt::occupancy)) {
				cerr << "Error - Invalid set of bloom filter parameters! b: "
						<< optarg << endl;
				return 0;
			}
			if (opt::fpr > 1) {
				cerr << "Error -b cannot be greater than 1 " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'p': {
			opt::prefix = optarg;
			break;
		}
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				cerr << "Error - Invalid parameter! t: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'S': {
			opt::sseeds = convertInputString(optarg);
			opt::kmerSize = opt::sseeds[0].size();
			break;
		}
		case 'F': {
			opt::idByFile = true;
			break;
		}
		case 'k': {
			stringstream convert(optarg);
			if (!(convert >> opt::kmerSize)) {
				cerr << "Error - Invalid set of bloom filter parameters! k: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'g': {
			stringstream convert(optarg);
			if (!(convert >> opt::hashNum)) {
				cerr << "Error - Invalid set of bloom filter parameters! g: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'n': {
			stringstream convert(optarg);
			if (!(convert >> opt::entryNum)) {
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
			opt::verbose++;
			break;
		}
//		case 'a': {
//			stringstream convert(optarg);
//			if (!(convert >> opt::allowMisses)) {
//				cerr << "Error - Invalid parameter! a: " << optarg << endl;
//				exit(EXIT_FAILURE);
//			}
//			break;
//		}
		case OPT_VERSION:{
			printVersion();
			exit(EXIT_SUCCESS);
		}
		default: {
			die = true;
			break;
		}
		}
	}

#if defined(_OPENMP)
	if (opt::threads > 0)
	omp_set_num_threads(opt::threads);
#endif

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
	if (opt::prefix.size() == 0) {
		cerr << "Need Filter Prefix ID" << endl;
		die = true;
	}
	if (opt::prefix.find('/') != string::npos) {
		cerr << "Prefix ID cannot have '/' characters" << endl;
		die = true;
	}
	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	//set number of hash functions used
	if (!opt::sseeds.empty()) {
		opt::hashNum = opt::sseeds.size();
	} else {
		cerr << "Please pick number of hash values (-g)\n";
		exit(1);
	}

	MIBFGen filterGen(inputFiles, opt::kmerSize, opt::entryNum);
	filterGen.generate(opt::prefix, opt::occupancy);
	if (opt::verbose) {
		cerr << "Multi Index Bloom Filter Creation Complete." << endl;
	}
	return 0;
}
