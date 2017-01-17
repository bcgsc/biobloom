/*
 * BioBloomRecruit.cpp
 *
 *  Created on: Jan 14, 2017
 *      Author: cjustin
 */

#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include "Common/SeqEval.h"
#include <getopt.h>
#include "config.h"
#include "Common/Options.h"
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "bbr"

void printVersion() {
	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	"Written by Justin Chu.\n"
	"\n"
	"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog() {
	static const char dialog[] =
		"Usage: bbr -p [PREFIX] -s [REPEATFILE] -b [score] [BAITFASTA] [FASTQ1] [FASTQ2] \n"
		"Recruits sequences from fastq files using bait sequences. Accepts fa and gzip.\n"
		"\n"
		"  -p, --prefix=STR       Filter prefix and filter ID. Required option.\n"
		"  -h, --help             Display this dialog.\n"
		"  -v, --version          Display version information.\n"
		"  -t, --threads=INT      The number of threads to use. Experimental. [1]\n"
		"                         Currently only active with the (-r) option.\n"
		"  -f, --fpr=NUM          Maximum false positive rate to use in filter.\n"
		"                         For Bloom Maps this value reference to the occupancy rate\n"
		"                         of the filter rather than the FPR [0.0078125]\n"
		"  -g, --hash_num=INT     Set number of hash functions to use in filter instead\n"
		"                         of automatically using calculated optimal number of\n"
		"                         functions.\n"
		"  -k, --kmer_size=INT    K-mer size to use to create filter. [64]\n"
		"  -s, --subtract=STR     Path to filter that you want to uses to prevent the\n"
		"                         addition of k-mers contained into new filter. You may\n"
		"                         only use filters with k-mer sizes equal the one you\n"
		"                         wish to create. Use this to minimize repeat propagation\n"
		"                         when generating progressive filters.\n"
		"  -n, --num_ele=INT      Set the number of expected elements. If set to 0 number\n"
		"                         is determined from sequences sizes within files. [0]\n"
		"  -r, --progressive=NUM  Progressive filter creation. The score threshold is\n"
		"                         specified by N, which may be either a floating point\n"
		"                         score between 0 and 1 or a positive integer.  If N is a\n"
		"                         positive integer, it is interpreted as the minimum\n"
		"                         number of contiguous matching bases required for a\n"
		"                         match. [0.5]\n"
		"  -b, --baitScore=NUM    Score threshold when considering only bait. [r]\n"
		"  -e, --iterations=INT   Pass through files N times if threshold is not met."
		"  -i, --inclusive        If one paired read matches, both reads will be included\n"
		"                         in the filter. Only active with the (-r) option.\n"
		"\n"
		"Report bugs to <cjustin@bcgsc.ca>.";
	cerr << dialog << endl;
	exit(0);
}

int main(int argc, char *argv[]) {

	bool die = false;

	//switch statement variable
	int c;

//	//command line variables
//	double fpr = 0.0075;
//	string filterPrefix = "";
//	string outputDir = "";
//	unsigned kmerSize = 25;
//	unsigned hashNum = 0;
//	string subtractFilter = "";
//	size_t entryNum = 0;
//	bool printReads = false;
//	double progressive = -1;
//	bool inclusive = false;
//	SeqEval::EvalMode evalMode = SeqEval::EVAL_STANDARD;

	//long form arguments
	static struct option long_options[] = {
		{
			"fpr", required_argument, NULL, 'f' }, {
			"prefix", required_argument, NULL, 'p' }, {
			"threads", required_argument, NULL, 't' }, {
			"version", no_argument, NULL, 'v' }, {
			"hash_num", required_argument, NULL, 'g' }, {
			"kmer_size", required_argument, NULL, 'k' }, {
			"subtract",	required_argument, NULL, 's' }, {
			"num_ele", required_argument, NULL, 'n' }, {
			"help", no_argument, NULL, 'h' }, {
			"print_reads", no_argument, NULL, 'P' }, {
			"progressive", required_argument, NULL, 'r' }, {
			"baitScore", required_argument, NULL, 'b' }, {
			"iterations", required_argument, NULL, 'e' }, {
			NULL, 0, NULL, 0 } };
	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:p:k:n:g:hvs:n:t:r:b:e:", long_options,
			&option_index)) != -1) {
		switch (c) {
		case 'f': {
			stringstream convert(optarg);
			if (!(convert >> opt::fpr)) {
				cerr << "Error - Invalid set of bloom filter parameters! f: "
						<< optarg << endl;
				return 0;
			}
			if (opt::fpr > 1) {
				cerr << "Error -f cannot be greater than 1 " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'p': {
			stringstream convert(optarg);
			convert >> opt::prefix;
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
			if (!(convert >> opt::numEle)) {
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
			if (!(convert >> opt::subtract)) {
				cerr << "Error - Invalid set of bloom filter parameters! s: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'r': {
			stringstream convert(optarg);
			unsigned matchLen;
			// if arg is a positive integer > 1, interpret as minimum match
			// length in bases
			if ((convert >> matchLen) && matchLen > 1) {
				opt::pScore = matchLen;
				opt::mode = SeqEval::EVAL_MIN_MATCH_LEN;
			} else {
				// not a positive integer > 1, so interpret as floating
				// point score between 0 and 1
				stringstream convert2(optarg);
				if (!(convert2 >> opt::pScore)) {
					cerr << "Error - Invalid set of bloom filter parameters! r: "
						<< optarg << endl;
					return 0;
				}
				if (opt::pScore < 0 || opt::pScore > 1) {
					cerr << "Error - r must be a positive integer or a floating "
						<< "point between 0 and 1. Input given:"
						<< optarg << endl;
					exit(EXIT_FAILURE);
				}
			}
			break;
		}
		case 'b': {
			stringstream convert(optarg);
			if (!(convert >> opt::baitThreshold)) {
				cerr << "Error - Invalid set of bloom filter parameters! b: "
						<< optarg << endl;
				return 0;
			}
			if (opt::baitThreshold < 0) {
				cerr << "Error - b must be a positive integer or a floating "
						<< "point between 0 and 1. Input given:" << optarg
						<< endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'e': {
			stringstream convert(optarg);
			if (!(convert >> opt::progItrns)) {
				cerr << "Error - Invalid set of bloom filter parameters! b: "
						<< optarg << endl;
				return 0;
			}
			if (opt::progItrns > 0) {
				cerr << "Error - e must be > 1" << optarg
						<< endl;
				exit(EXIT_FAILURE);
			}
			break;
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
	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	//set number of hash functions used
	if (opt::hashNum == 0) {
		//get optimal number of hash functions
		opt::hashNum = unsigned(-log(opt::fpr) / log(2));
	}

	//START HERE

	return 0;
}
