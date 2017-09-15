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
#include "BloomFilterGenerator.h"
#include "Common/BloomFilterInfo.h"
#include "Common/SeqEval.h"
#include <getopt.h>
#include "config.h"
#include "Common/Options.h"
#include <fstream>
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "biobloommaker"

void printVersion() {
	const char VERSION_MESSAGE[] =
	PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	"Written by Justin Chu.\n"
	"\n"
	"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog() {
	static const char dialog[] =
		"Usage: biobloommaker -p [FILTERID] [OPTION]... [FILE]...\n"
		"Usage: biobloommaker -p [FILTERID] -r 0.2 [FILE]... [FASTQ1] [FASTQ2] \n"
		"Creates a bf and txt file from a list of fasta files. The input sequences are\n"
		"cut into a k-mers with a sliding window and their hash signatures are inserted\n"
		"into a bloom filter.\n"
		"\n"
		"  -p, --file_prefix=N    Filter prefix and filter ID. Required option.\n"
		"  -o, --output_dir=N     Output location of the filter and filter info files.\n"
		"  -h, --help             Display this dialog.\n"
		"  -v  --version          Display version information.\n"
		"  -t, --threads=N        The number of threads to use.\n"
		"\nAdvanced options:\n"
		"  -f, --fal_pos_rate=N   Maximum false positive rate to use in filter. [0.0075]\n"
		"  -g, --hash_num=N       Set number of hash functions to use in filter instead\n"
		"                         of automatically using calculated optimal number of\n"
		"                         functions.\n"
		"  -k, --kmer_size=N      K-mer size to use to create filter. [25]\n"
		"  -s, --subtract=N       Path to filter that you want to uses to prevent the\n"
		"                         addition of k-mers contained into new filter. You may\n"
		"                         only use filters with k-mer sizes equal the one you\n"
		"                         wish to create. Use this to minimize repeat propagation\n"
		"                         when generating progressive filters.\n"
		"  -d, --no_rep_kmer      Remove all repeat k-mers from the resulting filter in\n"
		"                         progressive mode.\n"
		"  -n, --num_ele=N        Set the number of expected elements. If set to 0 number\n"
		"                         is determined from sequences sizes within files. [0]\n"
		"  -I, --interval         the interval to report file processing status [10000000]\n"
		"\nOptions for progressive filters:\n"
		"  -P, --print_reads      During progressive filter creation, print tagged reads\n"
		"                         to STDOUT in FASTQ format [disabled]\n"
		"  -r, --progressive=N    Progressive filter creation. The score threshold is\n"
		"                         specified by N, which may be either a floating point\n"
		"                         score between 0 and 1 or a positive integer.  If N is a\n"
		"                         positive integer, it is interpreted as the minimum\n"
		"                         number of contiguous matching bases required for a\n"
		"                         match.\n"
		"  -a, --streak=N         The number of hits tiling in second pass needed to jump\n"
		"                         Several tiles upon a miss. Progressive mode only. [3]\n"
		"  -l, --file_list=N      A file of list of file pairs to run in parallel.\n"
		"  -b, --baitScore=N      Score threshold when considering only bait. [r]\n"
		"  -e, --iterations=N     Pass through files N times if threshold is not met.\n"
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

	//command line variables
	string filterPrefix = "";
	string outputDir = "";
	unsigned kmerSize = 25;
	unsigned hashNum = 0;
	string subtractFilter = "";
	size_t entryNum = 0;
	bool printReads = false;
	double progressive = -1;
	bool inclusive = false;
	string fileListFilename = "";

	//long form arguments
	static struct option long_options[] = { { "fal_pos_rate", required_argument,
	NULL, 'f' }, { "file_prefix", required_argument, NULL, 'p' }, {
			"output_dir", required_argument, NULL, 'o' }, { "threads",
	required_argument, NULL, 't' }, { "inclusive", no_argument, NULL, 'i' }, {
			"version", no_argument, NULL, 'v' }, { "hash_num",
	required_argument, NULL, 'g' }, { "kmer_size", required_argument,
	NULL, 'k' }, { "subtract", required_argument, NULL, 's' }, { "num_ele",
	required_argument, NULL, 'n' }, { "interval",
	required_argument, NULL, 'I' }, { "file_list",
	required_argument, NULL, 'l' }, { "help", no_argument, NULL, 'h' }, {
			"print_reads", no_argument, NULL, 'P' }, { "progressive",
	required_argument, NULL, 'r' }, { "baitScore",
	required_argument, NULL, 'b' }, { "streak",
	required_argument, NULL, 'e' }, { "iterations",
	required_argument, NULL, 'a' }, { "no_rep_kmer",
	no_argument, NULL, 'd' }, {
	NULL, 0, NULL, 0 } };

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:p:o:k:n:g:hvs:n:t:Pr:ib:e:l:daI:",
			long_options, &option_index)) != -1) {
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
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				cerr << "Error - Invalid parameter! t: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'i': {
			inclusive = true;
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
		case 'I': {
			stringstream convert(optarg);
			if (!(convert >> opt::fileInterval)) {
				cerr << "Error - Invalid parameters! I: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'l': {
			stringstream convert(optarg);
			if (!(convert >> fileListFilename)) {
				cerr << "Error - Invalid set of bloom filter parameters! l: "
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
		case 'P': {
			printReads = true;
			break;
		}
		case 'd': {
			opt::noRep = true;
			break;
		}
		case 'r': {
			stringstream convert(optarg);
			unsigned matchLen;
			// if arg is a positive integer > 1, interpret as minimum match
			// length in bases
			if ((convert >> matchLen) && matchLen > 1) {
				progressive = matchLen;
				SeqEval::evalMode = SeqEval::EVAL_MIN_MATCH_LEN;
			} else {
				// not a positive integer > 1, so interpret as floating
				// point score between 0 and 1
				stringstream convert2(optarg);
				if (!(convert2 >> progressive)) {
					cerr
							<< "Error - Invalid set of bloom filter parameters! r: "
							<< optarg << endl;
					return 0;
				}
				if (progressive < 0 || progressive > 1) {
					cerr
							<< "Error - r must be a positive integer or a floating "
							<< "point between 0 and 1. Input given:" << optarg
							<< endl;
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
			if (opt::progItrns < 1) {
				cerr << "Error - e must be > 1" << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'a': {
			stringstream convert(optarg);
			if (!(convert >> opt::streakThreshold)) {
				cerr << "Error - Invalid parameter! a: " << optarg << endl;
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
		hashNum = BloomFilterInfo::calcOptimalHashNum(opt::fpr);
	}

	string file1 = "";
	string file2 = "";

	if (progressive != -1) {

		if (opt::baitThreshold == -1) {
			opt::baitThreshold = progressive;
		} else if ((opt::baitThreshold < 1 && progressive > 1)
				|| (opt::baitThreshold > 1 && progressive < 1)) {
			cerr
					<< "Both the bait threshold and progressive bloom filter must be either < 1 or > 1"
					<< endl;
			exit(1);
		}

		if (inputFiles.size() > 2 && fileListFilename == "") {
			file2 = inputFiles.back();
			inputFiles.pop_back();
			file1 = inputFiles.back();
			inputFiles.pop_back();
		} else if (fileListFilename == "") {
			cerr
					<< "at least 3 inputs are required when using progressive mode\nbiobloommaker <options> seed_seq reads1 reads2"
					<< endl;
			exit(1);
		}
		cerr << "Building Bloom filter in progressive mode. ";
		switch (SeqEval::evalMode) {
		case SeqEval::EVAL_MIN_MATCH_LEN:
			cerr << "Min match length = " << (unsigned) round(progressive)
					<< " bp" << endl;
			break;
		case SeqEval::EVAL_STANDARD:
		default:
			cerr << "Score threshold = " << progressive << endl;
			break;
		}
	}

	//create filter
	BloomFilterGenerator filterGen(inputFiles, kmerSize, hashNum, entryNum);

	if (entryNum == 0) {
		filterGen = BloomFilterGenerator(inputFiles, kmerSize, hashNum);
		entryNum = filterGen.getExpectedEntries();
	}

	BloomFilterInfo info(filterPrefix, kmerSize, hashNum, opt::fpr, entryNum,
			inputFiles);

	//get calculated size of Filter
	size_t filterSize = info.getCalcuatedFilterSize();
	cerr << "Allocating " << filterSize
			<< " bits of space for filter and will output filter this size (plus header)"
			<< endl;
	filterGen.setFilterSize(filterSize);

	size_t redundNum = 0;
	//output filter
	if (progressive != -1) {
		createMode mode = PROG_STD;
		if (inclusive) {
			mode = PROG_INC;
		}

		//load in file list
		if (fileListFilename != "") {
			string line;
			ifstream myfile(fileListFilename.c_str());
			if (myfile.is_open()) {
				getline(myfile, line);
				stringstream ss(line);
				string fileName1 = "";
				string fileName2 = "";
				ss >> fileName1;
				//Test file for pe or se reads
				if (ss >> fileName2) {
					opt::fileList1.push_back(fileName1);
					opt::fileList2.push_back(fileName2);
					while (getline(myfile, line)) {
						ss.str("");
						ss.clear();
						ss << line;
						ss >> fileName1;
						ss >> fileName2;
						opt::fileList1.push_back(fileName1);
						opt::fileList2.push_back(fileName2);
					}
					myfile.close();

					cerr << "Using file list paired mode" << endl;
					if (opt::baitThreshold == progressive) {
						redundNum = filterGen.generateProgressive(
								outputDir + filterPrefix + ".bf", progressive,
								opt::fileList1, opt::fileList2, mode,
								printReads, subtractFilter);
					} else {
						redundNum = filterGen.generateProgressiveBait(
								outputDir + filterPrefix + ".bf", progressive,
								opt::fileList1, opt::fileList2, mode,
								printReads, subtractFilter);
					}
				} else {
					opt::fileList1.push_back(fileName1);
					while (getline(myfile, line)) {
						ss.str("");
						ss.clear();
						ss << line;
						ss >> fileName1;
						opt::fileList1.push_back(fileName1);
					}
					myfile.close();

					cerr << "Using file list single end mode" << endl;
					if (opt::baitThreshold == progressive) {
						redundNum = filterGen.generateProgressive(
								outputDir + filterPrefix + ".bf", progressive,
								opt::fileList1, printReads, subtractFilter);
					} else {
						cerr
								<< "single end bait mode not implemented. If needed feature request cjustin@bcgsc.ca."
								<< endl;
						exit(1);
					}
				}
			} else
				cout << "Unable to open file";
		} else {
			if (opt::baitThreshold == progressive) {
				redundNum = filterGen.generateProgressive(
						outputDir + filterPrefix + ".bf", progressive, file1,
						file2, mode, printReads, subtractFilter);
			} else {
				redundNum = filterGen.generateProgressiveBait(
						outputDir + filterPrefix + ".bf", progressive, file1,
						file2, mode, printReads, subtractFilter);
			}
		}
	} else if (!subtractFilter.empty()) {
		redundNum = filterGen.generate(outputDir + filterPrefix + ".bf",
				subtractFilter);
	} else {
		redundNum = filterGen.generate(outputDir + filterPrefix + ".bf");
	}
	info.setTotalNum(filterGen.getTotalEntries());
	info.setRedundancy(redundNum);

	//code for redundancy checking
	//calculate redundancy rate
	double redunRate = double(redundNum) / double(filterGen.getTotalEntries())
			- info.getRedundancyFPR();
	if (redunRate > 0.25) {
		cerr
				<< "The ratio between redundant k-mers and unique k-mers is approximately: "
				<< redunRate << endl;
		cerr
				<< "Consider checking your files for duplicate sequences and adjusting them accordingly.\n"
						"High redundancy will cause filter sizes used overestimated, potentially resulting in a larger than needed filter.\n"
						"Alternatively you can set the number of elements wanted in the filter with (-n) and ignore this message."
				<< endl;
	}

	//output info
	info.printInfoFile(outputDir + filterPrefix + ".txt");
	cerr << "Filter Creation Complete." << endl;

	return 0;
}
