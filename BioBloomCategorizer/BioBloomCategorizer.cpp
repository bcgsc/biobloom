/*
 * BioBloomCategorizer.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: cjustin
 */
//todo: refactor into additional class
//todo: UNIT TESTS!
#include <sstream>
#include <string>
#include <getopt.h>
#include "boost/unordered/unordered_map.hpp"
#include "Common/HashManager.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include "BioBloomClassifier.h"
using namespace std;

/*
 * Parses input string into seperate strings, returning a vector.
 */
vector<string> convertInputString(const string &inputString)
{
	vector<string> currentInfoFile;
	string temp;
	stringstream converter(inputString);
	while (converter >> temp) {
		currentInfoFile.push_back(temp);
	}
	return currentInfoFile;
}

void printHelpDialog()
{
	static const char dialog[] =
			"Usage: BioBloomCategorizer [OPTION]... -f \"[FILTER1], [FILTER2]...\" [FILE]...\n"
					"Categorize Sequences. The input format may be FASTA, FASTQ, qseq,\n"
					"export, SAM or BAM format and compressed with gz, bz2 or xz and\n"
					"may be tarred.\n"
					"\n"
					"  -p, --prefix=N         Output prefix to use. Otherwise will output\n"
					"                         to current directory.\n"
					"  -t, --min_hit_thr=N    Minimum Hit Threshold Value. Uses absolute\n"
					"                         hit number of read to categorize. [2]\n"
					"  -f, --filter_files=N   List of filter files to use.\n"
					"                         Eg. \"filter1.bf filter2.bf\"\n"
					"  -m, --min_hit_pro=N    Minimum Hit Proportion Threshold Value. Uses\n"
					"                         Proportion of hits to categorize. [0.25]\n"
					"  -o, --output_fastq     Output categorized reads in FastQ files.\n"
					"  -h, --help             Display this dialog."
					"\n"
					"Report bugs to <cjustin@bcgsc.ca>.\n";
	cerr << dialog << endl;
	exit(0);
}

int main(int argc, char *argv[])
{

	//switch statement variable
	int c;

	//command line variables
	string rawInputFiles = "";
	string outputPrefix = "";
	string filtersFile = "";
	int16_t minHit = 2;
	double percentHit = 0.25;
	bool printReads = false;
	bool die = false;

	//long form arguments
	static struct option long_options[] = { { "prefix", 0, NULL, 'p' }, {
			"min_hit_thr", 0, NULL, 't' }, { "min_hit_per", 0, NULL, 'm' }, {
			"output_fastq", 0, NULL, 'o' }, { "filter_files", 1, NULL, 'f' }, {
			"help", 0, NULL, 'h' }, { NULL, 0, NULL, 0 } };

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:t:om:p:h", long_options,
			&option_index)) != -1)
	{
		switch (c) {
		case 'm': {
			stringstream convert(optarg);
			if (!(convert >> percentHit)) {
				cerr << "Error - Invalid parameter! m: " << optarg << endl;
				return 0;
			}
			if (percentHit > 1) {
				cerr << "Error -m cannot be greater than 1 " << optarg << endl;
				return 0;
			}
			break;
		}
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> minHit)) {
				cerr << "Error - Invalid parameter! t: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'f': {
			filtersFile = optarg;
			break;
		}
		case 'p': {
			outputPrefix = optarg;
			break;
		}
		case 'h': {
			printHelpDialog();
			break;
		}
		case 'o': {
			printReads = true;
			break;
		}
		default: {
			die = true;
			break;
		}
		}
	}

	vector<string> filterFilePaths = convertInputString(filtersFile);
	vector<string> inputFiles = convertInputString(rawInputFiles);

	while (optind < argc) {
		inputFiles.push_back(argv[optind]);
		optind++;
	}

	//Check needed options
	if (inputFiles.size() == 0) {
		cerr << "Need Input File" << endl;
		die = true;
	}

	if (filterFilePaths.size() == 0) {
		cerr << "Need Filter File (-f)" << endl;
		die = true;
	}
	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	//load filters
	BioBloomClassifier BBC(filterFilePaths, minHit, percentHit);

	//filtering step
	//create directory structure if it does not exist
	BBC.filterPair(inputFiles[0], inputFiles[1], outputPrefix);
//	if (printReads) {
//		BBC.filterPrintReads(inputFiles, outputPrefix);
//	} else {
//		BBC.filter(inputFiles, outputPrefix);
//	}
}
