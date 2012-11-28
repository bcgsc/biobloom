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

int main(int argc, char *argv[])
{

	//switch statement variable
	int c;

	//command line variables
	string rawInputFiles = "";
	string outputDir = "";
	string filtersFile = "";
	int16_t minHit = 2;
	bool printReads = false;

	//long form arguments
	//each option format { "optionName", necessary option or not, I have no idea, 'symbol'}
	static struct option long_options[] = { { "input_files", 1, NULL, 'i' }, {
			"output_dir", 0, NULL, 'o' }, { "filters", 1, NULL, 'f' }, {
			"min_hit", 0, NULL, 'm' }, { "print_fasta", 1, NULL, 'p' }, { NULL,
			0, NULL, 0 } };

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "i:f:o:k:m:p", long_options,
			&option_index)) != -1)
	{
		switch (c) {
		case 'i': {
			rawInputFiles = optarg;
			break;
		}
		case 'f': {
			filtersFile = optarg;
			break;
		}
		case 'o': {
			outputDir = optarg;
			if (outputDir.at(outputDir.length() - 1) != '/') {
				outputDir = outputDir + '/';
			}
			break;
		}
		case 'm': {
			stringstream convert(optarg);
			if (!(convert >> minHit)) {
				cout << "Error - Invalid parameter! m: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'p': {
			printReads = true;
			break;
		}
		default: {
			cout << (char) c << " - Option not recognized" << endl;
			exit(1);
		}
		}
	}

	//Check needed options
	if (rawInputFiles == "") {
		cout << "Need Input File (-i)" << endl;
		exit(1);
	}

	if (filtersFile == "") {
		cout << "Need Filter File (-i)" << endl;
		exit(1);
	}

	vector<string> filterFilePaths = convertInputString(filtersFile);
	vector<string> inputFiles = convertInputString(rawInputFiles);

	//load filters
	BioBloomClassifier BBC(filterFilePaths);

	//filtering step
	//create directory structure if it does not exist
	if (printReads) {
		BBC.filterPrintReads(inputFiles, outputDir);
	} else {
		BBC.filter(inputFiles, outputDir);
	}
}
