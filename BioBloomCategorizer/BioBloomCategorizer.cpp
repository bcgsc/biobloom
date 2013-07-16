/*
 * BioBloomCategorizer.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: cjustin
 */
//todo: UNIT TESTS!
#include <sstream>
#include <string>
#include <getopt.h>
#include "boost/unordered/unordered_map.hpp"
#include "Common/HashManager.h"
#include <vector>
#include <sys/stat.h>
#include "BioBloomClassifier.h"
#include "DataLayer/Options.h"
#include "config.h"

using namespace std;

#define PROGRAM "biobloomcategorizer"

void printVersion()
{
	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
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
	vector<string> currentInfoFile;
	string temp;
	stringstream converter(inputString);
	while (converter >> temp) {
		currentInfoFile.push_back(temp);
	}
	return currentInfoFile;
}

void folderCheck(const string &path)
{
	struct stat sb;

	if (stat(path.c_str(), &sb) == 0) {
		if (!S_ISDIR(sb.st_mode)) {
			cerr << "Error: Output folder - file exists with this name. "
					<< path << endl;
			exit(1);
		}
	} else {
		cerr << "Error: Output folder does not exist. " << path << endl;
		exit(1);
	}
}

void printHelpDialog()
{
	const char dialog[] =
			"Usage: biobloomcategorizer [OPTION]... -f \"[FILTER1]...\" [FILE]...\n"
					"Categorize Sequences. The input format may be FASTA, FASTQ, qseq, export, SAM or\n"
					"BAM format and compressed with gz, bz2 or xz and may be tarred.\n"
					"\n"
					"  -p, --prefix=N         Output prefix to use. Otherwise will output to current\n"
					"                         directory.\n"
					"  -f, --filter_files=N   List of filter files to use. Required option. \n"
					"                         eg. \"filter1.bf filter2.bf\"\n"
					"  -e, --paired_mode      Uses paired-end information. For BAM or SAM file if\n"
					"                         they are poorly ordered, memory usage will be much\n"
					"                         larger than normal. Sorting by read name may be needed.\n"
					"  -c, --counts=N         Outputs summary of raw counts of user specified hit\n"
					"                         counts to each filter of each read or read-pair. [0]\n"
					"  -g, --gz_output        Outputs all output files in compressed gzip.\n"
					"      --fa               Output categorized reads in Fasta files.\n"
					"      --fq               Output categorized reads in Fastq files.\n"
					"      --chastity         Discard and do not evaluate unchaste reads.\n"
					"      --no-chastity      Do not discard unchaste reads. [default]\n"
					"  -v  --version          Display version information.\n"
					"  -h, --help             Display this dialog.\n"
					"\nAdvanced options:\n"
					"  -t, --min_hit_thr=N    Minimum Hit Threshold Value. The absolute hit number\n"
					"                         needed for a hit to be considered a match. [2]\n"
					"  -m, --min_hit_pro=N    Minimum Hit Proportion Threshold Value. The proportion\n"
					"                         needed for a hit to be considered a match. [0.2]\n"
					"  -r, --redundant=N      The number of redundant tiles to use. Lowers effective\n"
					"                         false positive rate at the cost of increasing the\n"
					"                         effective kmer length by N. [0]\n"
					"\nOption presets:\n"
					"      --default          Run categorizer assuming default presets (ie. no\n"
					"                         advanced options toggled) [default]\n"
					"      --low_mem          Run categorizer assuming low memory presets.\n"
					"      --minimize_fpr     Run categorizer assuming minimized false positive rate\n"
					"                         presets.\n"
					"Report bugs to <cjustin@bcgsc.ca>.";
	cerr << dialog << endl;
	exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
	opt::chastityFilter = false;

	//switch statement variable
	int c;

	//control variables
	bool die = false;

	//command line variables
	string rawInputFiles = "";
	string outputPrefix = "";
	string filtersFile = "";
	string outputReadType = "";
	bool paired = false;
	size_t rawCounts = 0;

	int fastq = 0;
	int fasta = 0;
	string filePostfix = "";

	//advanced options
	size_t minHit = 2;
	double percentHit = 0.2;
	uint8_t tileModifier = 0;

	//preset options
	int defaultSettings = 0;
	int lowMem = 0;
	int minimize_fpr = 0;
	string presetType = "default";

	//long form arguments
	static struct option long_options[] = {
			{
					"prefix", optional_argument, NULL, 'p' }, {
					"min_hit_thr", optional_argument, NULL, 't' }, {
					"min_hit_pro", optional_argument, NULL, 'm' }, {
					"filter_files", required_argument, NULL, 'f' }, {
					"paired_mode", no_argument, NULL, 'e' }, {
					"counts", no_argument, NULL, 'c' }, {
					"help", no_argument, NULL, 'h' }, {
					"gz_output", no_argument, NULL, 'g' }, {
					"redundant", required_argument, NULL, 'r' }, {
					"chastity", no_argument, &opt::chastityFilter, 1 }, {
					"no-chastity", no_argument, &opt::chastityFilter, 0 }, {
					"fq", no_argument, &fastq, 0 }, {
					"fa", no_argument, &fasta, 0 }, {
					"default", no_argument, &defaultSettings, 0 }, {
					"low_mem", no_argument, &lowMem, 0 }, {
					"minimize_fpr", no_argument, &minimize_fpr, 0 }, {
					"version", no_argument, NULL, 0 }, {
					NULL, 0, NULL, 0 } };

	//check if only one preset was set
	if (defaultSettings || minimize_fpr || lowMem) {
		if (!(defaultSettings ^ minimize_fpr ^ lowMem)) {
			cerr << "Error: Cannot mix option presets" << endl;
			exit(1);
		}
	}

	//set presets

	if (lowMem) {
		tileModifier = 1;
		presetType = "low_mem";
	} else if (minimize_fpr) {
		tileModifier = 1;
		percentHit = 0.25;
		presetType = "minimum_fpr";
	}

	//actual checking step
	//Todo: add checks for duplicate options being set
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:t:m:p:hec:gr:v", long_options,
			&option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'm': {
			stringstream convert(optarg);
			if (!(convert >> percentHit)) {
				cerr << "Error - Invalid parameter! m: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			if (percentHit > 1) {
				cerr << "Error -m cannot be greater than 1 " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			presetType = "custom";
			break;
		}
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> minHit)) {
				cerr << "Error - Invalid parameter! t: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			presetType = "custom";
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
		case 'e': {
			paired = true;
			break;
		}
		case 'r': {
			stringstream convert(optarg);
			if (!(convert >> tileModifier)) {
				cerr << "Error - Invalid parameter! r: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			presetType = "custom";
			break;
		}
		case 'c': {
			stringstream convert(optarg);
			if (!(convert >> rawCounts)) {
				cerr << "Error - Invalid parameter! c: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'g': {
			filePostfix = ".gz";
			break;
		}
		case 'v': {
			printVersion();
			break;
		}
		case '?': {
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

	bool pairedBAMSAM = false;

	//check validity of inputs for paired end mode
	if (paired) {
		if (inputFiles.size() == 1
				&& (inputFiles[0].substr(inputFiles[0].size() - 4) != "bam"
						|| inputFiles[0].substr(inputFiles[0].size() - 4)
								!= "sam"))
		{
			pairedBAMSAM = true;
		} else if (inputFiles.size() == 2) {
			pairedBAMSAM = false;
		} else {
			cerr << "Usage of paired end mode:\n"
					<< "BioBloomCategorizer [OPTION]... -f \"[FILTER1]...\" [FILEPAIR1] [FILEPAIR2]\n"
					<< "or BioBloomCategorizer [OPTION]... -f \"[FILTER1]...\" [PAIREDBAMSAM]\n"
					<< endl;
			exit(1);
		}
	}

	//Check needed options
	if (inputFiles.size() == 0) {
		cerr << "Error: Need Input File" << endl;
		die = true;
	}
	if (filterFilePaths.size() == 0) {
		cerr << "Error: Need Filter File (-f)" << endl;
		die = true;
	}
	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}
	//check if output folder exists
	if (outputPrefix.find('/') != string::npos) {
		string tempStr = outputPrefix.substr(0, outputPrefix.find_last_of("/"));
		folderCheck(tempStr);
	}

	//set file output type
	if (fastq) {
		outputReadType = "fq";
	} else if (fasta) {
		outputReadType = "fa";
	}

	//load filters
	BioBloomClassifier BBC(filterFilePaths, minHit, percentHit, rawCounts,
			outputPrefix, filePostfix, tileModifier);

	//check filter preset type
	if (presetType != "custom") {
		if (!BBC.checkFilterPresetType(presetType)) {
			cerr
					<< "If you know what you are doing please ignore these warnings. Program will proceed."
					<< endl;
		}
	}

	//filtering step
	//create directory structure if it does not exist
	if (paired) {
		if (outputReadType != "") {
			if (pairedBAMSAM) {
				BBC.filterPairBAMPrint(inputFiles[0], outputReadType);
			} else {
				BBC.filterPairPrint(inputFiles[0], inputFiles[1],
						outputReadType);
			}
		} else {
			if (pairedBAMSAM) {
				BBC.filterPairBAM(inputFiles[0]);
			} else {
				BBC.filterPair(inputFiles[0], inputFiles[1]);
			}
		}
	} else {
		if (outputReadType != "") {
			BBC.filterPrint(inputFiles, outputReadType);
		} else {
			BBC.filter(inputFiles);
		}
	}
}
