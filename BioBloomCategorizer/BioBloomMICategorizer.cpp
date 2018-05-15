/*
 * BioBloomCategorizer.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: cjustin
 */
#include <sstream>
#include <string>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include "MIBFClassifier.hpp"
#include "config.h"
#include "Common/Options.h"
#include <zlib.h>
#include <fstream>
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "biobloommicategorizer"

void printVersion()
{
	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
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

/*
 * checks if file exists
 */
bool fexists(const string &filename) {
	ifstream ifile(filename.c_str());
	bool good = ifile.good();
	ifile.close();
	return good;
}


void printHelpDialog()
{
	const char dialog[] =
	"Usage: biobloommicategorizer [OPTION]... -f \"[FILTER1]...\" [FILE]...\n"
	"biobloomcategorizer [OPTION]... -e -f \"[FILTER1]...\" [FILE1.fq] [FILE2.fq]\n"
	"The input format may be FASTA, FASTQ, and compressed with gz.\n"
	"\n"
	"  -p, --prefix=N         Output prefix to use. Otherwise will output to current\n"
	"                         directory.\n"
	"  -f, --filter=N         Path of miBF '.bf' file\n"
	"  -e, --paired_mode      Uses paired-end information. For BAM or SAM files, if\n"
	"                         they are poorly ordered, the memory usage will be much\n"
	"                         larger than normal. Sorting by read name may be needed.\n"
	"  -s, --min_FPR=N        Minimum -10*log(FPR) threshold for a match. [60.0]\n"
	"  -t, --threads=N        The number of threads to use. [1]\n"
	"  -d, --stdout_filter    Outputs all matching reads to stdout in fastq\n"
	"                         and if paired will output will be interlaced.\n"
	"  -n, --inverse          Inverts the output of -d.\n"
	"      --fa               Output categorized reads in Fasta files.\n"
	"      --fq               Output categorized reads in Fastq files.\n"
	"      --version          Display version information.\n"
	"  -h, --help             Display this dialog.\n"
	"  -v, --verbose          Display verbose output\n"
	"  -I, --interval         the interval to report file processing status [10000000]\n"
	"Experimental options (may change in the future)\n"
	"  -a, --allowed_miss=N   Allowed misses in a bloom filter query, only works for\n"
	"                         miBFs.[0]\n"
	"  -m, --multi=N          Multi Match threshold. [3]\n"
	"  -r, --streak=N         The number of hits tiling in second pass needed to jump\n"
	"                         Several tiles upon a miss. Small values decrease\n"
	"                         runtime but decrease sensitivity. [3]\n"
	"  --debug                debug filter output mode.\n"
	"Report bugs to <cjustin@bcgsc.ca>.";

	cerr << dialog << endl;
	exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
	//switch statement variable
	int c;

	//control variables
	bool die = false;
	int FASTQ = 0;
	int FASTA = 0;
	int TSV = 0;
	int OPT_VERSION = 0;

	vector<string> inputFiles;

	//long form arguments
	static struct option long_options[] = { {
		"prefix", required_argument, NULL, 'p' }, {
		"filter", required_argument, NULL, 'f' }, {
		"paired", no_argument, NULL, 'e' }, {
		"min_FPR", required_argument, NULL, 's' }, {
		"help", no_argument, NULL, 'h' }, {
		"interval",	required_argument, NULL, 'I' }, {
		"threads", required_argument, NULL, 't' }, {
		"allowed_miss", required_argument, NULL, 'a' }, {
		"fq", no_argument, &FASTQ, 1 }, {
		"fa", no_argument, &FASTA, 1 }, {
		"tsv", no_argument, &TSV, 1 }, {
		"version", no_argument, &OPT_VERSION, 1 }, {
		"multi", required_argument, NULL, 'm' }, {
		"streak", required_argument, NULL, 'r' }, {
		"stdout_filter", no_argument, NULL, 'd' }, {
		"inverse", no_argument, NULL, 'n' }, {
		"debug", no_argument, &opt::debug, 1 }, {
		"verbose", no_argument, NULL, 'v' }, {
		NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "p:f:es:hI:t:a:m:r:dnv", long_options,
			&option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 's': {
			stringstream convert(optarg);
			if (!(convert >> opt::score)) {
				cerr << "Error - Invalid set of bloom filter parameters! s: "
						<< optarg << endl;
				exit(EXIT_FAILURE);
			}
			if (opt::score <= 0) {
				cerr << "Error - s must be a float > 0 Input given:" << optarg
						<< endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'f': {
			opt::filtersFile = optarg;
			break;
		}
		case 'p': {
			opt::outputPrefix = optarg;
			break;
		}
		case 'h': {
			printHelpDialog();
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
		case 'e': {
			opt::paired = true;
			break;
		}
		case 'a': {
			stringstream convert(optarg);
			if (!(convert >> opt::allowMisses)) {
				cerr << "Error - Invalid parameter! a: " << optarg << endl;
				exit(EXIT_FAILURE);
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
		case 'v': {
			opt::verbose++;
			break;
		}
		case 'r': {
			stringstream convert(optarg);
			if (!(convert >> opt::streakThreshold)) {
				cerr << "Error - Invalid parameter! r: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'm': {
			stringstream convert(optarg);
			if (!(convert >> opt::multiThresh)) {
				cerr << "Error - Invalid parameter! m: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'd': {
			opt::stdout = true;
			break;
		}
		case 'n': {
			opt::inverse = true;
			break;
		}
		case '?': {
			die = true;
			break;
		}
		}
	}

#if defined(_OPENMP)
	if (opt::threads > 0)
	omp_set_num_threads(opt::threads);
#endif

	if (OPT_VERSION) {
		printVersion();
	}

	while (optind < argc) {
		inputFiles.push_back(argv[optind]);
		optind++;
	}

	//check validity of inputs for paired end mode
	if (opt::paired && inputFiles.size() != 2) {
		cerr << "Usage of paired end mode:\n"
				<< "BioBloomCategorizer [OPTION]... -f \"[FILTER1]...\" [FILEPAIR1] [FILEPAIR2]\n"
				<< endl;
		exit(1);
	}

	//Check needed options
	if (inputFiles.size() == 0) {
		cerr << "Error: Need Input File" << endl;
		die = true;
	}
	if (opt::filtersFile.empty()) {
		cerr << "Error: Need Filter File (-f)" << endl;
		die = true;
	}
	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}
	//check if output folder exists
	if (opt::outputPrefix.find('/') != string::npos) {
		string tempStr = opt::outputPrefix.substr(0, opt::outputPrefix.find_last_of("/"));
		folderCheck(tempStr);
	}

	//set file output type
	if (FASTA && FASTQ) {
		cerr
				<< "Error: fasta (--fa) and fastq (--fq) outputs types cannot be both set"
				<< endl;
		exit(1);
	} else if (FASTQ) {
		opt::outputType = opt::FASTQ;
	} else if (FASTA) {
		opt::outputType = opt::FASTA;
	} else if (TSV) {
		opt::outputType = opt::TSV;
	}

	//check if files exist
	if (!fexists(opt::filtersFile)) {
		cerr << "Error: " + opt::filtersFile + " File cannot be opened" << endl;
		exit(1);
	}
	MIBFClassifier BMC(opt::filtersFile);
	if (opt::paired) {
		BMC.filterPair(inputFiles[0], inputFiles[1]);
	} else if (opt::debug) {
		BMC.filterOld(inputFiles);
	} else {
		BMC.filter(inputFiles);
	}
	return 0;
}
