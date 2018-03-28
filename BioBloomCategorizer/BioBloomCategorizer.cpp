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
#include "BioBloomClassifier.h"
#include "MIBFClassifier.hpp"
#include "config.h"
#include "Common/Options.h"
#include "Common/SeqEval.h"
#include <zlib.h>
#include <fstream>
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "biobloomcategorizer"

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
	"Usage: biobloomcategorizer [OPTION]... -f \"[FILTER1]...\" [FILE]...\n"
	"biobloomcategorizer [OPTION]... -e -f \"[FILTER1]...\" [FILE1.fq] [FILE2.fq]\n"
	"Categorize Sequences. The input format may be FASTA, FASTQ, qseq, export, SAM or\n"
	"BAM format and compressed with gz, bz2 or xz and may be tarred.\n"
	"\n"
	"  -p, --prefix=N         Output prefix to use. Otherwise will output to current\n"
	"                         directory.\n"
	"  -f, --filter_files=N   List of filter files to use. Required option. \n"
	"                         eg. \"filter1.bf filter2.bf\"\n"
	"  -e, --paired_mode      Uses paired-end information. For BAM or SAM files, if\n"
	"                         they are poorly ordered, the memory usage will be much\n"
	"                         larger than normal. Sorting by read name may be needed.\n"
	"  -i, --inclusive        If one paired read matches, both reads will be included\n"
	"                         in the filter. \n"
	"  -s, --score=N          Score threshold for matching. N may be either a\n"
	"                         floating point score between 0 and 1 or a positive\n"
	"                         integer representing the minimum match length in bases.\n"
	"                         If N is a floating point, the maximum threshold is any \n"
	"                         number less than 1, and the minimum is 0 (highest\n"
	"                         sensitivity). If set to 1, the best hit is used rather\n"
	"                         than the threshold and the score will be appended to the\n"
	"                         header of the output read. [0.15]\n"
	"  -w, --with_score       Output multimatches with scores in the order of filter.\n"
	"  -t, --threads=N        The number of threads to use. [1]\n"
	"  -g, --gz_output        Outputs all output files in compressed gzip.\n"
	"      --fa               Output categorized reads in Fasta files.\n"
	"      --fq               Output categorized reads in Fastq files.\n"
	"      --chastity         Discard and do not evaluate unchaste reads.\n"
	"      --no-chastity      Do not discard unchaste reads. [default]\n"
	"  -l, --file_list=N      A file of list of file pairs to run in parallel.\n"
	"  -v  --version          Display version information.\n"
	"  -h, --help             Display this dialog.\n"
	"  -I, --interval         the interval to report file processing status [10000000]\n"
	"Advanced options:\n"
//	"  -m, --min_hit=N        Minimum Hit Threshold Value. The absolute hit number\n"
//	"                         needed over initial tiling of read to continue. Higher\n"
//	"                         values decrease runtime but lower sensitivity.[0]\n"
	"  -r, --streak=N         The number of hits tiling in second pass needed to jump\n"
	"                         Several tiles upon a miss. Small values decrease\n"
	"                         runtime but decrease sensitivity. [3]\n"
//	"  -o, --min_hit_only     Use only initial pass filtering to evaluate reads. Fast\n"
//	"                         but low specificity, use only on long reads (>100bp).\n"
	"  -c, --ordered          Use ordered filtering. Order of filters matters\n"
	"                         (filters listed first have higher priority). Only taken\n"
	"                         advantage of when k-mer sizes and number of hash\n"
	"                         functions are the same.\n"
	"  -d, --stdout_filter    Outputs all matching reads to stdout for the first\n"
	"                         filter listed by -f. Reads are outputed in fastq,\n"
	"                         and if paired will output will be interlaced.\n"
	"Options for multi index bloom filters:\n"
	"  -a, --allowed_miss=N   Allowed misses in a bloom filter query, only works for\n"
	"                         miBFs.[0]\n"
	"  --debug                debug filter output mode."
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

	//command line variables
	string rawInputFiles = "";
	string filtersFile = "";
	bool paired = false;

	int fastq = 0;
	int fasta = 0;
	bool withScore = false;
	bool stdout = false;

	//advanced options
	bool collab = false;

	string fileListFilename = "";

	//long form arguments
	static struct option long_options[] = { {
		"prefix", required_argument, NULL, 'p' }, {
		"filter_files", required_argument, NULL, 'f' }, {
		"paired_mode", no_argument, NULL, 'e' }, {
		"inclusive", no_argument, NULL, 'i' }, {
		"score", required_argument, NULL, 's' }, {
		"help", no_argument, NULL, 'h' }, {
		"interval",	required_argument, NULL, 'I' }, {
		"threads", required_argument, NULL, 't' }, {
		"allowed_miss", required_argument, NULL, 'a' }, {
		"gz_output", no_argument, NULL, 'g' }, {
		"fq", no_argument, &fastq, 1 }, {
		"fa", no_argument, &fasta, 1 }, {
		"file_list", required_argument, NULL, 'l' }, {
		"version", no_argument, NULL, 'v' }, {
		"min_hit_thr", required_argument, NULL, 'm' }, {
		"streak", required_argument, NULL, 'r' }, {
		"min_hit_only", no_argument, NULL, 'o' }, {
		"ordered", no_argument, NULL, 'c' }, {
		"stdout_filter", no_argument, NULL, 'd' }, {
		"with_score", no_argument, NULL, 'w' }, {
		"debug", no_argument, &opt::debug, 1 }, {
		"verbose", no_argument, &opt::verbose, 1 }, {
		NULL, 0, NULL, 0 } };

	//actual checking step
	//Todo: add checks for duplicate options being set
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:m:p:hegl:vs:r:t:cdiwI:a:G:", long_options,
			&option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'm': {
			stringstream convert(optarg);
			if (!(convert >> opt::minHit)) {
				cerr << "Error - Invalid parameter! m: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 's': {
			stringstream convert(optarg);
			unsigned matchLen;
			// if arg is a positive integer > 1, interpret as minimum match
			// length in bases
			if ((convert >> matchLen) && matchLen > 1) {
				opt::score = (unsigned)matchLen;
				SeqEval::evalMode = SeqEval::EVAL_MIN_MATCH_LEN;
				cerr << "Min match length threshold: " << matchLen
					<<" bp" << endl;
			} else {
				// not a positive integer > 1, so interpret as floating
				// point score between 0 and 1
				stringstream convert2(optarg);
				if (!(convert2 >> opt::score)) {
					cerr << "Error - Invalid set of bloom filter parameters! s: "
						<< optarg << endl;
					exit(EXIT_FAILURE);
				}
				if (opt::score < 0 || opt::score > 1) {
					cerr << "Error - s must be a positive integer or a floating "
						<< "point between 0 and 1. Input given:"
						<< optarg << endl;
					exit(EXIT_FAILURE);
				}
				SeqEval::evalMode = SeqEval::EVAL_STANDARD;
				if (opt::score == 1)
					cerr << "Running in 'best match' mode (no score threshold)"
						<< endl;
				else
					cerr << "Min score threshold: " << opt::score << endl;
			}
			break;
		}
		case 'f': {
			filtersFile = optarg;
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
			paired = true;
			break;
		}
		case 'i': {
			opt::inclusive = true;
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
		case 'G': {
			stringstream convert(optarg);
			if (!(convert >> opt::maxGroupSize)) {
				cerr << "Error - Invalid parameter! G: " << optarg << endl;
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
		case 'g': {
			opt::filePostfix = ".gz";
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
		case 'v': {
			printVersion();
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
		case 'c': {
			collab = true;
			break;
		}
		case 'd': {
			stdout = true;
			break;
		}
		case 'w': {
			withScore = true;
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
				&& (inputFiles[0].substr(inputFiles[0].size() - 4) == ".bam"
						|| inputFiles[0].substr(inputFiles[0].size() - 4)
								== ".sam"))
		{
			pairedBAMSAM = true;
		} else if (inputFiles.size() == 2) {
			pairedBAMSAM = false;
		}
		else if (fileListFilename == "") {
			cerr << "Usage of paired end mode:\n"
					<< "BioBloomCategorizer [OPTION]... -f \"[FILTER1]...\" [FILEPAIR1] [FILEPAIR2]\n"
					<< "or BioBloomCategorizer [OPTION]... -f \"[FILTER1]...\" [PAIREDBAMSAM]\n"
					<< endl;
			exit(1);
		}
	}

	//Check needed options
	if (inputFiles.size() == 0 && fileListFilename == "") {
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
	if (opt::outputPrefix.find('/') != string::npos) {
		string tempStr = opt::outputPrefix.substr(0, opt::outputPrefix.find_last_of("/"));
		folderCheck(tempStr);
	}

	//set file output type
	if (fastq && fasta) {
		cerr
				<< "Error: fasta (--fa) and fastq (--fq) outputs types cannot be both set"
				<< endl;
		exit(1);
	} else if (fastq) {
		opt::outputType = "fq";
	} else if (fasta) {
		opt::outputType = "fa";
	}

	//-w option cannot be used without output method
	if (withScore && opt::outputType == "") {
		cerr << "Error: -w option cannot be used without output method" << endl;
		exit(1);
	}

	// -w option cannot be used with min match length arg to -s
	if (withScore && SeqEval::evalMode == SeqEval::EVAL_MIN_MATCH_LEN) {
		cerr << "Error: -w option is not supported when using "
			<< "minimum match length (positive integer arg "
			<< "to -s)" << endl;
		exit(1);
	}

	//check filters
	for (vector<string>::const_iterator it = filterFilePaths.begin();
			it != filterFilePaths.end(); ++it)
	{
		//check if files exist
		if (!fexists(*it)) {
			cerr << "Error: " + (*it) + " File cannot be opened" << endl;
			exit(1);
		}
		string infoFileName = (*it).substr(0, (*it).length() - 2) + "txt";
		//TODO: error if filter types are mixed
		if (!fexists(infoFileName) && opt::filterType != BLOOMMAP) {
			cerr << "Multi Index Mode detected from lack of filter.txt files" << endl;
			opt::filterType = BLOOMMAP;
		}
	}

	if(opt::filterType == BLOOMMAP){
		if(opt::outputType == ""){
			opt::outputType = "fa";
		}
		MIBFClassifier BMC(filterFilePaths[0]);
		if (paired) {
			cerr << "paired MIBF classification not supported yet" << endl;
			exit(1);
//			BMC.filterPair(inputFiles[0], inputFiles[1]);
		} else if (opt::debug) {
			BMC.filterOld(inputFiles);
		}
		else{
			BMC.filter(inputFiles);
		}
		//load filters
		return 0;
	}

	//load filters
	BioBloomClassifier bbc(filterFilePaths, opt::score, opt::outputPrefix, opt::filePostfix,
			withScore);

	if (stdout) {
		bbc.setStdout();
	}

	if (collab && opt::minHit) {
		cerr << "Error: -m -c outputs types cannot be both set" << endl;
		exit(1);
	} else if (collab) {
		bbc.setOrderedFilter();
	}

	//filtering step
	//create directory structure if it does not exist
	if (paired) {
		if (opt::inclusive) {
			bbc.setInclusive();
		}
		if (opt::outputType != "") {
			if (pairedBAMSAM) {
				bbc.filterPairPrint(inputFiles[0], opt::outputType);
			} else {
				bbc.filterPairPrint(inputFiles[0], inputFiles[1],
						opt::outputType);
			}
		} else {
			if (pairedBAMSAM) {
				bbc.filterPair(inputFiles[0]);
			} else if (fileListFilename != "") {
				string line;
				ifstream myfile(fileListFilename.c_str());
				if (myfile.is_open()) {
					while (getline(myfile, line)) {
						stringstream ss(line);
						string fileName1 = "";
						string fileName2 = "";
						ss >> fileName1;
						ss >> fileName2;
						opt::fileList1.push_back(fileName1);
						opt::fileList2.push_back(fileName2);
					}
					myfile.close();
					cerr << "Using file list" << endl;
					bbc.filterPair(opt::fileList1, opt::fileList2);
				} else
					cout << "Unable to open file";
			}
			else{
				bbc.filterPair(inputFiles[0], inputFiles[1]);
			}
		}
	} else {
		if (opt::outputType != "") {
			bbc.filterPrint(inputFiles, opt::outputType);
		} else {
			bbc.filter(inputFiles);
		}
	}
}
