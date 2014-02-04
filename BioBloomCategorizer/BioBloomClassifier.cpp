/*
 * BioBloomClassifier.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: cjustin
 */

#include "BioBloomClassifier.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "Common/Dynamicofstream.h"
#include "ResultsManager.h"
#if _OPENMP
# include <omp.h>
#endif

BioBloomClassifier::BioBloomClassifier(const vector<string> &filterFilePaths,
		double scoreThreshold, const string &prefix,
		const string &outputPostFix, uint16_t streakThreshold, uint16_t minHit,
		bool minHitOnly) :
		scoreThreshold(scoreThreshold), filterNum(filterFilePaths.size()), prefix(
				prefix), postfix(outputPostFix), streakThreshold(
				streakThreshold), minHit(minHit), minHitOnly(minHitOnly), noMatch(
				"noMatch"), multiMatch("multiMatch")
{
	loadFilters(filterFilePaths);
}

/*
 * Generic filtering function (single end, no fa or fq file outputs)
 */
void BioBloomClassifier::filter(const vector<string> &inputFiles)
{

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, scoreThreshold);

	size_t totalReads = 0;

	//print out header info and initialize variables

	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it)
	{
		FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
		for (FastqRecord rec;;) {
			bool good;
#pragma omp critical(sequence)
			{
				good = sequence >> rec;
				//track read progress
			}
			if (good) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 1000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				unordered_map<string, bool> hits(filterNum);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					if (minHitOnly) {
						evaluateRead(rec, *j, hits);
					} else {
						evaluateReadStd(rec, *j, hits);
					}
				}

				//Evaluate hit data and record for summary
				resSummary.updateSummaryData(hits);
			} else
				break;
		}
		assert(sequence.eof());
	}

	cerr << "Total Reads:" << totalReads << endl;

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

/*
 * Filters reads
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * Prints reads into seperate files
 *
 * outputType must be fa or fq
 */
void BioBloomClassifier::filterPrint(const vector<string> &inputFiles,
		const string &outputType)
{

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, scoreThreshold);

	size_t totalReads = 0;

	unordered_map<string, shared_ptr<Dynamicofstream> > outputFiles;
	shared_ptr<Dynamicofstream> no_match(
			new Dynamicofstream(
					prefix + "_" + noMatch + "." + outputType + postfix));
	shared_ptr<Dynamicofstream> multi_match(
			new Dynamicofstream(
					prefix + "_" + multiMatch + "." + outputType + postfix));
	outputFiles[noMatch] = no_match;
	outputFiles[multiMatch] = multi_match;

	//initialize variables
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			shared_ptr<Dynamicofstream> temp(
					new Dynamicofstream(
							prefix + "_" + *i + "." + outputType + postfix));
			outputFiles[*i] = temp;
		}
	}

	//print out header info and initialize variables

	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it)
	{
		FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
		for (FastqRecord rec;;) {
			bool good;
#pragma omp critical(sequence)
			{
				good = sequence >> rec;
			}
			if (good) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 1000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				unordered_map<string, bool> hits(filterNum);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					if (minHitOnly) {
						evaluateRead(rec, *j, hits);
					} else {
						evaluateReadStd(rec, *j, hits);
					}
				}

				//Evaluate hit data and record for summary
				const string &outputFileName = resSummary.updateSummaryData(
						hits);
				if (outputType == "fa") {
#pragma omp critical(outputFiles)
					{
						(*outputFiles[outputFileName]) << ">" << rec.id << "\n"
								<< rec.seq << "\n";
					}
				} else {
#pragma omp critical(outputFiles)
					{
						(*outputFiles[outputFileName]) << "@" << rec.id << "\n"
								<< rec.seq << "\n+\n" << rec.qual << "\n";
					}
				}
			} else
				break;
		}
		assert(sequence.eof());
	}

	//close sorting files
	for (unordered_map<string, shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j)
	{
		j->second->close();
	}
	cerr << "Total Reads:" << totalReads << endl;

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPair(const string &file1, const string &file2)
{

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, scoreThreshold);

	size_t totalReads = 0;

	cerr << "Filtering Start" << "\n";

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
	for (FastqRecord rec1;;) {
		FastqRecord rec2;
		bool good1;
		bool good2;

#pragma omp critical(sequence1)
		{
			good1 = sequence1 >> rec1;
			good2 = sequence2 >> rec2;
			//track read progress
		}

		if (good1 && good2) {
#pragma omp critical(totalReads)
			{
				++totalReads;
				if (totalReads % 1000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}
			}

			//hits results stored in hashmap of filter names and hits
			unordered_map<string, bool> hits1(filterNum);
			unordered_map<string, bool> hits2(filterNum);

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
				string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
				if (tempStr1 == tempStr2) {
					if (minHitOnly) {
						evaluateRead(rec1, *j, hits1);
						evaluateRead(rec2, *j, hits2);
					} else {
						evaluateReadStd(rec1, *j, hits1);
						evaluateReadStd(rec2, *j, hits2);
					}
				} else {
					cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
							<< tempStr2 << endl;
					exit(1);
				}
			}

			string readID = rec1.id.substr(0, rec1.id.length() - 2);

			//Evaluate hit data and record for summary
			resSummary.updateSummaryData(hits1, hits2);
		} else
			break;
	}
	if (!sequence1.eof() || !sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
	}

	cerr << "Total Reads:" << totalReads << endl;

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * prints reads
 */
void BioBloomClassifier::filterPairPrint(const string &file1,
		const string &file2, const string &outputType)
{

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, scoreThreshold);

	size_t totalReads = 0;

	unordered_map<string, shared_ptr<Dynamicofstream> > outputFiles;
	shared_ptr<Dynamicofstream> noMatch1(
			new Dynamicofstream(
					prefix + "_" + noMatch + "_1." + outputType + postfix));
	shared_ptr<Dynamicofstream> noMatch2(
			new Dynamicofstream(
					prefix + "_" + noMatch + "_2." + outputType + postfix));
	shared_ptr<Dynamicofstream> multiMatch1(
			new Dynamicofstream(
					prefix + "_" + multiMatch + "_1." + outputType + postfix));
	shared_ptr<Dynamicofstream> multiMatch2(
			new Dynamicofstream(
					prefix + "_" + multiMatch + "_2." + outputType + postfix));
	outputFiles[noMatch + "1"] = noMatch1;
	outputFiles[noMatch + "2"] = noMatch2;
	outputFiles[multiMatch + "1"] = multiMatch1;
	outputFiles[multiMatch + "2"] = multiMatch2;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			shared_ptr<Dynamicofstream> temp1(
					new Dynamicofstream(
							prefix + "_" + *i + "_1." + outputType + postfix));
			shared_ptr<Dynamicofstream> temp2(
					new Dynamicofstream(
							prefix + "_" + *i + "_2." + outputType + postfix));
			outputFiles[*i + "1"] = temp1;
			outputFiles[*i + "2"] = temp2;
		}
	}

	//for output files in consistent order

	cerr << "Filtering Start" << "\n";

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
	for (FastqRecord rec1;;) {
		FastqRecord rec2;
		bool good1;
		bool good2;

#pragma omp critical(sequence1)
		{
			good1 = sequence1 >> rec1;
			good2 = sequence2 >> rec2;
		}

		if (good1 && good2) {
#pragma omp critical(totalReads)
			{
				++totalReads;
				if (totalReads % 1000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}
			}

			//hits results stored in hashmap of filter names and hits
			unordered_map<string, bool> hits1(filterNum);
			unordered_map<string, bool> hits2(filterNum);

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
				string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
				if (tempStr1 == tempStr2) {
					if (minHitOnly) {
						evaluateRead(rec1, *j, hits1);
						evaluateRead(rec2, *j, hits2);
					} else {
						evaluateReadStd(rec1, *j, hits1);
						evaluateReadStd(rec2, *j, hits2);
					}
				} else {
					cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
							<< tempStr2 << endl;
					exit(1);
				}
			}

			string readID = rec1.id.substr(0, rec1.id.length() - 2);

			//Evaluate hit data and record for summary

			const string &outputFileName = resSummary.updateSummaryData(hits1,
					hits2);
#pragma omp critical(outputFiles)
			{
				if (outputType == "fa") {
					(*outputFiles[outputFileName + "1"]) << ">" << rec1.id
							<< "\n" << rec1.seq << "\n";
					(*outputFiles[outputFileName + "2"]) << ">" << rec2.id
							<< "\n" << rec2.seq << "\n";
				} else {
					(*outputFiles[outputFileName + "1"]) << "@" << rec1.id
							<< "\n" << rec1.seq << "\n+\n" << rec1.qual << "\n";
					(*outputFiles[outputFileName + "2"]) << "@" << rec2.id
							<< "\n" << rec2.seq << "\n+\n" << rec2.qual << "\n";
				}
			}
		} else
			break;
	}
	if (!sequence1.eof() || !sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
	}

	//close sorting files
	for (unordered_map<string, shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j)
	{
		j->second->close();
	}

	cerr << "Total Reads:" << totalReads << endl;

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPairBAM(const string &file)
{

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, scoreThreshold);

	unordered_map<string, FastqRecord> unPairedReads;

	size_t totalReads = 0;

	//print out header info and initialize variables for summary

	cerr << "Filtering Start" << "\n";

	FastaReader sequence(file.c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
	for (FastqRecord rec;;) {
		bool good;
		bool pairfound;
		FastqRecord rec1;
		FastqRecord rec2;
#pragma omp critical(unPairedReads)
		{
			good = sequence >> rec;
			//track read progress
			if (good) {
				string readID = rec.id.substr(0, rec.id.length() - 2);
				if (unPairedReads.find(readID) != unPairedReads.end()) {
					rec1 = rec.id.at(rec.id.length() - 1) == '1' ?
							rec : unPairedReads[readID];
					rec2 = rec.id.at(rec.id.length() - 1) == '2' ?
							rec : unPairedReads[readID];
					pairfound = true;
					unPairedReads.erase(readID);
				} else {
					pairfound = false;
					unPairedReads[readID] = rec;
				}
			}
		}

		if (good) {
			if (pairfound) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 1000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				unordered_map<string, bool> hits1(filterNum);
				unordered_map<string, bool> hits2(filterNum);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					if (minHitOnly) {
						evaluateRead(rec1, *j, hits1);
						evaluateRead(rec2, *j, hits2);
					} else {
						evaluateReadStd(rec1, *j, hits1);
						evaluateReadStd(rec2, *j, hits2);
					}
				}

				//Evaluate hit data and record for summary
				resSummary.updateSummaryData(hits1, hits2);
			}
		} else
			break;
	}
	assert(sequence.eof());

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * Prints reads into separate files
 */
void BioBloomClassifier::filterPairBAMPrint(const string &file,
		const string &outputType)
{

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, scoreThreshold);

	unordered_map<string, FastqRecord> unPairedReads;

	size_t totalReads = 0;

	unordered_map<string, shared_ptr<Dynamicofstream> > outputFiles;
	shared_ptr<Dynamicofstream> noMatch1(
			new Dynamicofstream(
					prefix + "_" + noMatch + "_1." + outputType + postfix));
	shared_ptr<Dynamicofstream> noMatch2(
			new Dynamicofstream(
					prefix + "_" + noMatch + "_2." + outputType + postfix));
	shared_ptr<Dynamicofstream> multiMatch1(
			new Dynamicofstream(
					prefix + "_" + multiMatch + "_1." + outputType + postfix));
	shared_ptr<Dynamicofstream> multiMatch2(
			new Dynamicofstream(
					prefix + "_" + multiMatch + "_2." + outputType + postfix));
	outputFiles[noMatch + "1"] = noMatch1;
	outputFiles[noMatch + "2"] = noMatch2;
	outputFiles[multiMatch + "1"] = multiMatch1;
	outputFiles[multiMatch + "2"] = multiMatch2;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			shared_ptr<Dynamicofstream> temp1(
					new Dynamicofstream(
							prefix + "_" + *i + "_1." + outputType + postfix));
			shared_ptr<Dynamicofstream> temp2(
					new Dynamicofstream(
							prefix + "_" + *i + "_2." + outputType + postfix));
			outputFiles[*i + "1"] = temp1;
			outputFiles[*i + "2"] = temp2;
		}
	}

	//print out header info and initialize variables for summary
	cerr << "Filtering Start" << "\n";

	FastaReader sequence(file.c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
	for (FastqRecord rec;;) {
		bool good;
		bool pairfound = false;
		FastqRecord rec1;
		FastqRecord rec2;
#pragma omp critical(unPairedReads)
		{
			good = sequence >> rec;
			//track read progress
			if (good) {
				string readID = rec.id.substr(0, rec.id.length() - 2);
				if (unPairedReads.find(readID) != unPairedReads.end()) {
					rec1 = rec.id.at(rec.id.length() - 1) == '1' ?
							rec : unPairedReads[readID];
					rec2 = rec.id.at(rec.id.length() - 1) == '2' ?
							rec : unPairedReads[readID];
					pairfound = true;
					unPairedReads.erase(readID);
				} else {
					unPairedReads[readID] = rec;
				}
			}
		}
		if (good) {
			if (pairfound) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 1000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				unordered_map<string, bool> hits1(filterNum);
				unordered_map<string, bool> hits2(filterNum);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					string tempStr1 = rec1.id.substr(0,
							rec1.id.find_last_of("/"));
					string tempStr2 = rec2.id.substr(0,
							rec2.id.find_last_of("/"));
					if (tempStr1 == tempStr2) {
						if (minHitOnly) {
							evaluateRead(rec1, *j, hits1);
							evaluateRead(rec2, *j, hits2);
						} else {
							evaluateReadStd(rec1, *j, hits1);
							evaluateReadStd(rec2, *j, hits2);
						}
					} else {
						cerr << "Read IDs do not match" << "\n" << tempStr1
								<< "\n" << tempStr2 << endl;
						exit(1);
					}
				}

				//Evaluate hit data and record for summary
				const string &outputFileName = resSummary.updateSummaryData(
						hits1, hits2);
#pragma omp critical(outputFiles)
				{

					if (outputType == "fa") {
						(*outputFiles[outputFileName + "1"]) << ">" << rec1.id
								<< "\n" << rec1.seq << "\n";
						(*outputFiles[outputFileName + "2"]) << ">" << rec2.id
								<< "\n" << rec2.seq << "\n";
					} else {
						(*outputFiles[outputFileName + "1"]) << "@" << rec1.id
								<< "\n" << rec1.seq << "\n+\n" << rec1.qual
								<< "\n";
						(*outputFiles[outputFileName + "2"]) << "@" << rec2.id
								<< "\n" << rec2.seq << "\n+\n" << rec2.qual
								<< "\n";
					}
				}
			}
		} else
			break;
	}
	assert(sequence.eof());

	//close sorting files
	for (unordered_map<string, shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j)
	{
		j->second->close();
	}

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

//helper methods

/*
 * Loads list of filters into memory
 * todo: Implement non-block I/O when loading multiple filters at once
 */
void BioBloomClassifier::loadFilters(const vector<string> &filterFilePaths)
{
	cerr << "Starting to Load Filters." << endl;
	//load up files
	for (vector<string>::const_iterator it = filterFilePaths.begin();
			it != filterFilePaths.end(); ++it)
	{
		//check if files exist
		if (!fexists(*it)) {
			cerr << "Error: " + (*it) + " File cannot be opened" << endl;
			exit(1);
		}
		string infoFileName = (*it).substr(0, (*it).length() - 2) + "txt";
		if (!fexists(infoFileName)) {
			cerr
					<< "Error: " + (infoFileName)
							+ " File cannot be opened. A corresponding info file is needed."
					<< endl;
			exit(1);
		}

		//info file creation
		shared_ptr<BloomFilterInfo> info(new BloomFilterInfo(infoFileName));
		//append kmer size to hash signature to insure correct kmer size is used
		stringstream hashSig;
		hashSig << info->getHashNum() << info->getKmerSize();

		//if hashSig exists add filter to list
		if (infoFiles.count(hashSig.str()) != 1) {
			hashSigs.push_back(hashSig.str());
			vector<shared_ptr<BloomFilterInfo> > tempVect;
			shared_ptr<MultiFilter> temp(
					new MultiFilter(info->getHashNum(), info->getKmerSize()));
			filters[hashSig.str()] = temp;
			infoFiles[hashSig.str()] = tempVect;
		}
		infoFiles[hashSig.str()].push_back(info);
		boost::shared_ptr<BloomFilter> filter(
				new BloomFilter(info->getCalcuatedFilterSize(),
						info->getHashNum(), info->getKmerSize(), *it));
		filters[hashSig.str()]->addFilter(info->getFilterID(), filter);
		filtersSingle[info->getFilterID()] = filter;
		cerr << "Loaded Filter: " + info->getFilterID() << endl;
	}
	cerr << "Filter Loading Complete." << endl;
}

/*
 * checks if file exists
 */
const bool BioBloomClassifier::fexists(const string &filename) const
{
	ifstream ifile(filename.c_str());
	return ifile;
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Updates hits value to number of hits (hashSig is used to as key)
 * Faster variant that assume there a redundant tile of 0
 */
void BioBloomClassifier::evaluateRead(const FastqRecord &rec,
		const string &hashSig, unordered_map<string, bool> &hits)
{
	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*filters[hashSig]).getFilterIds();

	//get kmersize for set of info files
	uint16_t kmerSize = infoFiles.at(hashSig).front()->getKmerSize();

	unordered_map<string, uint16_t> tempHits;

	//Establish tiling pattern
	uint16_t startModifier1 = (rec.seq.length() % kmerSize) / 2;
	size_t currentKmerNum = 0;

	for (vector<string>::const_iterator i = idsInFilter.begin();
			i != idsInFilter.end(); ++i)
	{
		tempHits[*i] = 0;
	}

	ReadsProcessor proc(kmerSize);
	//cut read into kmer size given
	while (rec.seq.length() >= (currentKmerNum + 1) * kmerSize) {

		const unsigned char* currentKmer = proc.prepSeq(rec.seq,
				currentKmerNum * kmerSize + startModifier1);

		//check to see if string is invalid
		if (currentKmer != NULL) {

			const unordered_map<string, bool> &results =
					filters[hashSig]->multiContains(currentKmer);

			//record hit number in order
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				if (results.at(*i)) {
					++tempHits[*i];
				}
			}
		}
		++currentKmerNum;
	}
	for (vector<string>::const_iterator i = idsInFilter.begin();
			i != idsInFilter.end(); ++i)
	{
		hits[*i] = tempHits.at(*i) >= minHit;
	}
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Updates hits value to number of hits (hashSig is used to as key)
 */
void BioBloomClassifier::evaluateReadStd(const FastqRecord &rec,
		const string &hashSig, unordered_map<string, bool> &hits)
{

	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*filters[hashSig]).getFilterIds();

	uint16_t kmerSize = infoFiles.at(hashSig).front()->getKmerSize();

	ReadsProcessor proc(kmerSize);

	double normalizationValue = rec.seq.length() - kmerSize + 1;
	double threshold = scoreThreshold * normalizationValue;

	for (vector<string>::const_iterator i = idsInFilter.begin();
			i != idsInFilter.end(); ++i)
	{
		bool pass = false;
		hits[*i] = false;
		if (minHit > 0) {
			uint16_t screeningHits = 0;
			size_t screeningLoc = rec.seq.length() % kmerSize / 2;
			//First pass filtering
			while (rec.seq.length() >= screeningLoc + kmerSize) {
				const unsigned char* currentKmer = proc.prepSeq(rec.seq,
						screeningLoc);
				if (currentKmer != NULL) {
					if (filtersSingle.at(*i)->contains(currentKmer)) {
						screeningHits++;
						if (screeningHits >= minHit) {
							pass = true;
							break;
						}
					}
				}
				screeningLoc += kmerSize;
			}
		} else {
			pass = true;
		}
		if (pass) {
			size_t currentLoc = 0;
			double score = 0;
			uint16_t streak = 0;
			while (rec.seq.length() >= currentLoc + kmerSize) {
				const unsigned char* currentKmer = proc.prepSeq(rec.seq,
						currentLoc);
				if (streak == 0) {
					if (currentKmer != NULL) {
						if (filtersSingle.at(*i)->contains(currentKmer)) {
							score += 0.5;
							++streak;
						}
						++currentLoc;
					} else {
						currentLoc += kmerSize + 1;
					}
				} else {
					if (currentKmer != NULL) {
						if (filtersSingle.at(*i)->contains(currentKmer)) {
							++streak;
							score += 1 - 1 / (2 * streak);
							++currentLoc;
							if (threshold <= score) {
								hits[*i] = true;
								break;
							}
							continue;
						}
					} else {
						currentLoc += kmerSize + 1;
					}
					if (streak < streakThreshold) {
						++currentLoc;
					} else {
						currentLoc += kmerSize;
					}
					streak = 0;
				}
			}
		}
	}
}

BioBloomClassifier::~BioBloomClassifier()
{
}

