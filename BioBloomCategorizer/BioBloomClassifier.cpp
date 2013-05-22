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

BioBloomClassifier::BioBloomClassifier(const vector<string> &filterFilePaths,
		size_t minHit, double percentMinHit, size_t maxHitValue,
		const string &prefix, const string &outputPostFix, uint8_t tileModifier) :
		minHit(minHit), percentMinHit(percentMinHit), filterNum(
				filterFilePaths.size()), maxHitValue(maxHitValue), noMatch(
				"noMatch"), multiMatch("multiMatch"), prefix(prefix), postfix(
				outputPostFix), tileModifier(tileModifier)
{
	loadFilters(filterFilePaths);
}

/*
 * Loads list of filters into memory
 * todo: Implement non-block I/O when loading multiple filters at once
 */
void BioBloomClassifier::filter(const vector<string> &inputFiles)
{
	ofstream readStatusOutput((prefix + "_status.tsv" + postfix).c_str(),
			ios::out);

	//variables for storing results summary
	unordered_map<string, size_t> aboveThreshold;
	unordered_map<string, size_t> belowThreshold;
	unordered_map<string, vector<size_t> > rawHits;

	size_t totalReads = 0;

	//print out header info and initialize variables
	readStatusOutput
			<< initSummaryVars(hashSigs, aboveThreshold, belowThreshold,
					rawHits);

	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it)
	{
		FastaReader sequence((*it).c_str(), FastaReader::NO_FOLD_CASE);
		FastqRecord rec;

		//stored out of loop so reallocation does not have to be done
		unordered_map<string, size_t> hits(filterNum);
		while (sequence >> rec) {

			//track read progress
			++totalReads;
			if (totalReads % 1000000 == 0) {
				cerr << "Currently Reading Read Number: " << totalReads << endl;
			}

			//print readID
			readStatusOutput << rec.id << "\t" << rec.seq.length();

			//initialize hits to zero
			initHits(hits);

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				evaluateRead(rec, *j, hits, tileModifier);
			}

			//print hit results to read status
			readStatusOutput << getReadStatStr(rec.id, rec.seq.length(), hits);

			//Evaluate hit data and record for summary
			updateSummaryData(rec.seq.length(), hits, aboveThreshold,
					belowThreshold, rawHits);
		}
	}

	readStatusOutput.flush();
	readStatusOutput.close();

	cerr << "Total Reads:" << totalReads << endl;
	printSummary(prefix, aboveThreshold, belowThreshold, totalReads);
	printCountSummary(prefix, rawHits, totalReads);
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * Prints reads into seperate files
 */
void BioBloomClassifier::filterPrint(const vector<string> &inputFiles)
{

	ofstream readStatusOutput((prefix + "_status.tsv" + postfix).c_str(),
			ios::out);

	bool printReads = 1;

	//variables for storing results summary
	unordered_map<string, size_t> aboveThreshold;
	unordered_map<string, size_t> belowThreshold;
	unordered_map<string, vector<size_t> > rawHits;

	size_t totalReads = 0;

	unordered_map<string, shared_ptr<ofstream> > outputFiles;
	shared_ptr<ofstream> no_match(
			new ofstream((prefix + "_" + noMatch + ".fastq" + postfix).c_str(),
					ios::out));
	shared_ptr<ofstream> multi_match(
			new ofstream(
					(prefix + "_" + multiMatch + ".fastq" + postfix).c_str(),
					ios::out));
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
			shared_ptr<ofstream> temp(
					new ofstream(
							(prefix + "_" + *i + ".fastq" + postfix).c_str(),
							ios::out));
			outputFiles[*i] = temp;
		}
	}

	//print out header info and initialize variables for summary
	readStatusOutput
			<< initSummaryVars(hashSigs, aboveThreshold, belowThreshold,
					rawHits);

	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it)
	{
		FastaReader sequence((*it).c_str(), FastaReader::NO_FOLD_CASE);
		FastqRecord rec;
		//hits results stored in hashmap of filternames and hits
		unordered_map<string, size_t> hits(filterNum);
		while (sequence >> rec) {
			++totalReads;
			if (totalReads % 1000000 == 0) {
				cerr << "Currently Reading Read Number: " << totalReads << endl;
			}

			//initialize hits to zero
			initHits(hits);

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				evaluateRead(rec, *j, hits, tileModifier);
			}

			//print hit results to read status
			readStatusOutput << getReadStatStr(rec.id, rec.seq.length(), hits);

			//Evaluate hit data and record for summary
			const string &outputFileName = updateSummaryData(rec.seq.length(),
					hits, aboveThreshold, belowThreshold, rawHits);

			(*outputFiles[outputFileName]) << "@" << rec.id << "\n" << rec.seq
					<< "\n+\n" << rec.qual << endl;
		}
	}

	readStatusOutput.flush();
	readStatusOutput.close();

	//close sorting files
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			outputFiles[*i]->flush();
			outputFiles[*i]->close();
		}
	}
	outputFiles[noMatch]->flush();
	outputFiles[noMatch]->close();
	outputFiles[multiMatch]->flush();
	outputFiles[multiMatch]->close();
	cerr << "Total Reads:" << totalReads << endl;
	printSummary(prefix, aboveThreshold, belowThreshold, totalReads);
	printCountSummary(prefix, rawHits, totalReads);
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPair(const string &file1, const string &file2)
{

	ofstream readStatusOutput((prefix + "_status.tsv" + postfix).c_str(),
			ios::out);

	//variables for storing results summary
	unordered_map<string, size_t> aboveThreshold;
	unordered_map<string, size_t> belowThreshold;
	unordered_map<string, vector<size_t> > rawHits;

	size_t totalReads = 0;

	//print out header info and initialize variables for summary
	readStatusOutput
			<< initSummaryVars(hashSigs, aboveThreshold, belowThreshold,
					rawHits);

	cerr << "Filtering Start" << "\n";

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
	FastqRecord rec1;
	FastqRecord rec2;
	//hits results stored in hashmap of filter names and hits
	unordered_map<string, size_t> hits1(filterNum);
	unordered_map<string, size_t> hits2(filterNum);

	while (sequence1 >> rec1 && sequence2 >> rec2) {
		++totalReads;
		if (totalReads % 1000000 == 0) {
			cerr << "Currently Reading Read Number: " << totalReads << endl;
		}

		//initialize hits to zero
		initHits(hits1);
		initHits(hits2);

		//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
		for (vector<string>::const_iterator j = hashSigs.begin();
				j != hashSigs.end(); ++j)
		{
			string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
			string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
			if (tempStr1 == tempStr2) {
				evaluateRead(rec1, *j, hits1, tileModifier);
				evaluateRead(rec2, *j, hits2, tileModifier);
			} else {
				cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
						<< tempStr2 << endl;
				exit(1);
			}
		}

		string readID = rec1.id.substr(0, rec1.id.length() - 2);

		//print hit results to read status
		readStatusOutput
				<< getReadStatStrPair(readID, rec1.seq.length(),
						rec2.seq.length(), hits1, hits2);

		//Evaluate hit data and record for summary
		updateSummaryData(rec1.seq.length(), hits1, aboveThreshold,
				belowThreshold, rawHits);
		updateSummaryData(rec1.seq.length(), hits1, aboveThreshold,
				belowThreshold, rawHits);

		for (vector<string>::const_iterator j = hashSigs.begin();
				j != hashSigs.end(); ++j)
		{
			//update summary
			const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				//print read status
				readStatusOutput << "\t" << hits1[*i] << '/' << hits2[*i];

				//pick threshold, by percent or by absolute value
				uint16_t kmerSize = (*(infoFiles[*j].front())).getKmerSize();
				size_t threshold1 = size_t(
						percentMinHit * (rec1.seq.length() / kmerSize));
				size_t threshold2 = size_t(
						percentMinHit * (rec2.seq.length() / kmerSize));
				if (minHit > threshold1) {
					threshold1 = minHit;
				}
				if (minHit > threshold2) {
					threshold2 = minHit;
				}

				if (hits1[*i] >= threshold1 && hits2[*i] >= threshold2) {
					++aboveThreshold[*i];
				} else if (hits1[*i] != 0 && hits2[*i] != 0) {
					++belowThreshold[*i];
				}
				//modify total reads
				if (rawHits[*i].size() > hits1[*i] + hits2[*i]) {
					rawHits[*i][hits1[*i] + hits2[*i]]++;
				}
			}
		}
		readStatusOutput << "\n";
	}
	if (sequence2 >> rec2 && sequence1.eof() && sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
	}

	readStatusOutput.flush();
	readStatusOutput.close();

	cerr << "Total Reads:" << totalReads << endl;
	printSummary(prefix, aboveThreshold, belowThreshold, totalReads);
	printCountSummary(prefix, rawHits, totalReads);
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * prints reads
 */
void BioBloomClassifier::filterPairPrint(const string &file1,
		const string &file2)
{
	ofstream readStatusOutput((prefix + "_status.tsv" + postfix).c_str(),
			ios::out);

	//variables for storing results summary
	unordered_map<string, size_t> aboveThreshold;
	unordered_map<string, size_t> belowThreshold;
	unordered_map<string, vector<size_t> > rawHits;

	size_t totalReads = 0;

	unordered_map<string, shared_ptr<ofstream> > outputFiles;
	shared_ptr<ofstream> noMatch1(
			new ofstream(
					(prefix + "_" + noMatch + "_1.fastq" + postfix).c_str(),
					ios::out));
	shared_ptr<ofstream> noMatch2(
			new ofstream(
					(prefix + "_" + noMatch + "_2.fastq" + postfix).c_str(),
					ios::out));
	shared_ptr<ofstream> multiMatch1(
			new ofstream(
					(prefix + "_" + multiMatch + "_1.fastq" + postfix).c_str(),
					ios::out));
	shared_ptr<ofstream> multiMatch2(
			new ofstream(
					(prefix + "_" + multiMatch + "_2.fastq" + postfix).c_str(),
					ios::out));
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
			shared_ptr<ofstream> temp1(
					new ofstream(
							(prefix + "_" + *i + "_1.fastq" + postfix).c_str(),
							ios::out));
			shared_ptr<ofstream> temp2(
					new ofstream(
							(prefix + "_" + *i + "_2.fastq" + postfix).c_str(),
							ios::out));
			outputFiles[*i + "1"] = temp1;
			outputFiles[*i + "2"] = temp2;
		}
	}

	//print out header info and initialize variables for summary
	readStatusOutput
			<< initSummaryVars(hashSigs, aboveThreshold, belowThreshold,
					rawHits);

	cerr << "Filtering Start" << "\n";

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
	FastqRecord rec1;
	FastqRecord rec2;
	//hits results stored in hashmap of filter names and hits
	unordered_map<string, size_t> hits1(filterNum);
	unordered_map<string, size_t> hits2(filterNum);

	while (sequence1 >> rec1 && sequence2 >> rec2) {
		++totalReads;
		if (totalReads % 1000000 == 0) {
			cerr << "Currently Reading Read Number: " << totalReads << endl;
		}

		//initialize hits to zero
		initHits(hits1);
		initHits(hits2);

		//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
		for (vector<string>::const_iterator j = hashSigs.begin();
				j != hashSigs.end(); ++j)
		{
			string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
			string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
			if (tempStr1 == tempStr2) {
				evaluateRead(rec1, *j, hits1, tileModifier);
				evaluateRead(rec2, *j, hits2, tileModifier);
			} else {
				cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
						<< tempStr2 << endl;
				exit(1);
			}
		}

		string readID = rec1.id.substr(0, rec1.id.length() - 2);

		//print hit results to read status
		readStatusOutput
				<< getReadStatStrPair(readID, rec1.seq.length(),
						rec1.seq.length(), hits1, hits2);

		//Evaluate hit data and record for summary
		const string &outputFileName = updateSummaryDataPair(rec1.seq.length(),
				rec2.seq.length(), hits1, hits2, aboveThreshold, belowThreshold,
				rawHits);

		(*outputFiles[outputFileName + "1"]) << "@" << rec1.id << "\n"
				<< rec1.seq << "\n+\n" << rec1.qual << endl;
		(*outputFiles[outputFileName + "2"]) << "@" << rec2.id << "\n"
				<< rec2.seq << "\n+\n" << rec2.qual << endl;
	}
	if (sequence2 >> rec2 && sequence1.eof() && sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
	}

	readStatusOutput.flush();
	readStatusOutput.close();

	//close sorting files
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			outputFiles[*i + "1"]->flush();
			outputFiles[*i + "2"]->close();
			outputFiles[*i + "1"]->flush();
			outputFiles[*i + "2"]->close();
		}
	}
	outputFiles[noMatch + "1"]->flush();
	outputFiles[noMatch + "1"]->close();
	outputFiles[noMatch + "2"]->flush();
	outputFiles[noMatch + "2"]->close();
	outputFiles[multiMatch + "1"]->flush();
	outputFiles[multiMatch + "1"]->close();
	outputFiles[multiMatch + "2"]->flush();
	outputFiles[multiMatch + "2"]->close();
	cerr << "Total Reads:" << totalReads << endl;
	printSummary(prefix, aboveThreshold, belowThreshold, totalReads);
	printCountSummary(prefix, rawHits, totalReads);
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * Prints reads into seperate files
 */
void BioBloomClassifier::filterPairBAMPrint(const string &file)
{
	ofstream readStatusOutput((prefix + "_status.tsv" + postfix).c_str(),
			ios::out);

	//variables for storing results summary
	unordered_map<string, size_t> aboveThreshold;
	unordered_map<string, size_t> belowThreshold;
	unordered_map<string, vector<size_t> > rawHits;

	unordered_map<string, FastqRecord> unPairedReads;

	size_t totalReads = 0;

	unordered_map<string, shared_ptr<ofstream> > outputFiles;
	shared_ptr<ofstream> noMatch1(
			new ofstream(
					(prefix + "_" + noMatch + "_1.fastq" + postfix).c_str(),
					ios::out));
	shared_ptr<ofstream> noMatch2(
			new ofstream(
					(prefix + "_" + noMatch + "_2.fastq" + postfix).c_str(),
					ios::out));
	shared_ptr<ofstream> multiMatch1(
			new ofstream(
					(prefix + "_" + multiMatch + "_1.fastq" + postfix).c_str(),
					ios::out));
	shared_ptr<ofstream> multiMatch2(
			new ofstream(
					(prefix + "_" + multiMatch + "_2.fastq" + postfix).c_str(),
					ios::out));
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
			shared_ptr<ofstream> temp1(
					new ofstream(
							(prefix + "_" + *i + "_1.fastq" + postfix).c_str(),
							ios::out));
			shared_ptr<ofstream> temp2(
					new ofstream(
							(prefix + "_" + *i + "_2.fastq" + postfix).c_str(),
							ios::out));
			outputFiles[*i + "1"] = temp1;
			outputFiles[*i + "2"] = temp2;
		}
	}

	//print out header info and initialize variables for summary
	readStatusOutput
			<< initSummaryVars(hashSigs, aboveThreshold, belowThreshold,
					rawHits);

	cerr << "Filtering Start" << "\n";

	FastaReader sequence(file.c_str(), FastaReader::NO_FOLD_CASE);
	//hits results stored in hashmap of filter names and hits
	unordered_map<string, size_t> hits1(filterNum);
	unordered_map<string, size_t> hits2(filterNum);

	while (!sequence.eof()) {
		FastqRecord rec;
		if (sequence >> rec) {
			string readID = rec.id.substr(0, rec.id.length() - 2);
			if (unPairedReads.find(readID) != unPairedReads.end()) {

				const FastqRecord &rec1 =
						rec.id.at(rec.id.length() - 1) == '1' ?
								rec : unPairedReads[readID];
				const FastqRecord &rec2 =
						rec.id.at(rec.id.length() - 1) == '2' ?
								rec : unPairedReads[readID];

				++totalReads;
				if (totalReads % 1000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}

				//initialize hits to zero
				initHits(hits1);
				initHits(hits2);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					string tempStr1 = rec1.id.substr(0,
							rec1.id.find_last_of("/"));
					string tempStr2 = rec2.id.substr(0,
							rec2.id.find_last_of("/"));
					if (tempStr1 == tempStr2) {
						evaluateRead(rec1, *j, hits1, tileModifier);
						evaluateRead(rec2, *j, hits2, tileModifier);
					} else {
						cerr << "Read IDs do not match" << "\n" << tempStr1
								<< "\n" << tempStr2 << endl;
						exit(1);
					}
				}

				//print hit results to read status
				readStatusOutput
						<< getReadStatStrPair(readID, rec1.seq.length(),
								rec2.seq.length(), hits1, hits2);

				//Evaluate hit data and record for summary
				const string &outputFileName = updateSummaryDataPair(
						rec1.seq.length(), rec2.seq.length(), hits1, hits2,
						aboveThreshold, belowThreshold, rawHits);

				(*outputFiles[outputFileName + "1"]) << "@" << rec1.id << "\n"
						<< rec1.seq << "\n+\n" << rec1.qual << endl;
				(*outputFiles[outputFileName + "2"]) << "@" << rec2.id << "\n"
						<< rec2.seq << "\n+\n" << rec2.qual << endl;

				//clean up reads
				unPairedReads.erase(readID);

			} else {
				unPairedReads[readID] = rec;
			}
		}
	}

	readStatusOutput.flush();

	//close sorting files
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			outputFiles[*i + "1"]->flush();
			outputFiles[*i + "2"]->close();
			outputFiles[*i + "1"]->flush();
			outputFiles[*i + "2"]->close();
		}
	}
	outputFiles[noMatch + "1"]->flush();
	outputFiles[noMatch + "1"]->close();
	outputFiles[noMatch + "2"]->flush();
	outputFiles[noMatch + "2"]->close();
	outputFiles[multiMatch + "1"]->flush();
	outputFiles[multiMatch + "1"]->close();
	outputFiles[multiMatch + "2"]->flush();
	outputFiles[multiMatch + "2"]->close();
	cerr << "Total Reads:" << totalReads << endl;
	printSummary(prefix, aboveThreshold, belowThreshold, totalReads);
	printCountSummary(prefix, rawHits, totalReads);
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPairBAM(const string &file)
{
	ofstream readStatusOutput((prefix + "_status.tsv" + postfix).c_str(),
			ios::out);

	//variables for storing results summary
	unordered_map<string, size_t> aboveThreshold;
	unordered_map<string, size_t> belowThreshold;
	unordered_map<string, vector<size_t> > rawHits;

	unordered_map<string, FastqRecord> unPairedReads;

	size_t totalReads = 0;

	//print out header info and initialize variables for summary
	readStatusOutput
			<< initSummaryVars(hashSigs, aboveThreshold, belowThreshold,
					rawHits);

	cerr << "Filtering Start" << "\n";

	FastaReader sequence(file.c_str(), FastaReader::NO_FOLD_CASE);
	//hits results stored in hashmap of filter names and hits
	unordered_map<string, size_t> hits1(filterNum);
	unordered_map<string, size_t> hits2(filterNum);

	while (!sequence.eof()) {
		FastqRecord rec;
		if (sequence >> rec) {
			string readID = rec.id.substr(0, rec.id.length() - 2);
			if (unPairedReads.find(readID) != unPairedReads.end()) {

				const FastqRecord &rec1 =
						rec.id.at(rec.id.length() - 1) == '1' ?
								rec : unPairedReads[readID];
				const FastqRecord &rec2 =
						rec.id.at(rec.id.length() - 1) == '2' ?
								rec : unPairedReads[readID];

				++totalReads;
				if (totalReads % 1000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}

				//initialize hits to zero
				initHits(hits1);
				initHits(hits2);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					string tempStr1 = rec1.id.substr(0,
							rec1.id.find_last_of("/"));
					string tempStr2 = rec2.id.substr(0,
							rec2.id.find_last_of("/"));
					if (tempStr1 == tempStr2) {
						evaluateRead(rec1, *j, hits1, tileModifier);
						evaluateRead(rec2, *j, hits2, tileModifier);
					} else {
						cerr << "Read IDs do not match" << "\n" << tempStr1
								<< "\n" << tempStr2 << endl;
						exit(1);
					}
				}

				//print hit results to read status
				readStatusOutput
						<< getReadStatStrPair(readID, rec1.seq.length(),
								rec2.seq.length(), hits1, hits2);

				//Evaluate hit data and record for summary
				updateSummaryDataPair(rec1.seq.length(), rec2.seq.length(),
						hits1, hits2, aboveThreshold, belowThreshold, rawHits);

				//clean up reads
				unPairedReads.erase(readID);

			} else {
				unPairedReads[readID] = rec;
			}
		}
	}

	readStatusOutput.flush();
	readStatusOutput.close();

	cerr << "Total Reads:" << totalReads << endl;
	printSummary(prefix, aboveThreshold, belowThreshold, totalReads);
	printCountSummary(prefix, rawHits, totalReads);
}

/*
 * Loads list of filters into memory
 * todo: Implement non-block I/O when loading multiple filters at once
 */
void BioBloomClassifier::loadFilters(const vector<string> &filterFilePaths)
{
	cout << "Starting to Load Filters." << endl;
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
		hashSig << info->getKmerSize() << info->getSeedHashSigniture();

		//if hashSig exists add filter to list
		if (infoFiles.count(hashSig.str()) == 1) {
			infoFiles[hashSig.str()].push_back(info);
			filters[hashSig.str()]->addFilter(info->getCalcuatedFilterSize(),
					info->getFilterID(), *it);
		} else {
			vector<shared_ptr<BloomFilterInfo> > tempVect;
			tempVect.push_back(info);
			vector<size_t>::const_iterator seeds =
					info->getSeedValues().begin();
			//Create HashManager for MultiFilter
			HashManager hashMan;
			for (vector<string>::const_iterator hashFn =
					info->getHashFunctionNames().begin();
					hashFn != info->getHashFunctionNames().end(); ++hashFn)
			{
				hashMan.addHashFunction(*hashFn, *seeds);
				++seeds;
			}
			hashSigs.push_back(hashSig.str());
			shared_ptr<MultiFilter> temp(new MultiFilter(hashMan));
			filters[hashSig.str()] = temp;
			filters[hashSig.str()]->addFilter(info->getCalcuatedFilterSize(),
					info->getFilterID(), *it);
			infoFiles[hashSig.str()] = tempVect;
		}
		cerr << "Loaded Filter: " + info->getFilterID() << endl;
	}
	cerr << "Filter Loading Complete." << endl;
}

/*
 * Prints summary information:
 * -total reads over/under threshold
 * -percent reads over threshold
 * -total reads that don't hit filter at all
 */
void BioBloomClassifier::printSummary(const string &prefix,
		unordered_map<string, size_t> &aboveThreshold,
		unordered_map<string, size_t> &belowThreshold, size_t totalReads)
{
	ofstream summaryOutput((prefix + "_summary.tsv" + postfix).c_str(),
			ios::out);
	summaryOutput << "type";
	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> &idsInFilter = filters[*j]->getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << *i << "_"
					<< (*(infoFiles[*j].front())).getKmerSize();
		}
	}
	summaryOutput << "\n";

	//print summary information and close filehandles
	summaryOutput << "\nHits";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(aboveThreshold[*i]) / double(totalReads);
		}
	}
	summaryOutput << "\nMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(belowThreshold[*i]) / double(totalReads);
		}
	}
	summaryOutput << "\nConfidentMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(
							totalReads - belowThreshold[*i]
									- aboveThreshold[*i]) / double(totalReads);
		}
	}

	summaryOutput << "\n\nHits";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << aboveThreshold[*i];
		}
	}
	summaryOutput << "\nMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << belowThreshold[*i];
		}
	}
	summaryOutput << "\nConfidentMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< totalReads - belowThreshold[*i] - aboveThreshold[*i];
		}
	}
	summaryOutput.close();
}

/*
 * Prints raw read counts file
 */
void BioBloomClassifier::printCountSummary(const string &prefix,
		unordered_map<string, vector<size_t> > &rawHits, size_t total)
{
	if (maxHitValue > 0) {
		ofstream summaryOutput((prefix + "_rawCounts.tsv" + postfix).c_str(),
				ios::out);
		summaryOutput << "type";
		//initialize variables and print filter ids
		for (vector<string>::const_iterator j = hashSigs.begin();
				j != hashSigs.end(); ++j)
		{
			const vector<string> &idsInFilter = filters[*j]->getFilterIds();
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				summaryOutput << "\t" << *i << "_"
						<< (*(infoFiles[*j].front())).getKmerSize();
			}
		}
		summaryOutput << "\n";

		uint16_t currentHitVal = 0;
		size_t runningTotal = 0;

		while (currentHitVal < maxHitValue) {
			summaryOutput << currentHitVal;
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				const vector<string> &idsInFilter = filters[*j]->getFilterIds();
				for (vector<string>::const_iterator i = idsInFilter.begin();
						i != idsInFilter.end(); ++i)
				{
					if (rawHits[*i].size() < currentHitVal) {
						summaryOutput << "\t0";
					} else {
						runningTotal += rawHits[*i][currentHitVal];
						summaryOutput << "\t" << rawHits[*i][currentHitVal];
					}
				}
			}
			summaryOutput << "\n";
			currentHitVal++;
		}
		summaryOutput << ">" << maxHitValue << "\t" << total - runningTotal
				<< "\n";
		summaryOutput.flush();
		summaryOutput.close();
	}
}

//helper methods

/*
 * checks if file exists
 */
bool BioBloomClassifier::fexists(const string &filename) const
{
	ifstream ifile(filename.c_str());
	return ifile;
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Updates hits value to number of hits (hashSig is used to as key)
 */
//@Todo: Implement in more efficient way (if read has a miss no
void BioBloomClassifier::evaluateRead(const FastqRecord &rec,
		const string &hashSig, unordered_map<string, size_t> &hits,
		uint8_t tileModifier)
{
	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*filters[hashSig]).getFilterIds();

	//get kmersize for set of info files
	uint16_t kmerSize = (*(infoFiles[hashSig].front())).getKmerSize();

	//Establish tiling pattern
	uint16_t startModifier = (rec.seq.length() % (kmerSize + tileModifier))
			/ 2;

	ReadsProcessor proc(kmerSize);

	uint8_t currentLoc = 0;
	//cut read into kmer size + tilemodifier given
	while (rec.seq.length() >= currentLoc + kmerSize + tileModifier)
	{
		unordered_map<string, bool> tempResults;
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			tempResults[*i] = true;
		}

		for (uint8_t j = 0; j <= tileModifier; ++j) {
			const string &currentKmer = proc.prepSeq(rec.seq, currentLoc + startModifier + j);

			//check to see if string is invalid
			if (!currentKmer.empty()) {
				const unordered_map<string, bool> &results =
						filters[hashSig]->multiContains(currentKmer);

				for (vector<string>::const_iterator i = idsInFilter.begin();
						i != idsInFilter.end(); ++i)
				{
					if (!results.find(*i)->second) {
						tempResults[*i] = false;
					}
				}
			} else {
				for (vector<string>::const_iterator i = idsInFilter.begin();
						i != idsInFilter.end(); ++i)
				{
					tempResults[*i] = false;
				}
			}
		}
		//record hit number in order
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			if (tempResults.find(*i)->second) {
				++hits[*i];
			}
		}
		currentLoc += kmerSize + tileModifier;
	}
}

///*
// * For a single read evaluate hits for a single hash signature
// * Sections with ambiguity bases are treated as misses
// * Updates hits value to number of hits (hashSig is used to as key)
// */
//void BioBloomClassifier::evaluateRead(const FastqRecord &rec,
//		const string &hashSig, unordered_map<string, size_t> &hits)
//{
//	//get filterIDs to iterate through has in a consistent order
//	const vector<string> &idsInFilter = (*filters[hashSig]).getFilterIds();
//
//	//get kmersize for set of info files
//	uint16_t kmerSize = (*(infoFiles[hashSig].front())).getKmerSize();
//
//	//Establish tiling pattern
//	uint16_t startModifier1 = (rec.seq.length() % kmerSize) / 2;
//	size_t currentKmerNum = 0;
//
//	ReadsProcessor proc(kmerSize);
//	//cut read into kmer size given
//	while (rec.seq.length() >= (currentKmerNum + 1) * kmerSize) {
//
//		const string &currentKmer = proc.prepSeq(rec.seq,
//				currentKmerNum * kmerSize + startModifier1);
//
//		//check to see if string is invalid
//		if (!currentKmer.empty()) {
//			const unordered_map<string, bool> &results =
//					filters[hashSig]->multiContains(currentKmer);
//
//			//record hit number in order
//			for (vector<string>::const_iterator i = idsInFilter.begin();
//					i != idsInFilter.end(); ++i)
//			{
//				if (results.find(*i)->second) {
//					++hits[*i];
//				}
//			}
//		}
//		++currentKmerNum;
//	}
//}

/*
 * Initializes Summary Variables. Also prints heads for read status.
 */
//Todo: move printing of read status heads to other function?
const string BioBloomClassifier::initSummaryVars(vector<string> &hashSigs,
		unordered_map<string, size_t> &aboveThreshold,
		unordered_map<string, size_t> &belowThreshold,
		unordered_map<string, vector<size_t> > &rawHits)
{
	stringstream readStatusOutput;
	readStatusOutput << "readID\tseqSize";

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			readStatusOutput << "\t" << *i << "_"
					<< (*(infoFiles[*j].front())).getKmerSize();
			aboveThreshold[*i] = 0;
			belowThreshold[*i] = 0;

			//initialize rawHits
			vector<size_t> temp(maxHitValue);
			uint16_t counter = 0;
			rawHits[*i] = temp;
			while (counter < maxHitValue) {
				rawHits[*i][counter] = 0;
				counter++;
			}
		}
	}
	readStatusOutput << "\n";
	return readStatusOutput.str();
}

/*
 * Initializes hits results to zero
 */
void BioBloomClassifier::initHits(unordered_map<string, size_t> &hits)
{
	//initialize hits to zero
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			hits[*i] = 0;
		}
	}
}

/*
 * return results of hits to be output for read status
 */
const string BioBloomClassifier::getReadStatStr(string const &readID,
		size_t readLength, unordered_map<string, size_t> &hits)
{
	stringstream str;
	str << readID << "\t" << readLength;
	//print readID
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		//update summary
		const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			//print to file
			str << str << "\t" << hits[*i];
		}
	}
	str << "\n";
	return str.str();
}

/*
 * return results of hits to be output for read status for paired end mode
 */
const string BioBloomClassifier::getReadStatStrPair(string const &readID,
		size_t readLength1, size_t readLength2,
		unordered_map<string, size_t> &hits1,
		unordered_map<string, size_t> &hits2)
{
	stringstream str;
	str << readID << "\t" << readLength1 << "|" << readLength2;
	//print readID
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		//update summary
		const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			//print to file
			str << "\t" << hits1[*i] << "|" << hits2[*i];
		}
	}
	str << "\n";
	return str.str();
}

/*
 * Records data for read status file based on thresholds
 * Returns qualifying read IDs that meet threshold
 */
const string BioBloomClassifier::updateSummaryData(size_t seqLen,
		unordered_map<string, size_t> &hits,
		unordered_map<string, size_t> &aboveThreshold,
		unordered_map<string, size_t> &belowThreshold,
		unordered_map<string, vector<size_t> > &rawHits)
{
	string filterID = noMatch;

	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		//update summary
		const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			//pick threshold, by percent or by absolute value
			uint16_t kmerSize = (*(infoFiles[*j].front())).getKmerSize();
			size_t threshold = size_t(percentMinHit * (seqLen / kmerSize));
			if (minHit > threshold) {
				threshold = minHit;
			}

			if (hits[*i] >= threshold) {
				++aboveThreshold[*i];
				if (filterID == noMatch) {
					filterID = *i;
				} else {
					filterID = multiMatch;
				}
			} else if (hits[*i] != 0) {
				++belowThreshold[*i];
			}

			//modify total reads
			if (rawHits[*i].size() > hits[*i]) {
				rawHits[*i][hits[*i]]++;
			}
		}
	}
	return filterID;
}

/*
 * Records data for read status file based on thresholds
 * Returns qualifying read IDs that meet threshold
 * both reads must qualify
 */
const string BioBloomClassifier::updateSummaryDataPair(size_t seqLen1,
		size_t seqLen2, unordered_map<string, size_t> &hits1,
		unordered_map<string, size_t> &hits2,
		unordered_map<string, size_t> &aboveThreshold,
		unordered_map<string, size_t> &belowThreshold,
		unordered_map<string, vector<size_t> > &rawHits)
{
	string filterID = noMatch;

	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		//update summary
		const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			//pick threshold, by percent or by absolute value
			uint16_t kmerSize = (*(infoFiles[*j].front())).getKmerSize();
			size_t threshold1 = size_t(percentMinHit * (seqLen1 / kmerSize));
			size_t threshold2 = size_t(percentMinHit * (seqLen2 / kmerSize));
			if (minHit > threshold1) {
				threshold1 = minHit;
			}
			if (minHit > threshold2) {
				threshold2 = minHit;
			}

			if (hits1[*i] >= threshold1 && hits2[*i] >= threshold2) {
				++aboveThreshold[*i];
				if (filterID == noMatch) {
					filterID = *i;
				} else {
					filterID = multiMatch;
				}
			} else if (hits1[*i] != 0 && hits2[*i] != 0) {
				++belowThreshold[*i];
			}
			//modify total reads
			if (rawHits[*i].size() > hits1[*i] + hits2[*i]) {
				rawHits[*i][hits1[*i] + hits2[*i]]++;
			}
		}
	}
	return filterID;
}

BioBloomClassifier::~BioBloomClassifier()
{
}

