/*
 * BloomMapClassifier.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: cjustin
 */

#include <BloomMapClassifier.h>
#include "DataLayer/FastaReader.h"
#include "Common/Options.h"
#include "Options.h"
#include <iostream>
#include "Common/Dynamicofstream.h"
#include "ResultsManager.h"

BloomMapClassifier::BloomMapClassifier(const string &filterFile) :
		m_filter(BloomMapSSBitVec<ID>(filterFile)) {
	cerr << "FPR given allowed misses: "
			<< m_filter.getFPR(
					m_filter.getSeedValues().size() - opt::allowMisses)
			<< endl;
	//load in ID file
	string idFile = (filterFile).substr(0, (filterFile).length() - 3)
			+ "_ids.txt";
	m_ids.set_empty_key(opt::EMPTY);
	if (!fexists(idFile)) {
		cerr
				<< "Error: " + (idFile)
						+ " File cannot be opened. A corresponding id file is needed."
				<< endl;
		exit(1);
	}
	else{
		cerr << "loading ID file: " << idFile << endl;
	}

	ifstream idFH(idFile.c_str(), ifstream::in);
	string line;
	getline(idFH, line);
	while (idFH.good()) {
		ID id;
		string fullID;
		stringstream converter(line);
		converter >> id;
		converter >> fullID;
		m_ids[id] = fullID;
		m_fullIDs.push_back(fullID);
		getline(idFH, line);
	}
	idFH.close();
}

void BloomMapClassifier::filter(const vector<string> &inputFiles) {
	size_t totalReads = 0;

	assert(opt::outputType != "");

	google::dense_hash_map<string, boost::shared_ptr<Dynamicofstream> > outputFiles;
	outputFiles.set_empty_key("");

	boost::shared_ptr<Dynamicofstream> no_match(
			new Dynamicofstream(
					opt::outputPrefix + "_" + NO_MATCH + "." + opt::outputType
							+ opt::filePostfix));
	boost::shared_ptr<Dynamicofstream> multi_match(
			new Dynamicofstream(
					opt::outputPrefix + "_" + MULTI_MATCH + "."
							+ opt::outputType + opt::filePostfix));
	outputFiles[NO_MATCH] = no_match;
	outputFiles[MULTI_MATCH] = multi_match;

	//print out header info and initialize variables
	ResultsManager resSummary(m_fullIDs, false);

	//initialize variables
	for (vector<string>::const_iterator i = m_fullIDs.begin();
			i != m_fullIDs.end(); ++i) {
		boost::shared_ptr<Dynamicofstream> temp(
				new Dynamicofstream(
						opt::outputPrefix + "_" + *i + "." + opt::outputType
								+ opt::filePostfix));
		outputFiles[*i] = temp;
	}

	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it) {
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
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				google::dense_hash_map<ID, unsigned> hitCounts;
				hitCounts.set_empty_key(opt::EMPTY);
				unsigned score = evaluateRead(rec, hitCounts);
				unsigned threshold = opt::score
						* (rec.seq.size() - m_filter.getKmerSize() + 1);
				vector<ID> hits;
				if(score > threshold){
					convertToHits(hitCounts, hits);
				}
				const string &outputFileName = resSummary.updateSummaryData(hits);
				printSingleToFile(outputFileName, rec, outputFiles);
			} else
				break;
		}
	}

	//close sorting files
	for (google::dense_hash_map<string, boost::shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j)
	{
		j->second->close();
		cerr << "File written to: "
				<< opt::outputPrefix + "_" + j->first + "." + opt::outputType + opt::filePostfix
				<< endl;
	}
	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << opt::outputPrefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

void BloomMapClassifier::filterPair(const string &file1, const string &file2) {
	size_t totalReads = 0;

	assert(opt::outputType != "");

	google::dense_hash_map<string, boost::shared_ptr<Dynamicofstream> > outputFiles;
	outputFiles.set_empty_key("");
	boost::shared_ptr<Dynamicofstream> noMatch1(
			new Dynamicofstream(
					opt::outputPrefix + "_" + NO_MATCH + "_1." + opt::outputType
							+ opt::filePostfix));
	boost::shared_ptr<Dynamicofstream> noMatch2(
			new Dynamicofstream(
					opt::outputPrefix + "_" + NO_MATCH + "_2." + opt::outputType
							+ opt::filePostfix));
	boost::shared_ptr<Dynamicofstream> multiMatch1(
			new Dynamicofstream(
					opt::outputPrefix + "_" + MULTI_MATCH + "_1." + opt::outputType
							+ opt::filePostfix));
	boost::shared_ptr<Dynamicofstream> multiMatch2(
			new Dynamicofstream(
					opt::outputPrefix + "_" + MULTI_MATCH + "_2." + opt::outputType
							+ opt::filePostfix));
	outputFiles[NO_MATCH + "_1"] = noMatch1;
	outputFiles[NO_MATCH + "_2"] = noMatch2;
	outputFiles[MULTI_MATCH + "_1"] = multiMatch1;
	outputFiles[MULTI_MATCH + "_2"] = multiMatch2;

	//print out header info and initialize variables
	ResultsManager resSummary(m_fullIDs, false);

	//initialize variables
	for (vector<string>::const_iterator i = m_fullIDs.begin();
			i != m_fullIDs.end(); ++i) {
		boost::shared_ptr<Dynamicofstream> temp1(
				new Dynamicofstream(
						opt::outputPrefix + "_" + *i + "_1." + opt::outputType
								+ opt::filePostfix));
		boost::shared_ptr<Dynamicofstream> temp2(
				new Dynamicofstream(
						opt::outputPrefix + "_" + *i + "_2." + opt::outputType
								+ opt::filePostfix));
		outputFiles[*i + "_1"] = temp1;
		outputFiles[*i + "_2"] = temp2;
	}

	cerr << "Filtering Start" << endl;

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
				if (totalReads % 10000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}
			}
			google::dense_hash_map<ID, unsigned> hitCounts1;
			google::dense_hash_map<ID, unsigned> hitCounts2;
			hitCounts1.set_empty_key(opt::EMPTY);
			hitCounts2.set_empty_key(opt::EMPTY);
			unsigned score1 = evaluateRead(rec1, hitCounts1);
			unsigned score2 = evaluateRead(rec1, hitCounts1);
			unsigned threshold1 = opt::score
					* (rec1.seq.size() - m_filter.getKmerSize() + 1);
			//TODO currently assuming both reads are the same length - May need change
//				unsigned threshold2 = opt::score
//						* (rec2.seq.size() - m_filter.getKmerSize() + 1);

			vector<ID> hits;
			//TODO default is currently inclusive mode(BOTH)
			if (score1 > threshold1 && score2 > threshold1) {
				convertToHitsBoth(hitCounts1, hitCounts2, hits);
			}
			const string &outputFileName = resSummary.updateSummaryData(hits);
			printPairToFile(outputFileName, rec1, rec2, outputFiles);
		} else
			break;
	}
	if (!sequence1.eof() || !sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
	}


	//close sorting files
	for (google::dense_hash_map<string, boost::shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j)
	{
		j->second->close();
		cerr << "File written to: "
				<< opt::outputPrefix + "_" + j->first + "." + opt::outputType + opt::filePostfix
				<< endl;
	}
	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << opt::outputPrefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

BloomMapClassifier::~BloomMapClassifier() {
}

