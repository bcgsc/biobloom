/*
 * BloomMapClassifier.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: cjustin
 */

#include <BloomMapClassifier.h>
#include "DataLayer/FastaReader.h"
#include "bloomfilter/RollingHashIterator.h"
#include "Common/Options.h"
#include "Options.h"
#include <iostream>
#include "Common/Dynamicofstream.h"
#include "ResultsManager.h"

BloomMapClassifier::BloomMapClassifier(const string &filterFile) :
		m_filter(BloomMap<ID>(filterFile)) {
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

					//for current read print out ids
					RollingHashIterator itr(rec.seq, m_filter.getHashNum(),
							m_filter.getKmerSize());
					google::dense_hash_map<ID, unsigned> tmpHash;
					tmpHash.set_empty_key(0);
					unsigned maxCount = 0;
					ID value = 0;
					unsigned nonZeroCount = 0;
					unsigned totalCount = 0;

					vector<ID> hits;

					while (itr != itr.end()) {
						ID id = m_filter.atBest(*itr);

						if (id != 0) {
							if (tmpHash.find(id) != tmpHash.end()) {
								++tmpHash[id];
							}
							if (maxCount == tmpHash[id]) {
								value = numeric_limits<ID>::max();
							} else if (maxCount < tmpHash[id]) {
								value = id;
								maxCount = tmpHash[id];
							}
							++nonZeroCount;
						}
						++totalCount;
						itr++;
					}
					double score = double(nonZeroCount) / double(totalCount);

					if (opt::score < score && value != opt::COLLI) {
						//currently only one hit (no multimatches)
						hits.push_back(value);
					}

					const string &outputFileName = resSummary.updateSummaryData(hits);
					printSingleToFile(outputFileName, rec, outputFiles);
				}
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

BloomMapClassifier::~BloomMapClassifier() {
}
