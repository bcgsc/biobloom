/*
 * BloomMapClassifier.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: cjustin
 */

#include <BloomMapClassifier.h>
#include <ResultsManager.hpp>
#include "Common/Options.h"
#include "Options.h"
#include <iostream>
#include <zlib.h>
#include "Common/kseq.h"
//KSEQ_INIT(gzFile, gzread)

BloomMapClassifier::BloomMapClassifier(const string &filterFile) :
		m_filter(MIBloomFilter<ID>(filterFile)) {
	cerr << "FPR given allowed misses: "
			<< m_filter.getFPR(opt::allowMisses)
			<< endl;
	//load in ID file
	string idFile = (filterFile).substr(0, (filterFile).length() - 3)
			+ "_ids.txt";
	google::dense_hash_map<ID, string> ids;
	google::dense_hash_map<ID, boost::shared_ptr< vector<ID> > > colliIDs;
	colliIDs.set_empty_key(opt::EMPTY);
	ids.set_empty_key(opt::EMPTY);
	if (!fexists(idFile)) {
		cerr << "Error: " << idFile.c_str()
				<< " File cannot be opened. A corresponding id file is needed."
				<< endl;
		exit(1);
	}
	else{
		cerr << "Loading ID file: " << idFile.c_str() << endl;
	}

	ifstream idFH(idFile.c_str(), ifstream::in);
	string line;
	getline(idFH, line);
	while (idFH.good()) {
		ID id;
		stringstream converter(line);
		converter >> id;

		colliIDs[id] = boost::shared_ptr<vector<ID> >(new vector<ID>());
		string fullID;
		unsigned idCount = 0;
		for (;converter >> fullID; ++idCount) {
			ID tempID = atoi(fullID.c_str());
			colliIDs[id]->push_back(tempID);
		}
		if (idCount == 1) {
			colliIDs[id]->clear();
			colliIDs[id]->push_back(id);
			ids[id] = fullID;
		}
		getline(idFH, line);
	}
	m_fullIDs = vector<string>(ids.size() + 1); //first element will be empty
	m_fullIDs[0] = "missing"; //set first element to missing
	for (google::dense_hash_map<ID, string>::const_iterator itr = ids.begin();
			itr != ids.end(); ++itr) {
		m_fullIDs[itr->first] = itr->second;
	}
	m_colliIDs = vector< boost::shared_ptr<vector<ID> > >(colliIDs.size() + 1);
	for (google::dense_hash_map<ID, boost::shared_ptr< vector<ID> > >::const_iterator itr =
			colliIDs.begin(); itr != colliIDs.end(); ++itr) {
		if(itr->second->size() <= opt::maxGroupSize)
			m_colliIDs[itr->first] = itr->second;
		else
			m_colliIDs[itr->first] = NULL;
	}
	idFH.close();
}

void BloomMapClassifier::filter(const vector<string> &inputFiles) {
	size_t totalReads = 0;

	string outputName = opt::outputPrefix + "_reads.tsv";

	//TODO output fasta?
	if(opt::outputType == "fq"){
		outputName = opt::outputPrefix + "_reads.fq";
	}

	Dynamicofstream readsOutput(outputName);

	//print out header info and initialize variables
	ResultsManager<ID> resSummary(m_fullIDs, false);

	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it) {
		gzFile fp;
		fp = gzopen(it->c_str(), "r");
		if (fp == NULL) {
			cerr << "Cannot open file" << it->c_str() << endl;
			exit(1);
		}
		kseq_t *seq = kseq_init(fp);
		int l;
#pragma omp parallel private(l)
		for (string sequence;;) {
			string name;
			string qual;
#pragma omp critical(sequence)
			{
				l = kseq_read(seq);
				sequence = string(seq->seq.s, seq->seq.l);
				name = string(seq->name.s, seq->name.l);
				qual = string(seq->qual.s, seq->qual.l);
			}
			if (l >= 0) {
#pragma omp atomic
				++totalReads;
#pragma omp critical(totalReads)
				if (totalReads % 1000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}

				google::dense_hash_map<ID, unsigned> hitCounts;
				hitCounts.set_empty_key(opt::EMPTY);
				unsigned score = evaluateRead(sequence, hitCounts);
				unsigned threshold = opt::score
						* ((l - m_filter.getKmerSize() + 1) * (opt::allowMisses + 1));
				vector<ID> hits;
				bool aboveThreshold = score > threshold;
				if (aboveThreshold) {
					convertToHits(hitCounts, hits);
				}
				ID idIndex = resSummary.updateSummaryData(hits);
				if (idIndex != opt::EMPTY) {
					const string &fullID =
							idIndex == opt::COLLI ?
									UNKNOWN : m_fullIDs.at(idIndex);
					if (opt::outputType == "fq") {
#pragma omp critical(outputFiles)
						{
							readsOutput << "@" << name << " " << fullID << "\n"
									<< sequence << "\n+\n" << qual << "\n";
						}
					}else {
#pragma omp critical(outputFiles)
						{
							readsOutput << fullID << "\t" << name << "\t"
									<< score;
							for (google::dense_hash_map<ID, unsigned>::iterator itr =
									hitCounts.begin(); itr != hitCounts.end();
									++itr) {
								readsOutput << "\t" << m_fullIDs.at(itr->first)
										<< ":" << itr->second;
							}
							readsOutput << "\n";
						}
					}
				}
			} else {
				break;
			}
		}
		kseq_destroy(seq);
		gzclose(fp);
	}

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

void BloomMapClassifier::filterPair(const string &file1, const string &file2) {
	size_t totalReads = 0;
	string outputName = opt::outputPrefix + "_reads.tsv";

	//TODO output fasta?
	if(opt::outputType == "fq"){
		outputName = opt::outputPrefix + "_reads.fq";
	}

	Dynamicofstream readsOutput(outputName);

	//print out header info and initialize variables
	ResultsManager<ID> resSummary(m_fullIDs, false);

	cerr << "Filtering Start" << endl;

	gzFile fp1;
	gzFile fp2;
	fp1 = gzopen(file1.c_str(), "r");
	if (fp1 == NULL) {
		cerr << "Cannot open file " << file1.c_str() << endl;
		exit(1);
	}
	fp2 = gzopen(file2.c_str(), "r");
	if (fp2 == NULL) {
		cerr << "Cannot open file " << file2.c_str() << endl;
		exit(1);
	}
	kseq_t *seq1 = kseq_init(fp1);
	kseq_t *seq2 = kseq_init(fp2);
	int l1;
	int l2;
#pragma omp parallel private(l1, l2)
	for (string sequence1;;) {
		string sequence2;
		//TODO sanity check if names to make sure both names are the same
		string name1;
		string name2;
		string qual1;
		string qual2;
#pragma omp critical(sequence)
		{
			l1 = kseq_read(seq1);
			sequence1 = string(seq1->seq.s, seq1->seq.l);
			name1 = string(seq1->name.s, seq1->name.l);
			qual1 = string(seq1->qual.s, seq1->qual.l);

			l2 = kseq_read(seq2);
			sequence2 = string(seq2->seq.s, seq2->seq.l);
			name2 = string(seq2->name.s, seq2->name.l);
			qual2 = string(seq2->qual.s, seq2->qual.l);
		}
		if (l1 >= 0 && l2 >= 0) {
#pragma omp atomic
			++totalReads;
#pragma omp critical(totalReads)
			if (totalReads % 1000000 == 0) {
				cerr << "Currently Reading Read Number: " << totalReads << endl;
			}

			google::dense_hash_map<ID, unsigned> hitCounts1;
			google::dense_hash_map<ID, unsigned> hitCounts2;
			hitCounts1.set_empty_key(opt::EMPTY);
			hitCounts2.set_empty_key(opt::EMPTY);

			unsigned score1 = evaluateRead(sequence1, hitCounts1);
			unsigned score2 = evaluateRead(sequence2, hitCounts2);
			unsigned threshold1 = opt::score > 1 ? opt::score : opt::score
					* ((l1 - m_filter.getKmerSize() + 1)
							* (opt::allowMisses + 1));
			unsigned threshold2 = opt::score > 1 ? opt::score : opt::score
					* ((l2 - m_filter.getKmerSize() + 1)
							* (opt::allowMisses + 1));

			vector<ID> hits;
			vector<ID> solidHits;
			bool aboveThreshold = false;
			if (opt::inclusive) {
				aboveThreshold = score1 + score2 >= threshold1 + threshold2;
				if(aboveThreshold)
					convertToHitsOnlyOne(hitCounts1, hitCounts2, hits);
			} else {
				aboveThreshold = score1 >= threshold1 && score2 >= threshold2;
				if(aboveThreshold)
					convertToHitsBoth(hitCounts1, hitCounts2, hits);
			}

			ID idIndex = opt::EMPTY;

			idIndex = resSummary.updateSummaryData(hits);

			if (idIndex != opt::EMPTY) {
				const string &fullID =
						idIndex == opt::COLLI ? UNKNOWN : m_fullIDs.at(idIndex);
				if (opt::outputType == "fq") {
#pragma omp critical(outputFiles)
					{
						readsOutput << "@" << name1 << " " << fullID << "\n"
								<< sequence1 << "\n+\n" << qual1 << "\n" << "@"
								<< name2 << " " << fullID << "\n" << sequence2
								<< "\n+\n" << qual2 << "\n";
					}
				} else {
#pragma omp critical(outputFiles)
					{
						readsOutput << fullID << "\t" << name1 << "\t" << score1
								<< "\t" << score2 << "\n";
					}
				}
			}
		} else {
			if (l1 != -1 || l2 != -1) {
				cerr << "Terminated without getting to eof at read "
						<< totalReads << endl;
				exit(1);
			}
			break;
		}
	}
	kseq_destroy(seq1);
	kseq_destroy(seq2);
	gzclose(fp1);
	gzclose(fp2);

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

BloomMapClassifier::~BloomMapClassifier() {
}

