/*
 * BloomMapClassifier.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: cjustin
 */

#include <BloomMapClassifier.h>
#include "DataLayer/FastaReader.h"
#include "bloomfilter/RollingHashIterator.h"
#include "Options.h"
//#include "ResultsManager.h"

BloomMapClassifier::BloomMapClassifier(const string &filterFile) :
		m_filter(BloomMap < ID > (filterFile)) {
}

void BloomMapClassifier::filter(const vector<string> &inputFiles) {
	size_t totalReads = 0;

	//print out header info and initialize variables
//	ResultsManager resSummary(m_filterOrder, m_inclusive);

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
					double score = double(nonZeroCount)/double(totalCount);

					if(opt::score < score && value != opt::COLLI)
					{
						cout << rec.id << " " << score << " " << value << endl;
					}
					else{
						cout << rec.id << " " << score << " NO_MATCH" << endl;
					}
				}
			} else
				break;
		}
	}
}

BloomMapClassifier::~BloomMapClassifier() {
}

