/*
 * BloomMapClassifier.h
 *
 *  Created on: Mar 23, 2016
 *      Author: cjustin
 */

#ifndef BLOOMMAPCLASSIFIER_H_
#define BLOOMMAPCLASSIFIER_H_

//#include <vector>
#include <string>
#include "boost/shared_ptr.hpp"
#include "Common/Options.h"
#include "Options.h"
#include <google/dense_hash_map>
#include "Common/Dynamicofstream.h"
#include <vector>
#include "DataLayer/FastaReader.h"
#include "BioBloomClassifier.h"
#include "bloomfilter/BloomMapSSBitVec.hpp"

using namespace std;

class BloomMapClassifier {
public:
	explicit BloomMapClassifier(const string &filterFile);
	void filter(const vector<string> &inputFiles);

	virtual ~BloomMapClassifier();
private:
	BloomMapSSBitVec<ID> m_filter;
	vector<string> m_fullIDs;
	google::dense_hash_map<ID, string> m_ids;

	//TODO: REFACTOR WITH BioBloomClassifier
	inline void printSingleToFile(const string &outputFileName,
			const FastqRecord &rec,
			google::dense_hash_map<string, boost::shared_ptr<Dynamicofstream> > &outputFiles) {
		if (opt::outputType == "fa") {
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
	}

//	void loadFilters(const vector<string> &filterFilePaths);
};

#endif /* BLOOMMAPCLASSIFIER_H_ */
