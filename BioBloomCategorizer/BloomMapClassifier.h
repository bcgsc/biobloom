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
#include "bloomfilter/BloomMap.hpp"
#include "boost/shared_ptr.hpp"
#include "Common/Options.h"

using namespace std;

class BloomMapClassifier {
public:
	explicit BloomMapClassifier(const string &filterFile);
	void filter(const vector<string> &inputFiles);

	virtual ~BloomMapClassifier();
private:
	BloomMap<ID> m_filter;
//	void loadFilters(const vector<string> &filterFilePaths);
};

#endif /* BLOOMMAPCLASSIFIER_H_ */
