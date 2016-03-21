/*
 * BloomMapGenerator.h
 *
 *  Created on: Mar 17, 2016
 *      Author: cjustin
 */

#ifndef BLOOMMAPGENERATOR_H_
#define BLOOMMAPGENERATOR_H_

#include <vector>
#include <string>
#include <stdint.h>
//#include <google/dense_hash_map>

using namespace std;

class BloomMapGenerator {
public:
	explicit BloomMapGenerator(vector<string> const &filenames,
			unsigned kmerSize, unsigned hashNum, size_t numElements);

	size_t generate(const string &filename);

	virtual ~BloomMapGenerator();
private:

	unsigned m_kmerSize;
	unsigned m_hashNum;
	size_t m_expectedEntries;
	size_t m_totalEntries;
	vector<string> m_fileNames;
	//google::dense_hash_map<string,ID> m_headerIDs;
};

#endif /* BLOOMMAPGENERATOR_H_ */
