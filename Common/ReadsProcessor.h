/*
 * ReadsProcessor.h
 *
 *  Created on: Aug 8, 2012
 *      Author: cjustin
 */

#ifndef READSPROCESSOR_H_
#define READSPROCESSOR_H_
#include <string>
#include <deque>

using namespace std;

class ReadsProcessor {
public:

	ReadsProcessor(int16_t windowSize);
	const string &prepSeq(string const &sequence, size_t position);
	virtual ~ReadsProcessor();
private:
	string outputFwd; //containers preventing reallocation of mem
	string outputRev; //containers preventing reallocation of mem
	string emptyResult;
	int16_t kmerSize;
	int16_t halfKmerSize;
};

#endif /* READSPROCESSOR_H_ */
