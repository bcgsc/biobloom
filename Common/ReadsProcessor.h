/*
 * ReadsProcessor.h
 *
 *  Created on: Aug 8, 2012
 *      Author: cjustin
 */

#ifndef READSPROCESSOR_H_
#define READSPROCESSOR_H_
#include <string>
#include <stdint.h>

using namespace std;

class ReadsProcessor {
public:

	ReadsProcessor(uint16_t windowSize);
	const char* prepSeq(string const &sequence, size_t position);
//	const string &prepSeqAmbigPos(string const &sequence, vector<size_t> &ambigPos);
	virtual ~ReadsProcessor();
private:

	uint16_t kmerSize;
	uint16_t kmerSizeInBytes;
	uint16_t halfSizeOfKmerInBytes;
};

#endif /* READSPROCESSOR_H_ */
