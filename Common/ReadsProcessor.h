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
	virtual ~ReadsProcessor();
private:

	//so reallocation does not have to be done
	char* fw;
	char* rv;
	const uint16_t kmerSize;
	char* emptyResult;
	uint16_t kmerSizeInBytes;
	uint8_t hangingBases; // used if k-mer is indivisible by 4
	uint16_t halfSizeOfKmerInBytes;
};

#endif /* READSPROCESSOR_H_ */
