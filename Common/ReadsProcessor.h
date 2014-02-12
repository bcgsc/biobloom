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
	ReadsProcessor(unsigned windowSize);
	const unsigned char* prepSeq(string const &sequence, size_t position);
	const string getBases(const unsigned char* c); //for debuging purposes
	virtual ~ReadsProcessor();
private:
	//so reallocation does not have to be done
	unsigned char* fw;
	unsigned char* rv;
	const unsigned kmerSize;
	unsigned kmerSizeInBytes;
	unsigned halfSizeOfKmerInBytes;
	unsigned hangingBases; // used if k-mer is indivisible by 4
	unsigned hangingBasesExist;
};

#endif /* READSPROCESSOR_H_ */
