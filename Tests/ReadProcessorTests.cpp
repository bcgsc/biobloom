/*
 * readProcessorTests.cpp
 *
 *  Created on: Sep 26, 2012
 *      Author: cjustin
 */

#include "Common/ReadsProcessor.h"
#include <assert.h>
#include <string>
#include <iostream>
#include "Common/city.h"
#include <stdio.h>
#include <string.h>

using namespace std;



int main(int argc, char **argv) {

	int16_t kmerSize = 4;

	ReadsProcessor proc(kmerSize);
	ReadsProcessor proc0(kmerSize);
	ReadsProcessor proc1(10);
	ReadsProcessor proc2(10);

	assert(strcmp(proc1.prepSeq("GTACATAAAT",0), proc2.prepSeq("ATTTATGTAC",0)) == 0);
	assert(strcmp(proc.prepSeq("AAAT",0), proc0.prepSeq("AAAT",0)) == 0);

	assert(!strcmp(proc1.prepSeq("AAAACGTTTT",0), proc2.prepSeq("TTTTGCAAAA",0)) == 0);
	assert(strcmp(proc1.prepSeq("AAAACGTTTT",0), proc2.prepSeq("AAAACGTTTT",0)) == 0);
	assert(!strcmp(proc1.prepSeq("GTACATAAAT",0), proc2.prepSeq("GTTTATGTAC",0)) == 0);

	assert(strcmp(proc.prepSeq("ATTT",0), proc0.prepSeq("AAAT",0)) == 0);
	assert(strcmp(proc.prepSeq("AAAT",0), proc0.prepSeq("AAAT",0)) == 0);

	//check hash func consistency
	assert(
			CityHash64WithSeed(proc.prepSeq("ATTT", 0), 4, 0) == CityHash64WithSeed(proc0.prepSeq("AAAT", 0), 4, 0));
	assert((string)proc.prepSeq("ATTT",0) == (string)proc0.prepSeq("AAAT",0));

	assert(strcmp(proc.prepSeq("tAGA",0), proc0.prepSeq("TaGA",0)) == 0);
	assert(!strcmp(proc.prepSeq("CTAA",0), proc0.prepSeq("CTAC",0)) == 0);

	cout << "Read Processor Tests Done." << endl;
	return 0;
}
