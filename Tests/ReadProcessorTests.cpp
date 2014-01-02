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
	ReadsProcessor proc3(15);
	ReadsProcessor proc4(15);

	//check unhandleable sequences (should return empty string)
	assert(proc.prepSeq("NNNA",0) == NULL);
	assert(proc1.prepSeq("AATANNNNNN",0) == NULL);

//	assert(compare(proc1.prepSeq("GTACATAAAT",0), proc2.prepSeq("ATTTATGTAC",0), 3) == 0);
//	assert(strcmp(proc.prepSeq("AAAT",0), proc0.prepSeq("AAAT",0)) == 0);
//
//	assert(!strcmp(proc1.prepSeq("AAAACGTTTT",0), proc2.prepSeq("TTTTGCAAAA",0)) == 0);
//	assert(strcmp(proc1.prepSeq("AAAACGTTTT",0), proc2.prepSeq("AAAACGTTTT",0)) == 0);
//	assert(!strcmp(proc1.prepSeq("GTACATAAAT",0), proc2.prepSeq("GTTTATGTAC",0)) == 0);
//
//	assert(strcmp(proc.prepSeq("ATTT",0), proc0.prepSeq("AAAT",0)) == 0);
//	assert(strcmp(proc.prepSeq("AAAT",0), proc0.prepSeq("AAAT",0)) == 0);

//	cout << proc.getBases(proc.prepSeq("AAAA",0)) << endl;
//	cout << proc0.getBases(proc0.prepSeq("TTTT",0)) << endl;

//	assert(!strcmp(proc3.prepSeq("GTACATAAATTAAAT",0), proc4.prepSeq("GTTTATGTACTAAAT",0)) == 0);
//	assert(!strcmp(proc3.prepSeq("GTACATAAATTAAAT",0), proc4.prepSeq("GTTTATGTACTAAAA",0)) == 0);
//	assert(!strcmp(proc3.prepSeq("GTACATAAATTAAAT",0), proc4.prepSeq("GTTTATGTACTAAAG",0)) == 0);
	cout << proc3.getBases(proc3.prepSeq("CCCCCCCCCCCCCCC",0)) << endl;
	cout << proc3.getBases(proc3.prepSeq("GGGGGGGGGGGGGGG",0)) << endl;
	cout << proc3.getBases(proc3.prepSeq("GTACATAAATTAAAT",0)) << endl;
	cout << proc3.getBases(proc3.prepSeq("GTTTATGTACTAAAT",0)) << endl;
	cout << proc3.getBases(proc3.prepSeq("ATTTAGTACATAAAC",0)) << endl;
	cout << proc3.getBases(proc3.prepSeq("AAAAAAAAAAAAAAA",0)) << endl;
	cout << proc3.getBases(proc3.prepSeq("AAAAAAAAAAAAAAA",0)) << endl;


	//check hash func consistency
//	assert(
//			CityHash64WithSeed(proc.prepSeq("ATTT", 0), 4, 0) == CityHash64WithSeed(proc0.prepSeq("AAAT", 0), 4, 0));
//
//	assert(strcmp(proc.prepSeq("tAGA",0), proc0.prepSeq("TaGA",0)) == 0);
//	assert(!strcmp(proc.prepSeq("CTAA",0), proc0.prepSeq("CTAC",0)) == 0);

	cout << "Read Processor Tests Done." << endl;
	return 0;
}
