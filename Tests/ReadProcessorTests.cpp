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

using namespace std;

int main(int argc, char **argv) {

	int16_t kmerSize = 4;

	ReadsProcessor proc(kmerSize);

	assert(proc.prepSeq("GAAT",0) == "ATTC");
	assert(proc.prepSeq("AAAT",0) == proc.prepSeq("AAAT",0));

	assert(proc.prepSeq("ATTT",0) == "AAAT");
	assert(proc.prepSeq("AAAT",0) == "AAAT");
	assert("AAAT" == proc.prepSeq("AAAT",0));
	assert("AAAT" == proc.prepSeq("AAAT",0));

	string string1 = proc.prepSeq("ATTT", 0);
	string string2 = proc.prepSeq("AAAT", 0);

	assert(string1 == string2);
	//check hash func consistency
	assert(
			CityHash64WithSeed(proc.prepSeq("ATTT", 0).c_str(), 4, 0) == CityHash64WithSeed(proc.prepSeq("AAAT", 0).c_str(), 4, 0));
	assert((string)proc.prepSeq("ATTT",0) == (string)proc.prepSeq("AAAT",0));

	assert(proc.prepSeq("tAGA",0) == proc.prepSeq("TaGA",0));
	assert((string)proc.prepSeq("CTAA",0) != (string)proc.prepSeq("CTAC",0));

	cout << "Read Processor Tests Done." << endl;
	return 0;
}
