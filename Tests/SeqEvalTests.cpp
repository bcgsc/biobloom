#include "Common/SeqEval.h"
#include <string>
#include <iostream>

using namespace std;

int main()
{
	const size_t filterBits = 8000;
	const unsigned k = 4;
	const unsigned numHashes = 2;
	string seq1 = "NNNNNAAAAA";
	string seq2 = "AAAAANNNNN";

	cerr << "Loading 'NNNNNAAAAA' into Bloom filter with k=4..." << endl;

	BloomFilter bloom(filterBits, numHashes, k);

	for (ntHashIterator i = ntHashIterator(seq1, numHashes, k);
			i != i.end(); ++i) {
		bloom.insert(*i);
	}

	cerr << "Query 'AAAAANNNNN' with min_match_len=5 returns TRUE... ";

	bool hit = SeqEval::evalRead(seq2, bloom, 5);

	if (hit)
		cerr << "PASSED" << endl;
	else
		cerr << "FAILED" << endl;

	cerr << "Query 'AAAAANNNNN' with min_match_len=6 returns FALSE... ";

	hit = SeqEval::evalRead(seq2, bloom, 6);

	if (!hit)
		cerr << "PASSED" << endl;
	else
		cerr << "FAILED" << endl;

}
