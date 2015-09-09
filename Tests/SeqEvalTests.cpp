#include "Common/SeqEval.h"
#include "Common/BloomFilter.h"
#include "Common/ReadsProcessor.h"
#include <string>
#include <iostream>

using namespace std;

int main()
{
	const size_t filterBits = 8000;
	const unsigned k = 4;
	const unsigned numHashes = 2;
	FastqRecord seq1, seq2;
	seq1.id = "seq1"; seq1.seq = "NNNNNAAAAA";
	seq2.id = "seq2"; seq2.seq = "AAAAANNNNN";

	cerr << "Loading 'NNNNNAAAAA' into Bloom filter with k=4..." << endl;

	BloomFilter bloom(filterBits, numHashes, k);
	ReadsProcessor proc(k);

	for (size_t i = 0; i < seq1.size(); ++i) {
		//string kmerStr = seq1.seq.substr(i, k);
		const unsigned char* kmer = proc.prepSeq(seq1.seq, i);
		if (kmer != NULL)
			bloom.insert(kmer);
	}

	cerr << "Query 'AAAAANNNNN' with min_match_len=5 returns TRUE... ";

	bool hit = evalRead(seq2, k, bloom, 5, 0, SeqEval::EVAL_MIN_MATCH_LEN);

	if (hit)
		cerr << "PASSED" << endl;
	else
		cerr << "FAILED" << endl;

	cerr << "Query 'AAAAANNNNN' with min_match_len=6 returns FALSE... ";

	hit = evalRead(seq2, k, bloom, 6, 0, SeqEval::EVAL_MIN_MATCH_LEN);

	if (!hit)
		cerr << "PASSED" << endl;
	else
		cerr << "FAILED" << endl;

}
