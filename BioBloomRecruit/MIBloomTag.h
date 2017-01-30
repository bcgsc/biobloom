/*
 * MIBloomTag.h
 *
 *  Created on: Jan 15, 2017
 *      Author: cjustin
 */

#ifndef MIBLOOMTAG_H_
#define MIBLOOMTAG_H_

#include "bloomfilter/MIBloomFilter.hpp"
#include "bloomfilter/BloomFilterN.hpp"
#include <vector>
#include <string>
#include <google/dense_hash_map>
#if _OPENMP
# include <omp.h>
#endif
#include "Common/Options.h"
#include "DataLayer/kseq.h"
//#include "WindowedFileParser.h"
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

#include <boost/shared_ptr.hpp>
#include "bloomfilter/ntHashIterator.hpp"

using namespace std;

class MIBloomTag {
public:
	typedef typename google::dense_hash_map<string, ID> StrIDMap;
	typedef typename google::dense_hash_map<ID, ID> IDMap;
	explicit MIBloomTag(vector<string> const &filenames,
			size_t numElements);

	MIBloomFilter<ID> generate(const string &filePrefix, double fpr,
			BloomFilterN &filterSub, const string &file1, const string &file2);

private:
	size_t m_expectedEntries;
	size_t m_totalEntries;
	vector<string> m_filenames;
	vector<boost::shared_ptr< IDMap > > m_colliIDs;
	StrIDMap m_ids;

	//helper methods
	inline uint64_t loadSeq(MIBloomFilter<ID> &filter, const string& seq,
			ID value) {
		uint64_t count = 0;
		for (ntHashIterator itr(seq, opt::kmerSize, opt::kmerSize);
				itr != itr.end(); ++itr) {
			count += filter.insertAndCheck(*itr, value, m_colliIDs);
		}
		return count;
	}

	inline string extractBarcode(const string &header){
		//TODO remove magic number
		return header.substr(header.find_last_of("_"), 16);
	}

	inline void loadFilter(MIBloomFilter<ID> &bf){
		unsigned threadNum = omp_get_max_threads();
		vector< boost::shared_ptr<vector<size_t> > > tempHashValues(threadNum);

		for(unsigned i = 0; i < threadNum; ++i){
			tempHashValues[i] = boost::shared_ptr<vector<size_t> >(new vector<size_t>(opt::hashNum));
		}

		for (unsigned i = 0; i < m_filenames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_filenames[i].c_str(), "r");
			if (fp == NULL) {
				cerr << "file " << m_filenames[i] << " cannot be opened" << endl;
				exit(1);
			}
			kseq_t *seq = kseq_init(fp);
			int l;
			string tempStr;
			string header;
#pragma omp parallel private(l, tempStr, header)
			for (;;) {
#pragma omp critical(kseq_read)
				{
					l = kseq_read(seq);
					if (l >= 0) {
						tempStr = string(seq->seq.s, seq->seq.l);
						header = string(seq->name.s, seq->name.l);
					}
				}

				if (l >= 0) {
					ID barcodeID = m_ids[extractBarcode(header)];
#pragma omp atomic update
					m_totalEntries += loadSeq(bf, tempStr, barcodeID);
				} else {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
	}

	inline void writeIDs(std::ofstream &file, google::dense_hash_map<string,ID> headerIDs) {
		assert(file);
		for (google::dense_hash_map<string, ID>::iterator itr =
				headerIDs.begin(); itr != headerIDs.end(); ++itr) {
			file << itr->second << "\t" << itr->first << "\n";
			assert(file);
		}
	}

};

#endif /* MIBLOOMTAG_H_ */
