/*
 * MIBloomTag.cpp
 *
 *  Created on: Jan 15, 2017
 *      Author: cjustin
 */

#include <BioBloomRecruit/MIBloomTag.h>

/*
 * If size is set to 0 (or not set) a filter size will be estimated based
 * on the number of k-mers in the file
 */
MIBloomTag::MIBloomTag(vector<string> const &filenames, size_t numElements = 0, IDMap idMap) :
		m_expectedEntries(numElements), m_totalEntries(0), m_fileNames(
				filenames) {
	//Instantiate dense hash map
	m_headerIDs.set_empty_key("");

	ID value = 0;
	size_t counts = 0;

	if (opt::idByFile) {
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			if (fp == NULL) {
				cerr << "file " << m_fileNames[i] << " cannot be opened" << endl;
				exit(1);
			}
			kseq_t *seq = kseq_init(fp);
			m_headerIDs[m_fileNames[i].substr(
					m_fileNames[i].find_last_of("/") + 1)] = ++value;
			int l;
			for (;;) {
				l = kseq_read(seq);
				if (l >= 0) {
					counts += seq->seq.l - m_kmerSize + 1;
				} else {
					kseq_destroy(seq);
					break;
				}
			}
			gzclose(fp);
		}
	} else {
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			kseq_t *seq = kseq_init(fp);
			int l;
			for (;;) {
				l = kseq_read(seq);
				if (l >= 0) {
					m_headerIDs[seq->name.s] = ++value;
					counts += seq->seq.l - m_kmerSize + 1;
				} else {
					kseq_destroy(seq);
					break;
				}
			}
			gzclose(fp);
		}
	}

	//estimate number of k-mers
	if (m_expectedEntries == 0) {
		m_expectedEntries = counts;
	}
	//assume each file one line per file
	cerr << "Expected number of elements: " << m_expectedEntries << endl;
}

MIBloomTag::~MIBloomTag() {
	// TODO Auto-generated destructor stub
}

