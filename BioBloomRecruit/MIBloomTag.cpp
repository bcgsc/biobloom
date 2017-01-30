/*
 * MIBloomTag.cpp
 *
 *  Created on: Jan 15, 2017
 *      Author: cjustin
 */

#include <BioBloomRecruit/MIBloomTag.h>

MIBloomTag::MIBloomTag(vector<string> const &filenames,
		size_t numElements = 0) :
		m_expectedEntries(numElements), m_totalEntries(0), m_filenames(
				filenames), m_colliIDs(vector<boost::shared_ptr< IDMap > >()) {
	//Instantiate dense hash map
	m_ids.set_empty_key("");

	ID value = 0;
	size_t counts = 0;

	for (unsigned i = 0; i < m_filenames.size(); ++i) {
		gzFile fp;
		fp = gzopen(m_filenames[i].c_str(), "r");
		kseq_t *seq = kseq_init(fp);
		int l;
		for (;;) {
			l = kseq_read(seq);
			if (l >= 0) {
				string barcode = extractBarcode(seq->name.s);
				if (m_ids.find(barcode) == m_ids.end()) {
					m_ids[barcode] = ++value;
					m_colliIDs[m_ids[barcode]] = boost::shared_ptr<IDMap>(new IDMap());
					m_colliIDs[m_ids[barcode]]->set_empty_key(opt::EMPTY);
				}
				counts += seq->seq.l - opt::kmerSize + 1;
			} else {
				kseq_destroy(seq);
				break;
			}
		}
		gzclose(fp);
	}

	//estimate number of k-mers
	if (m_expectedEntries == 0) {
		m_expectedEntries = counts;
	}

//assume each file one line per file
	cerr << "Expected number of elements: " << m_expectedEntries << endl;
}

MIBloomFilter<ID> MIBloomTag::generate(const string &filePrefix, double fpr,
		BloomFilterN &filterSub, const string &file1, const string &file2) {

	//TODO REMOVE
	cerr<< filePrefix << endl;
	cerr<< file1 << endl;
	cerr<< file2 << endl;

	//setup bloom filter
	MIBloomFilter<ID> filter(m_expectedEntries, fpr, opt::hashNum, opt::kmerSize);

	size_t baitFilterElements = 0;

	for (unsigned i = 0; i < m_filenames.size(); ++i) {
		gzFile fp;
		fp = gzopen(m_filenames[i].c_str(), "r");
		kseq_t *seq = kseq_init(fp);
		int l;
#pragma omp parallel private(l)
		for (;;) {
#pragma omp critical(kseq_read)
			{
				l = kseq_read(seq);
			}
			if (l >= int(opt::kmerSize)) {
#pragma omp atomic
				baitFilterElements += l - opt::kmerSize + 1;
			} else {
				break;
			}
		}
		kseq_destroy(seq);
		gzclose(fp);
	}

	if (filterSub.getKmerSize() != opt::kmerSize) {
		cerr
				<< "Error: Subtraction filter's different from current filter's k-mer size."
				<< endl;
		exit(1);
	}

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	//extract index from header
	loadFilter(filter);

	exit(0);

//	vector<boost::shared_ptr<ReadsProcessor> > procs(opt::threads);
//	//each thread gets its own thread processor
//	for (unsigned i = 0; i < opt::threads; ++i) {
//		procs[i] = boost::shared_ptr<ReadsProcessor>(
//				new ReadsProcessor(opt::kmerSize));
//	}
//
//	for (unsigned i = 0; i < opt::progItrns; ++i) {
//		size_t totalReads = 0;
//		cerr << "Iteration " << i + 1 << endl;
//		gzFile fp1;
//		fp1 = gzopen(m_filenames[i].c_str(), "r");
//		kseq_t *seq1 = kseq_init(fp1);
//
//		gzFile fp2;
//		fp2 = gzopen(m_filenames[i].c_str(), "r");
//		kseq_t *seq2 = kseq_init(fp2);
//
//#pragma omp parallel
//		for (int l1;;) {
//			int l2;
//			//TODO optimize potential minimize reallocaiton of tempStrs.
//			char * tempStr1;
//			char * tempStr2;
//
//			if (m_totalEntries >= m_expectedEntries) {
//#pragma omp critical(breakClose)
//				{
//					cerr << "K-mer threshold reached at read " << totalReads
//							<< " iteration " << i + 1 << endl;
//				}
//			}
//
//#pragma omp critical(sequence)
//			{
//				l1 = kseq_read(seq1);
//				if (l1 >= 0) {
//					tempStr1 = new char[l1 + 1];
//					strcpy(tempStr1, seq1->seq.s);
//				}
//				l2 = kseq_read(seq2);
//				if (l2 >= 0) {
//					tempStr2 = new char[l2 + 1];
//					strcpy(tempStr2, seq2->seq.s);
//				}
//			}
//
//
//#pragma omp critical(totalReads)
//			{
//				if (totalReads++ % 1000000 == 0) {
//					cerr << "Currently Reading Read Number: " << totalReads
//							<< " ; Unique k-mers Added: " << m_totalEntries
//							<< endl;
//				}
//			}
//
//			if (l1 >= opt::kmerSize || l2 >= opt::kmerSize) {
//					size_t numKmers1 = l1 > opt::kmerSize ? l1 - opt::kmerSize + 1 : 0;
//					size_t numKmers2 = l2 > opt::kmerSize ? l2 - opt::kmerSize + 1 : 0;
//					vector<vector<size_t> > hashValues1(numKmers1);
//					vector<vector<size_t> > hashValues2(numKmers2);
//					if (numKmers1 > opt::score
//							&& (SeqEval::evalRead(rec1, opt::kmerSize, filter,
//									opt::score, 1.0 - opt::score, opt::hashNum,
//									hashValues1, filterSub, opt::mode))) {
//						//load remaining sequences
//						for (unsigned i = 0; i < numKmers1; ++i) {
//							if (hashValues1[i].empty()) {
//								const unsigned char* currentSeq =
//										procs[omp_get_thread_num()]->prepSeq(
//												rec1.seq, i);
//								insertKmer(currentSeq, filter);
//							} else {
//								insertKmer(hashValues1[i], filter);
//							}
//						}
//						//load store second read
//						for (unsigned i = 0; i < numKmers2; ++i) {
//							const unsigned char* currentSeq =
//									procs[omp_get_thread_num()]->prepSeq(
//											rec2.seq, i);
//							insertKmer(currentSeq, filter);
//						}
//					} else if (numKmers2 > opt::score
//							&& (SeqEval::evalRead(rec2, opt::kmerSize, filter,
//									opt::score, 1.0 - opt::score, opt::hashNum,
//									hashValues2, filterSub, opt::mode))) {
//						//load remaining sequences
//						for (unsigned i = 0; i < numKmers1; ++i) {
//							if (hashValues1[i].empty()) {
//								const unsigned char* currentSeq =
//										procs[omp_get_thread_num()]->prepSeq(
//												rec1.seq, i);
//								insertKmer(currentSeq, filter);
//							} else {
//								insertKmer(hashValues1[i], filter);
//							}
//						}
//						//load store second read
//						for (unsigned i = 0; i < numKmers2; ++i) {
//							if (hashValues2[i].empty()) {
//								const unsigned char* currentSeq =
//										procs[omp_get_thread_num()]->prepSeq(
//												rec2.seq, i);
//								insertKmer(currentSeq, filter);
//							} else {
//								insertKmer(hashValues2[i], filter);
//							}
//						}
//					}
//			} else
//				break;
//		}
//		cerr << "Reads Read: " << totalReads << endl;
//	}
	return filter;
}

//MIBloomTag::~MIBloomTag() {
//	// TODO Auto-generated destructor stub
//}

