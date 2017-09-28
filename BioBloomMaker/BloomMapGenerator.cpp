/*
 * BloomMapGenerator.cpp
 *
 *  Created on: Mar 17, 2016
 *      Author: cjustin
 */

#include "BloomMapGenerator.h"
#include <zlib.h>
#include <stdio.h>
#include "Common/kseq.h"
#include "SpacedSeedIndex.h"
#include "SimMat.h"
KSEQ_INIT(gzFile, gzread)

/*
 * If size is set to 0 (or not set) a filter size will be estimated based
 * on the number of k-mers in the file
 */
BloomMapGenerator::BloomMapGenerator(vector<string> const &filenames,
		unsigned kmerSize, size_t numElements = 0) :
		m_kmerSize(kmerSize), m_expectedEntries(numElements), m_totalEntries(0), m_fileNames(
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
				cerr << "file " << m_fileNames[i] << " cannot be opened"
						<< endl;
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

/* Generate the bloom filter to the output filename
 */
void BloomMapGenerator::generate(const string &filePrefix, double fpr) {

	std::ofstream idFile;
	idFile.open((filePrefix + "_ids.txt").c_str());
	cerr << "Outputting IDs file: " << filePrefix + "_ids.txt" << endl;
	writeIDs(idFile, m_headerIDs);
	if (opt::colliIDs) {
		//if using prebuilt tree
		//load in newick file
		// assume tree format is as follows:
		// full ID is a string corresponding to fullIDs used
		// internal nodes do not yet have collision IDs

		if (fexists(filePrefix + ".newick")) {
			std::shared_ptr<BiRC::treelib::Tree> tree(
					BiRC::treelib::parse_newick_file(filePrefix + ".newick"));
			cerr << "Loading tree for collision IDs" << endl;

			//perform breadth first search of tree
			//assign largest IDs to top of tree
			assignIDBFS(*tree, tree->size());
			//starting at each child node assign collision ID groups
			setColliIds(*tree, tree->size(), idFile);

			//load ordered list of pairwise collisionsIDs
			cerr << "Loading ordered list of collision pairs" << endl;
			assert(fexists(filePrefix + "_orderedIDs.txt"));
			std::ifstream orderedIDFile;
			orderedIDFile.open((filePrefix + "_orderedIDs.txt").c_str(),
					ifstream::in);
			vector<pair<ID, ID> > orderedPairs;
			string line;
			getline(orderedIDFile, line);
			while (orderedIDFile.good()) {
				stringstream converter(line);
				ID id1, id2;
				converter >> id1;
				converter >> id2;
				orderedPairs.push_back(pair<ID, ID>(id1, id2));
				getline(orderedIDFile, line);
			}
			//assign remaining IDs to pairwise collisionIDs and assign groups
			setPairs(*tree, orderedPairs, tree->size(), idFile);

		} else {
			cerr << "Computing Similarity" << endl;
			m_colliIDs = generateGroups(idFile);
		}
	}
	idFile.close();

	MIBloomFilter<ID> * bloomMapBV;

	if (opt::sseeds.empty()) {
		bloomMapBV = generateBV(fpr);
	} else {
		cerr << "Spaced Seeds Detected" << endl;
		vector<vector<unsigned> > ssVal = parseSeedString(opt::sseeds);
		bloomMapBV = generateBV(fpr, &ssVal);
	}
	if(opt::colliIDs){
		bloomMapBV->setType(MIBloomFilter<ID>::MIBFCOLL);
	}

	cerr << "Populating values of multi index Bloom filter" << endl;

	//populate values in bitvector (bloomFilter)
	if (opt::idByFile) {
		//overwrite thread usage
		if (opt::colliAnalysis) {
			cerr << "Setting thread number to: " << m_fileNames.size() << endl;
			omp_set_num_threads(m_fileNames.size());
			Matrix colliMat(m_fileNames.size(), m_fileNames.size(), 0);
#pragma omp parallel for
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				int l;
				for (;;) {
					l = kseq_read(seq);
					if (l >= 0) {
						//k-merize with rolling hash insert into multi index Bloom filter
						loadSeq(*bloomMapBV, string(seq->seq.s, seq->seq.l),
								m_headerIDs[m_fileNames[i].substr(
										m_fileNames[i].find_last_of("/") + 1)],
								colliMat);
					} else {
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
			}
			omp_set_num_threads(opt::threads);
			cerr << "Outputting matrix" << endl;
			std::ofstream matFile;
			matFile.open((filePrefix + "_matrix.tsv").c_str());
			for (unsigned i = 0; i < colliMat.size1(); i++) {
				for (unsigned j = 0; j < colliMat.size2(); j++) {
					if (j == 0) {
						matFile << colliMat(i, j);
					} else {
						matFile << "\t" << colliMat(i, j);
					}
				}
				matFile << endl;
			}
			matFile.close();
			exit(0);
		} else {
#pragma omp parallel for
			for (unsigned i = 0; i < m_fileNames.size(); ++i) {
				gzFile fp;
				fp = gzopen(m_fileNames[i].c_str(), "r");
				kseq_t *seq = kseq_init(fp);
				int l;
				for (;;) {
					l = kseq_read(seq);
					if (l >= 0) {
						//k-merize with rolling hash insert into multi index Bloom filter
						loadSeq(*bloomMapBV, string(seq->seq.s, seq->seq.l),
								m_headerIDs[m_fileNames[i].substr(
										m_fileNames[i].find_last_of("/") + 1)]);
					} else {
						break;
					}
				}
				kseq_destroy(seq);
				gzclose(fp);
			}
		}
	} else {
#pragma omp parallel for
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			kseq_t *seq = kseq_init(fp);
			int l;
			for (;;) {
				l = kseq_read(seq);
				if (l >= 0) {
					//k-merize with rolling hash insert into multi index Bloom filter
					loadSeq(*bloomMapBV, string(seq->seq.s, seq->seq.l),
							m_headerIDs[seq->name.s]);
				} else {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
	}

	cerr << "Storing filter" << endl;

	//save filter
	bloomMapBV->store(filePrefix + ".bf");
	delete (bloomMapBV);
}

/*
 * Populates bit vector and returns a compressed multi index Bloom filter without values set
 */
inline MIBloomFilter<ID> * BloomMapGenerator::generateBV(double fpr,
		const vector<vector<unsigned> > * ssVal) {
	size_t filterSize = 0;
	if (ssVal != NULL)
		filterSize = calcOptimalSize(m_expectedEntries, opt::sseeds.size(),
				fpr);
	else
		filterSize = calcOptimalSize(m_expectedEntries, opt::hashNum, fpr);

	cerr << "Bit vector Size: " << filterSize << endl;
	size_t uniqueCounts = 0;

	sdsl::bit_vector bv(filterSize);

	cerr << "Populating initial bit vector" << endl;

	//populate sdsl bitvector (bloomFilter)
	if (opt::idByFile) {
#pragma omp parallel for
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			kseq_t *seq = kseq_init(fp);
			int l;
			size_t colliCounts = 0;
			size_t totalCount = 0;
			for (;;) {
				l = kseq_read(seq);
				if (l >= 0) {
					const string &seqStr = string(seq->seq.s, seq->seq.l);
					//k-merize with rolling hash insert into multi index Bloom filter
					if (ssVal != NULL)
						colliCounts += loadSeq(bv, seqStr, *ssVal);
					else
						colliCounts += loadSeq(bv, seqStr);
					totalCount += seq->seq.l - m_kmerSize + 1;
				} else {
					break;
				}
			}
#pragma omp atomic
			uniqueCounts += totalCount - colliCounts;
			kseq_destroy(seq);
			gzclose(fp);
		}
	} else {
#pragma omp parallel for
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			kseq_t *seq = kseq_init(fp);
			int l;
			for (;;) {
				l = kseq_read(seq);
				if (l >= 0) {
					const string &seqStr = string(seq->seq.s, seq->seq.l);
					//k-merize with rolling hash insert into multi index Bloom filter
					size_t colliCounts = 0;
					if (ssVal != NULL)
						colliCounts += loadSeq(bv, seqStr, *ssVal);
					else
						colliCounts += loadSeq(bv, seqStr);
#pragma omp atomic
					uniqueCounts += seq->seq.l - m_kmerSize + 1 - colliCounts;
				} else {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
	}
	cerr << "Approximate number of unique entries in filter: " << uniqueCounts
			<< endl;
	return new MIBloomFilter<ID>(m_expectedEntries, fpr, opt::hashNum,
			opt::kmerSize, bv, uniqueCounts, opt::sseeds);
}

inline vector<boost::shared_ptr<google::dense_hash_map<ID, ID> > > BloomMapGenerator::generateGroups(
		std::ofstream &file) {
	//compute binary tree based off of reference sequences
	//TODO make mini-spaced seed value not hard-coded
	string miniSeed = "11110000100010010101001000100001111";
	SpacedSeedIndex<ID> ssIdx(miniSeed, m_headerIDs.size());

#pragma omp parallel for
	for (unsigned i = 0; i < m_fileNames.size(); ++i) {
		//table for flagging occupancy of each readID/thread being inserted for each spaced seed
		//faster than using table directly
		uint64_t *occ = (uint64_t *) malloc(
				((ssIdx.size() + 63) / 64) * sizeof(uint64_t));
		for (unsigned i = 0; i < (ssIdx.size() + 63) / 64; ++i)
			occ[i] = 0;
		gzFile fp;
		fp = gzopen(m_fileNames[i].c_str(), "r");
		kseq_t *seq = kseq_init(fp);
		int l;
		for (;;) {
			l = kseq_read(seq);
			if (l >= 0) {
				const string &seqStr = string(seq->seq.s, seq->seq.l);
				if (opt::idByFile) {
					ssIdx.insert(
							m_headerIDs[m_fileNames[i].substr(
									m_fileNames[i].find_last_of("/") + 1)] - 1,
							seqStr, occ);
				} else {
					for (unsigned i = 0; i < (ssIdx.size() + 63) / 64; ++i)
						occ[i] = 0;
					ssIdx.insert(m_headerIDs[seq->name.s] - 1, seqStr, occ);
				}
			} else {
				kseq_destroy(seq);
				break;
			}
		}
		free(occ);
		gzclose(fp);
	}
	cerr << "Total seeds added: " << ssIdx.getCounts() << endl;
	cerr << "Index finished. Beginning similarity comparison." << endl;
	//construct similarity matrix
	SimMat<uint32_t, ID> calcSim(ssIdx, m_headerIDs.size());
	cerr << "Similarity Matrix finished. Obtaining grouping." << endl;
	return (calcSim.getGroupings(m_colliThresh, file));
}

BloomMapGenerator::~BloomMapGenerator() {
}

