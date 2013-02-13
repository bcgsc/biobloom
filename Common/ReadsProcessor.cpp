/*
 * ReadsProcessor.cpp
 * Contains methods for formatting sequences to place into bloom filter
 *
 *  Created on: Aug 8, 2012
 *      Author: cjustin
 */
#include "ReadsProcessor.h"
#include <cassert>
#include <iostream>

//TODO: write unit tests for this class

/*
 * Needed for use of optimized char* returning prepSeq
 */
ReadsProcessor::ReadsProcessor(uint16_t windowSize) :
		kmerSize(windowSize), emptyResult(""), halfKmerSize(
				windowSize / 2 + windowSize % 2) {
	outputFwd.resize(windowSize);
	outputRev.resize(windowSize);
}

/* Prepares DNA sequence for insertion into bloom filter by:
 * - Turning all lower-case sequences to upper-case
 * - Also looks into reverse compliment version and returns consistently
 *   that which is smaller (convention used A<C<G<T, i.e. alphabetical)
 * - Converts input to empty string if any character other than ATCG is found
 *
 * requires a start position
 */
const string &ReadsProcessor::prepSeq(string const &sequence, size_t position) {

	uint16_t outputIndex = 0;
	size_t revIndex = position + kmerSize - 1;

	// determines which compliment to use
	// parse through string converting and checking for lower-case and non ATCG characters
	// half the length of seq because after this point will be palindromic and is not worth checking
	for (size_t index = position; outputIndex < halfKmerSize; ++index) {
		//modify the forward
		switch (sequence[index]) {
		case 'A':
		case 'a':
			outputFwd[outputIndex] = 'A';
			break;
		case 'C':
		case 'c':
			outputFwd[outputIndex] = 'C';
			break;
		case 'G':
		case 'g':
			outputFwd[outputIndex] = 'G';
			break;
		case 'T':
		case 't':
			outputFwd[outputIndex] = 'T';
			break;
		default:
			return emptyResult;
		}

		//modify the reverse
		switch (sequence[revIndex]) {
		case 'A':
		case 'a':
			outputRev[outputIndex] = 'T';
			break;
		case 'C':
		case 'c':
			outputRev[outputIndex] = 'G';
			break;
		case 'G':
		case 'g':
			outputRev[outputIndex] = 'C';
			break;
		case 'T':
		case 't':
			outputRev[outputIndex] = 'A';
			break;
		default:
			return emptyResult;
		}

		//compare and convert if not already established
		//forward is smaller
		if (outputFwd[outputIndex] < outputRev[outputIndex]) {
			//finish off sequence
			++index;
			++outputIndex;
			for (; outputIndex < kmerSize; ++index) {
				switch (sequence[index]) {
				case 'A':
				case 'a':
					outputFwd[outputIndex] = 'A';
					break;
				case 'C':
				case 'c':
					outputFwd[outputIndex] = 'C';
					break;
				case 'G':
				case 'g':
					outputFwd[outputIndex] = 'G';
					break;
				case 'T':
				case 't':
					outputFwd[outputIndex] = 'T';
					break;
				default:
					return emptyResult;
				}
				++outputIndex;
			}
			return outputFwd;
		}
		//reverse is smaller
		else if (outputFwd[outputIndex] > outputRev[outputIndex]) {
			//finish off sequence
			--revIndex;
			++outputIndex;
			for (; outputIndex < kmerSize; --revIndex) {
				switch (sequence[revIndex]) {
				case 'A':
				case 'a':
					outputRev[outputIndex] = 'T';
					break;
				case 'C':
				case 'c':
					outputRev[outputIndex] = 'G';
					break;
				case 'G':
				case 'g':
					outputRev[outputIndex] = 'C';
					break;
				case 'T':
				case 't':
					outputRev[outputIndex] = 'A';
					break;
				default:
					return emptyResult;
				}
				++outputIndex;
			}
			return outputRev;
		}
		++outputIndex;
		--revIndex;
	}
	//palamdromic
	return outputFwd;
}

//const string &ReadsProcessor::prepSeqAmbigPos(string const &sequence, vector<size_t> &ambigPos) {
//
//	size_t outputIndex = 0;
//	size_t revIndex = kmerSize - 1;
//
//	// determines which compliment to use
//	// parse through string converting and checking for lower-case and non ATCG characters
//	// half the length of seq because after this point will be palindromic and is not worth checking
//	for (size_t index = 0; outputIndex < halfKmerSize; ++index) {
//		//modify the forward
//		switch (sequence[index]) {
//		case 'A':
//		case 'a':
//			outputFwd[outputIndex] = 'A';
//			break;
//		case 'C':
//		case 'c':
//			outputFwd[outputIndex] = 'C';
//			break;
//		case 'G':
//		case 'g':
//			outputFwd[outputIndex] = 'G';
//			break;
//		case 'T':
//		case 't':
//			outputFwd[outputIndex] = 'T';
//			break;
//		default:
//			outputFwd[outputIndex] = 'N';
//			ambigPos.push_back(outputIndex);
//			break;
//		}
//
//		//modify the reverse
//		switch (sequence[revIndex]) {
//		case 'A':
//		case 'a':
//			outputRev[outputIndex] = 'T';
//			break;
//		case 'C':
//		case 'c':
//			outputRev[outputIndex] = 'G';
//			break;
//		case 'G':
//		case 'g':
//			outputRev[outputIndex] = 'C';
//			break;
//		case 'T':
//		case 't':
//			outputRev[outputIndex] = 'A';
//			break;
//		default:
//			outputFwd[outputIndex] = 'N';
//			ambigPos.push_back(outputIndex);
//			break;
//		}
//
//		//compare and convert if not already established
//		//forward is smaller
//		if (outputFwd[outputIndex] < outputRev[outputIndex]) {
//			//finish off sequence
//			++index;
//			++outputIndex;
//			for (; outputIndex < kmerSize; ++index) {
//				switch (sequence[index]) {
//				case 'A':
//				case 'a':
//					outputFwd[outputIndex] = 'A';
//					break;
//				case 'C':
//				case 'c':
//					outputFwd[outputIndex] = 'C';
//					break;
//				case 'G':
//				case 'g':
//					outputFwd[outputIndex] = 'G';
//					break;
//				case 'T':
//				case 't':
//					outputFwd[outputIndex] = 'T';
//					break;
//				default:
//					outputFwd[outputIndex] = 'N';
//					ambigPos.push_back(outputIndex);
//					break;
//				}
//				++outputIndex;
//			}
//			return outputFwd;
//		}
//		//reverse is smaller
//		else if (outputFwd[outputIndex] > outputRev[outputIndex]) {
//			//finish off sequence
//			--revIndex;
//			++outputIndex;
//			for (; outputIndex < kmerSize; --revIndex) {
//				switch (sequence[revIndex]) {
//				case 'A':
//				case 'a':
//					outputRev[outputIndex] = 'T';
//					break;
//				case 'C':
//				case 'c':
//					outputRev[outputIndex] = 'G';
//					break;
//				case 'G':
//				case 'g':
//					outputRev[outputIndex] = 'C';
//					break;
//				case 'T':
//				case 't':
//					outputRev[outputIndex] = 'A';
//					break;
//				default:
//					outputFwd[outputIndex] = 'N';
//					ambigPos.push_back(outputIndex);
//					break;
//				}
//				++outputIndex;
//			}
//			return outputRev;
//		}
//		++outputIndex;
//		--revIndex;
//		ambigPos.pop_back();
//	}
//	//palamdromic
//	return outputFwd;
//}

ReadsProcessor::~ReadsProcessor() {
}

