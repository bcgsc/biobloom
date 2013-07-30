/*
 * WindowedFileParser.cpp
 *
 * Currently for fasta files only and needs fasta index. Parse file line by line.
 *
 *  Created on: Jul 18, 2012
 *      Author: cjustin
 */
#include "WindowedFileParser.h"
#include <sstream>
#include "DataLayer/FastaIndex.h"

WindowedFileParser::WindowedFileParser(string const &fileName,
		uint16_t windowSize) :
		windowSize(windowSize), proc(ReadsProcessor(windowSize))
{
	fastaFileHandle.open(fileName.c_str(), ifstream::in);
	assert(fastaFileHandle);

//create in memory index
	WindowedFileParser::initializeIndex(fileName);
	currentHeader = "";
	setLocationByHeader(headers[0]);
}

const vector<string> WindowedFileParser::getHeaders() const
{
	return headers;
}

//sets the location in the file to the start of the sequence given a header
void WindowedFileParser::setLocationByHeader(string const &header)
{
	sequenceNotEnd = true;
	currentHeader = header;
	fastaFileHandle.seekg(fastaIndex[header].start, ios::beg);
	currentCharNumber = 0;
	string bufferString;
	getline(fastaFileHandle, currentString);
	while ((currentString.length() < windowSize)
			&& (currentCharNumber < fastaIndex[currentHeader].size)
			&& getline(fastaFileHandle, bufferString))
	{
		currentString += bufferString;
	}
	currentLinePos = 0;
}

const size_t WindowedFileParser::getSequenceSize(string const &header)
{
	return fastaIndex[header].size;
}

//Todo: Optimize to skip sections when finding a non ATCG character
// Also upper case conversion should occur only once.
/*
 * Return the next string in sliding window, also cleans and formats
 * sequences using ReadProcessor
 */
const string &WindowedFileParser::getNextSeq()
{
	if (currentString.length() < windowSize + currentLinePos) {
		currentString.erase(0, currentLinePos);
		currentLinePos = 0;
		//grow the sequence to match the correct window size
		//stop if there are no more lines left in fasta file
		while (fastaFileHandle.is_open()
				&& (currentString.length() < windowSize)
				&& (currentCharNumber < fastaIndex[currentHeader].size)
				&& getline(fastaFileHandle, bufferString))
		{
			currentString += bufferString;
		}
		//if there is not enough sequence for a full kmer
		if (currentString.length() < windowSize) {
			sequenceNotEnd = false;
			return currentString;
		}
	}

	return proc.prepSeq(currentString, currentLinePos++);
}

const bool WindowedFileParser::notEndOfSeqeunce()
{
	return sequenceNotEnd;
}

/*
 * Initializes fasta index in memory. Will output index if not present.
 * Input file refers to input fasta file not the index file.
 */
void WindowedFileParser::initializeIndex(string const &fileName)
{
	string faiFile = fileName + ".fai";
	//look for fasta index file
	//append .fai to end of filename
	ifstream indexFile;
	indexFile.open(faiFile.c_str(), ifstream::in);
	if (!indexFile)
	{
		cerr << "Fasta files must be indexed. Use samtools faidx." << endl;
		exit(1);
	}

	string line;
	//read file and populate index struct
	if (indexFile.is_open()) {
		while (getline(indexFile, line)) {
			stringstream ss(line);
			string header;
			FastaIndexValue value;
			ss >> header >> value.size >> value.start >> value.bpPerLine
					>> value.charsPerLine;
			value.index = headers.size();
			fastaIndex[header] = value;
			headers.push_back(header);
		}
		indexFile.close();
	}
}

WindowedFileParser::~WindowedFileParser()
{
	fastaFileHandle.close();
}

