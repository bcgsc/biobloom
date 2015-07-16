BioBloomtools User Manual
------
BioBloom Tools (BBT) provides the means to create filters for a given reference and then to categorizes sequences. This methodology is faster than alignment but does not provide mapping locations. BBT was initially intended to be used for pre-processing and QC applications like contamination detection, but is flexible to accommodate other purposes. This tool is intended to be a pipeline component to replace costly alignment steps.

This tool is free for academic (BCCA licence).

We have a comercial licence also, please contact us (prebstein@bccancer.bc.ca) if you wish to use the tools for comercial uses.

Any comments and suggestions can be directed to [JustinChu](https://github.com/JustinChu) or emailed to cjustin@bcgsc.ca.

######Affliations:
Canada’s Michael Smith Genome Sciences Centre, BC Cancer Agency, Vancouver BC Canada V5Z 4S6

Department of Bioinformatics, University of British Columbia, Vancouver BC V6T 1Z4

Table of Contents
======
1. [Compiling and Installing BioBloomTools](#1)
  * Dependencies
  * How to Install
2. [Generating Bloom Filters from Reference Sequence with Biobloommaker](#2)
3. [Classifying and Analyzing Sequences with Biobloomcategorizer](#3)
4. [Program Outputs](#4)
  * Biobloommaker
    * Bloom Filter File (filterID.bf)
    * Bloom Filter Info File (filterID.txt)
  * Biobloomcategorizer
    * Summary File (summary.tsv)
    * Categorized Sequence FastA/FastQ Files
5. [Program Options](#5)
  * Biobloommaker
  * Biobloomcategorizer
6. [Understanding BioBloomTools](#6)
  * About Bloom Filters
  * How false positive rates correlates to memory usage
  * How many hash functions should be used?
  * K-mer Tiling and Reduction of the effect of false positive k-mers
  * What is inside the Bloom Filter Info File?
    * Obtaining the number of unique k-mers in the reference
    * Obtaining the number of redundant k-mers in the reference
7. [Advanced options and Best Practices](#7)
  * How can I reduce my memory usage?
  * How can I make my results more sensitive?
  * How can I make my results more specific?
  * How can I make the program faster?
 
<a name="1"></a>
1. Compiling and Installing BioBloomTools
======
######Dependencies:
* GCC (tested on 4.8.4)
* Boost (tested on 1.54)
* zlibdev
* Autotools (if directly cloning from repo)

######Compilation:
If cloning directly from the repository run:
```
./autogen.sh
```
Compiling BioBloomTools should be as easy as:
```
./configure && make
```
To install BBT in a specified directory:
```
./configure --prefix=/BBT/PATH && make install
```
If your boost library headers are not in your PATH you can specify its location:
```
./configure –-with-boost=/boost/path --prefix=/BBT/PATH && make install
```

<a name="2"></a>
2. Generating Bloom Filters from Reference Sequences with Biobloommaker
======
To create bloom filters from a FastA file, the FastA file must by indexed. Indexing can
be done by programs like [samtools](https://github.com/samtools/samtools) (faidx) or [fastahack](https://github.com/ekg/fastahack).

After you have your FastA file and index, a .bf file with corresponding information text
file can be created by running the command:
```
./biobloommaker –p input input1.fasta input2.fasta
```
-p is the prefix for the files being created, it also acts as an ID for the filter.

The options above are the bare minimum options you must use to run the program, but it is possible to customize many aspects of your filter that can drastically change performance depending on your needs. See section 5 for advanced options. You can also use the -h command for a listing of the options.

The optimal size of the filter will be calculated based on the maximum false positive rate (default is 0.075) and the number of hash functions (can be set but is optimized based on FPR).

Two files will be generated binary Bloom filter file (.bf) and an information file in INI format (.txt). The information file must be kept with the .bf file to provide all the needed information to run the categorization.





