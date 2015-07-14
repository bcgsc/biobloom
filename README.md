BioBloomtools User Manual
------
BioBloom Tools (BBT) provides the means to create filters for a given reference and then to categorizes sequences. This methodology is faster than alignment but does not provide mapping locations. BBT was initially intended to be used for pre-processing and QC applications like contamination detection, but is flexible to accommodate other purposes. This tool is intended to be a pipeline component to replace costly alignment steps.

This tool is free for academic use, please contact us if you wish to use the tools for comercial uses.

Any comments and suggestions can be directed to [JustinChu](https://github.com/JustinChu) or emailed to cjustin@bcgsc.ca.

######Affliations:
Canada’s Michael Smith Genome Sciences Centre, BC Cancer Agency, Vancouver BC Canada V5Z 4S6
Department of Bioinformatics, University of British Columbia, Vancouver BC V6T 1Z4

Table of Contents
======
1. Compiling and Installing BioBloomTools
  * Dependencies
  * How to Install
2. Generating Bloom Filters from Reference Sequence with Biobloommaker
3. Classifying and Analyzing Sequences with Biobloomcategorizer
4. Program Outputs
  * Biobloommaker
    * Bloom Filter File (filterID.bf)
    * Bloom Filter Info File (filterID.txt)
  * Biobloomcategorizer
    * Summary File (summary.tsv)
    * Categorized Sequence FastA/FastQ Files
5. Program Options
  * Biobloommaker
  * Biobloomcategorizer
6. Understanding BioBloomTools
  * About Bloom Filters
  * How false positive rates correlates to memory usage
  * How many hash functions should be used?
  * K-mer Tiling and Reduction of the effect of false positive k-mers
  * What is inside the Bloom Filter Info File?
    * Obtaining the number of unique k-mers in the reference
    * Obtaining the number of redundant k-mers in the reference
7. Advanced options and Best Practices
  * How can I reduce my memory usage?
  * How can I make my results more sensitive?
  * How can I make my results more specific?
  * How can I make the program faster?
