# BioBloom Tools User Manual

BioBloom Tools (BBT) provides the means to create filters for a given reference and then to categorize sequences. This methodology is faster than alignment but does not provide mapping locations. BBT was initially intended to be used for pre-processing and QC applications like contamination detection, but is flexible to accommodate other purposes. This tool is intended to be a pipeline component to replace costly alignment steps.

Relevant [paper](http://bioinformatics.oxfordjournals.org/content/30/23/3402.long):

> **BioBloom tools: fast, accurate and memory-efficient host species sequence screening using bloom filters.**<br/>
> _Justin Chu, Sara Sadeghi, Anthony Raymond, Shaun D. Jackman, Ka Ming Nip, Richard Mar, Hamid Mohamadi, Yaron S. Butterfield, A. Gordon Robertson, Inanç Birol_<br/>
> Bioinformatics 2014; **30** (23): 3402-3404.<br/>
> doi: [10.1093/bioinformatics/btu558](http://dx.doi.org/10.1093/bioinformatics/btu558)<br/>
> PMID: [25143290](http://www.ncbi.nlm.nih.gov/pubmed/25143290)

This tool is free for academic use ([BCCA licence](COPYING)).

We have a commercial licence also, please contact us (prebstein at bccancer dot bc dot ca) if you wish to use the tools for commercial uses.

Any questions, comments and suggestions can be directed to [@JustinChu](https://github.com/JustinChu) or emailed to cjustin at bcgsc dot ca.

### Affiliations:
Canada’s Michael Smith Genome Sciences Centre, BC Cancer Agency, Vancouver BC Canada V5Z 4S6

Department of Bioinformatics, University of British Columbia, Vancouver BC V6T 1Z4

Table of Contents
------
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
5. [Understanding BioBloomTools](#5)
    * About Bloom Filters
    * How false positive rates correlates to memory usage
    * How many hash functions should be used?
    * K-mer Tiling and Reduction of the effect of false positive k-mers
    * What is inside the Bloom Filter Info File?
        * Obtaining the number of unique k-mers in the reference
        * Obtaining the number of redundant k-mers in the reference
    * Specifications on memory, cpu and storage requirements
6. [Advanced options and Best Practices](#6)

<a name="1"></a>
## 1. Compiling and Installing BioBloomTools

### Dependencies:
* GCC (tested on 4.8.4)
* Boost (tested on 1.54)
* zlibdev
* Autotools (if directly cloning from repo)

### Compilation:
If cloning directly from the repository make sure you get the require submodules:
```bash
git submodule update --init
```

If cloning directly from the repository run:
```bash
./autogen.sh
```
Compiling BioBloomTools should be as easy as:
```bash
./configure && make
```
To install BBT in a specified directory:
```bash
./configure --prefix=/BBT/PATH && make install
```
If your boost library headers are not in your PATH you can specify their location:
```bash
./configure –-with-boost=/boost/path --prefix=/BBT/PATH && make install
```

<a name="2"></a>
## 2. Generating Bloom Filters from Reference Sequences with Biobloommaker

To create bloom filters from a FastA file, the FastA file must by indexed. Indexing can
be done by programs like [samtools](https://github.com/samtools/samtools) (faidx) or [fastahack](https://github.com/ekg/fastahack). (As of version 2.0.13 this is no longer required)

After you have your FastA file and index, a `.bf` file with corresponding information text
file can be created by running the command:
```bash
./biobloommaker –p input input1.fasta input2.fasta
```
`-p` is the prefix for the files being created, it also acts as an ID for the filter.

The options above are the bare minimum options you must use to run the program, but it is possible to customize many aspects of your filter that can drastically change performance depending on your needs. See section 5 for advanced options. You can also use the `-h` command for a listing of the options.

The optimal size of the filter will be calculated based on the maximum false positive rate (default is 0.075) and the number of hash functions (can be set but is optimized based on FPR).

Two files will be generated binary Bloom filter file (`.bf`) and an information file in INI format (`.txt`). The information file must be kept with the `.bf` file to provide all the needed information to run the categorization.

<a name="3"></a>
## 3. Classifying and Analyzing Sequences with Biobloomcategorizer

Once you have filters created, you can use them with Biobloomcategorizer to categorize sequences. The file formats that can be used are the following: SAM, BAM, FastQ, FastA and qseq. To read BAM files samtools must be in your path. Gzip and Bz2 compression is also handled if your system has gzip and `bunzip2` installed.

Before starting make sure the listed `.bf` file is in the same directory as its corresponding information `.txt` file.

A `summary.tsv` file, `readStatus.tsv` summary file, and FastA files containing categorized reads can be generated with the following command:
```bash
./biobloomcategorizer –-fa –p /output/prefix –f "filter1.bf filter2.bf filter3.bf" inputReads1.bam.bz2 inputreads2_qseq.txt
```
`-p` is the output directory for the output files.
`-f` is the filter(s) you used to categorize the sequences. You can specify as many as
you need.

There are some advanced options open can use outlined in section 5. Notable option
one can use is the paired end mode `-e`:
```bash
./biobloomcategorizer -e –p /output/prefix –f "filter1.bf filter2.bf filter3.bf" inputReads1_1.fq inputreads1_2.fq
```
`-e` will require that both reads match when making the call about what reference they belong in.

By default `-e` will only count a read if both reads match a filter. If you want only it to count situations where only one read matches the filter then the `-i` (`--inclusive`) option can also be used.

These are general use cases you can use to run the program, but it is possible to customize many aspects of your filter that can drastically change performance depending on your needs. See [section 6](#6) for advanced options. You can also using the `-h` command for a listing on the options.

<a name="4"></a>
## 4. Program Output

### A. Biobloommaker
#### i. Bloom Filter File (`filterID.bf`)

* Simply a bit array representing the bloom filter dumped as a file. It is simple, so that it can be used in almost any system format. It is useless without its paired info file however.

#### ii. Bloom Filter Info File (`filterID.txt`)

* This is the information file of the bloom filter, containing the information like false positive rate and hash functions used. It is in human readable INI format. It is intended to be read by Biobloomcategorizer in tandem with its paired `.bf` file to perform categorization.

### B. Biobloomcategorizer
#### i. Summary File (`summary.tsv`)

* Tab separated file. Contains proportion information about reads mapping to each filter. Gives a overview of your results.
    * The `hits` column is the number of hits to this filter, regardless if unique.
    * The `shared` column a subset of hits that also hit other filter.
    * The `multiMatch` entry refers to any hits are shared at least 2 genomes.
    * The `noMatch` entry is a subset of reads that do not match any of the filters used.

#### ii. Categorized Sequence FastA/FastQ Files

* In the output directory there will be files for every filter used in addition to `multiMatch` and `noMatch` files. The reads will be categorized in these locations based on the threshold (`-m` and `-t`) values used.
* Reads outputted will have a value (e.g. `/1`) appended to the end of each ID to denote pair information about the read.

<a name="5"></a>
## 5. Understanding BioBloomTools


### A. About Bloom Filters
The whole idea of using bloom filter centers on getting the time complexity of a hash (i.e. a `O(1)` time complexity for look-ups) with a lower space requirement. This is resolved in a bloom filter by not storing the entry, but rather storing the entry’s bit signature (determined via hashing) in a bit array. However, this means there is no collision detection and all bloom filters will have some sort of false positive rate associated with them. The false positive rate must be carefully considered as it determines the expected size of the filter used. Too small of a false positive rate can mean a large bloom filter, but too large of a false positive rate could introduce too much error for the filter to be practical.

### B. How false positive rates correlates to memory usage
![false discovery rate and bits per entry](Doc/FDR_vs_Size.png)

This figure shows the relationship (assuming optimal number of hash functions have been used) between false positive rate and the space cost per entry. To use this chart divide your amount of space in bits you have to work with to the number base pairs you have in your reference.

For example, say I want the human genome (~ `3.4×10^9`) filter to fit into~3GB of memory. `(8×3×230)/4×10^9` = ~8bits per entry, meaning a filter with 2% FPR at maximum should be used.

### C. How many hash functions should be used?

![false discovery rate and bits per entry](Doc/FDR_vs_HashFunct.png)

The number of hash functions refers to the number of hash functions used by a single filter per element. In practice the approximate optimal number of hash functions will be calculated automatically by our program.

We give users the ability to change the number of hash functions because very low false positive rates will have an optimal number of associated hash functions that may be very large and may slow down classification. Also it is recommended that if running multiple filters at the same time that they all use the same number of hash functions.

### D. K-mer Tiling
We use a sliding window across each read of size k to categorize sequences. Single base overlaps that both hit a filter are unlikely to be false positives. This information is used to reduce the effect of any false positives. This concept is also used to further improve speed, we also employ a jumping k-mer (rather than sliding) heuristic that skips k k-mers when a miss is detected after a long series of adjacent hits.

### E. What is inside the Bloom Filter Info File?

#### i. Obtaining the number of unique k-mers in the reference:
Within the information txt file for each bloom filter there is a "num_entries" entry that lets you know how many unique k-mers have been added to the filter. It is a lower bound estimate due to possible false positives.

#### ii. Obtaining the number of redundant k-mers in the reference:
Within the information txt file for each bloom filter there is a `redundant_sequences` entry that lets you know how many redundant k-mers have been added to the filter. It is an upper bound estimate due to possible false positives. The `redundant_fpr` represents the probability a unique entry could be mis-classified as redundant. Thus, to get the approximate number of unique k-mers take the `redundant_fpr` value and multiply it with the `redundant_sequences` + `num_entries` and add that to the `num_entries`.

Expected total number of k-mers = (`redundant_sequences` + `num_entries`) * `redundant_fpr` + `num_entries`

### F. Specifications on memory, cpu and storage requirements

#### Memory:
Memory usage is determined by the size of the database given so there is no "optimal" memory requirements. There is a bit of overhead but the formula roughly matches this: `m=(-ln(p)/((ln(2))^2)) * n` where `m` is size of the filter in bits, `p` is the FPR and `n` is the number of elements (number of bases in reference fasta file).

If used in an job based automated cluster environment either the user can specify memory of you can infer memory usage based on the input bloom filter sizes (e.g. size of the input bf files + 100mb overhead).

#### Storage:
When creating filters:
Again, as alluded to above, storage is dependant of the input reference fasta file, following the same formula. Here is a nice sanity check: memory usage is roughly equal the sum of sizes of the raw bloom filter file used and vice versa.

#### When filtering reads:
The size of the output is proportional to the input since the results need to be stored. The contents of the output fastq files can be compressed directly as needed with the `--gz` option however (which is what `zlib` was needed in the installation).

BBT does not create temporary files so scratch space is not needed.

If used in an job based automated cluster environment where users have their own allocated storage they should make sure they have space for the output bloom filter. When categorizing reads they should make sure they have space for the output (if they want the reads `--fa` or `--fq`) which will be roughly the size of the input files since all they are doing is partitioning the reads the reads.

#### CPU:
There is no cpu minimum speed or number of cores, though it will run faster with more and faster cpus. In terms of a maximum, speed can become I/O bound quickly. When using only a few bloom filters(<5) in BBC the number of cores (>4) may not matter too much, but you will get better performance with multiple threads if more bloom filters are used at the same time. Also BBM does not yet fully support threads (though if critically needed I think I could implement this fully - contact me).

<a name="6"></a>
## 6. Advanced options and Best Practices

### General
Help dialog from the `--help` options has the most up to date information about all the options in this tool and consulting it should give you good information about the effects of some options.

### biobloommaker
```
Usage: biobloommaker -p [FILTERID] [OPTION]... [FILE]...
Usage: biobloommaker -p [FILTERID] -r 0.2 [FILE]... [FASTQ1] [FASTQ2] 
Creates a bf and txt file from a list of fasta files. The input sequences are
cut into a k-mers with a sliding window and their hash signatures are inserted
into a bloom filter.

  -p, --file_prefix=N    Filter prefix and filter ID. Required option.
  -o, --output_dir=N     Output location of the filter and filter info files.
  -h, --help             Display this dialog.
  -v  --version          Display version information.
  -t, --threads=N        The number of threads to use.

Advanced options:
  -f, --fal_pos_rate=N   Maximum false positive rate to use in filter. [0.0075]
  -g, --hash_num=N       Set number of hash functions to use in filter instead
                         of automatically using calculated optimal number of
                         functions.
  -k, --kmer_size=N      K-mer size to use to create filter. [25]
  -s, --subtract=N       Path to filter that you want to uses to prevent the
                         addition of k-mers contained into new filter. You may
                         only use filters with k-mer sizes equal the one you
                         wish to create.
  -n, --num_ele=N        Set the number of expected elements. If set to 0 number
                         is determined from sequences sizes within files. [0]
  -P, --print_reads      During progressive filter creation, print tagged reads
                         to STDOUT in FASTQ format [disabled]
  -r, --progressive=N    Progressive filter creation. The score threshold is
                         specified by N, which may be either a floating point score
                         between 0 and 1 or a positive integer.  If N is a
                         positive integer, it is interpreted as the minimum
                         number of contiguous matching bases required for a
                         match.
  -i, --inclusive        If one paired read matches, both reads will be included
                         in the filter. Only active with the (-r) option.

Report bugs to <cjustin@bcgsc.ca>.
```

### biobloomcategorizer
```
Usage: biobloomcategorizer [OPTION]... -f "[FILTER1]..." [FILE]...
biobloomcategorizer [OPTION]... -e -f "[FILTER1]..." [FILE1.fq] [FILE2.fq]
Categorize Sequences. The input format may be FASTA, FASTQ, qseq, export, SAM or
BAM format and compressed with gz, bz2 or xz and may be tarred.

  -p, --prefix=N         Output prefix to use. Otherwise will output to current
                         directory.
  -f, --filter_files=N   List of filter files to use. Required option. 
                         eg. "filter1.bf filter2.bf"
  -e, --paired_mode      Uses paired-end information. For BAM or SAM files, if
                         they are poorly ordered, the memory usage will be much
                         larger than normal. Sorting by read name may be needed.
  -i, --inclusive        If one paired read matches, both reads will be included
                         in the filter. 
  -s, --score=N          Score threshold for matching. N may be either a
                         floating point score between 0 and 1 or a positive
                         integer representing the minimum match length in bases.
                         If N is a floating point, the maximum threshold is any 
                         number less than 1, and the minimum is 0 (highest
                         sensitivity). If set to 1, the best hit is used rather
                         than the threshold and the score will be appended to the
                         header of the output read. [0.15]
  -w, --with_score       Output multimatches with scores in the order of filter.
  -t, --threads=N        The number of threads to use. [1]
  -g, --gz_output        Outputs all output files in compressed gzip.
      --fa               Output categorized reads in Fasta files.
      --fq               Output categorized reads in Fastq files.
      --chastity         Discard and do not evaluate unchaste reads.
      --no-chastity      Do not discard unchaste reads. [default]
  -l  --length_cutoff=N  Discard reads shorter that the cutoff N. [0]
  -v  --version          Display version information.
  -h, --help             Display this dialog.
Advanced options:
  -m, --min_hit=N        Minimum Hit Threshold Value. The absolute hit number
                         needed over initial tiling of read to continue. Higher
                         values decrease runtime but lower sensitivity.[0]
  -r, --streak=N         The number of hit tiling in second pass needed to jump
                         Several tiles upon a miss. Small values decrease
                         runtime but decrease sensitivity. [3]
  -o, --min_hit_only     Use only initial pass filtering to evaluate reads. Fast
                         but low specificity, use only on long reads (>100bp).
  -c, --ordered          Use ordered filtering. Order of filters matters
                         (filters listed first have higher priority). Only taken
                         advantage of when k-mer sizes and number of hash
                         functions are the same.
  -d, --stdout_filter=N  Outputs all matching reads to stdout for the specified
                         filter. N is the filter ID without file extension.
                         Reads are outputed in fastq, and if paired will output
                         will be interlaced.
```

### A. How can distinguish between organisms that share lots of k-mer content?

Using the `--with_score` (`-w`) in biobloomcategorizer will output the score of each filter (in the order specified with `-f`) in the header of the multimatch filter.

Using this option will make the program run a bit slower but will allow users to see the scores assigned to each filter, so as to make a more informed decision about how the read should be binned.

### B. How can I reduce my memory usage?

Memory usage is directly dependent on the filter size, which is in turn a function of the false positive rate. In biobloommaker reducing memory increases the false positive rate (`-f`) until the memory usage is acceptable. You may need to increase score threshold (`-s`) in biobloomcategorizer to keep the specificity high.

### C. How can I make my results more sensitive?

In biobloomcategorizer try to decrease the score threshold (`-s`). If that still does not work, in biobloommaker try reducing the k-mer (`-k`) size to allow more tiles, which can help with sensitivity.

### D. How can I make my results more specific?

In biobloomcategorizer you can increase score threshold (`-s`).

In biobloommaker decreasing the false positive rate (`-f`) and increasing the k-mer (`-k`) size to allow more tiles can help with specificity. Decreasing the filter false positive rate will increase memory usage.

### E. How can I make the program faster?
There are multiple ways to speed up biobloomcategorizer. Here are a few options:

The `--ordered` option, other than priotizing the first filters in the list (specified by `-f`), will have an added benefit of speeding up the program by avoiding some evaluations if a match is already found. Furthermore, because of this speed up, this option maybe appropriate even in situations where no hierarchy is desired (filters must be unrelated in this case).

In biobloomcategorizer set a min hit threshold (`-m`) >0. This will use a faster rescreening categorization algorithm that uses jumping k-mer tiles to prescreen reads. This will decrease sensitivity but will increase speed. Large values will further decrease sensitivity.

Finally if speed is still an issue, using the min hit threshold only (`-o`) option will use only this screening method and not use the standard sliding tiles algorithm at all. This will greatly increase speed at the expense of sensitivity and specificity. This may be appropriate if your reads are long (>150bp), paired and have minimal read errors. If this method is used, it is recommended that you use an `-m` of at least 2 or 3.


