# BioBloom Tools User Manual - Multi-Index Bloom Filter sequence classifier

Multi-Index Bloom Filters (miBF) allow for efficent classification of a many targets. The use of miBFs is recommend if classifying against >100 targets as `biobloomcategorizer` cannot scale to this many targets. It also also for higher sensitivity through the use of spaced seeds.

## A preprint is available [here](https://www.biorxiv.org/content/biorxiv/early/2018/10/05/434795.full.pdf):
Data and filter used for Read Binning experiments are found [here](http://www.bcgsc.ca/downloads/btl/bbt/Metagenomic_miBF)
Data and filter used for Metagenomic classification experiments is found [here](http://www.bcgsc.ca/downloads/btl/bbt/CancerCensusGenes_miBF/)

Table of Contents
------
1. [Generating Multi-Index Bloom Filters from Reference Sequences with biobloommimaker](#1)
2. [Classifying and Analyzing Sequences with biobloommicategorizer](#2)
3. [Output files](#3)
4. [Multiple spaced spaced seeds](#4)
5. [Discussion and FAQs](#5)

<a name="1"></a>
## 1. Generating Multi-Index Bloom Filters from Reference Sequences with biobloommimaker

Before you can index your fasta files a reverse complement version of them must first be generated with the `util/genRV.sh` script:
```bash
util/genRV.sh input1.fa input2.fa 
```
it will output `<fasta>.rv` files for each fasta file listed.

A `<prefix>.bf`, `<prefix>.bf.sdsl` and `<prefix>_ids.txt` can be created by running the command:
```bash
./biobloommimaker -p filterName -S "01111011100000011111101110000111100100000101001001011100011000010111101001100011 10101101001111010011100101110001111111111010010000110011001000100100011100001000 00010000111000100100010011001100001001011111111110001110100111001011110010110101 11000110010111101000011000111010010010100000100111100001110111111000000111011110" input1.fasta ...
```
`-p` is the prefix for the files being created.

The options above are the bare minimum options you must use to run the program, but it is possible to customize many aspects of your filter that can drastically change performance depending on your needs. You can also use the `--help` parameter for a listing of the options. Indeed, the multiple space seeds used in this example may not be optimal for your purposes. By default, each header in the file will be indexed as a seperate entry.

### Help dialog of biobloommimaker:

```
Usage: biobloommimaker -p [FILTERID] [OPTION]... [FILE]...
Usage: biobloommimaker -p [FILTERID] -r 0.2 [FILE]... [FASTQ1] [FASTQ2] 
Creates a multi-index Bloom Filter from a list of fasta files.

  -p, --file_prefix=N    Filter prefix and filter ID. Required option.
  -h, --help             Display this dialog.
      --version          Display version information.
  -v  --verbose          Display verbose output.
  -t, --threads=N        The number of threads to use.
  -b, --occupancy=N      Occupancy of Bloom filter.[0.5]
  -n, --num_ele=N        Set the number of expected distinct k-mer frames.
                         If set to 0 number is determined from sequences
                         sizes within files. [0]
  -S, --seed_str=N       Generate a miBF using multiple spaced seeds. Expects
                         list of seed 1s & 0s separated by spaces.
  -F, --by_file          Assign IDs by file rather than by fasta header
k-mer mode options (disabled when using spaced seeds):
  -g, --hash_num=N       Set number of hash functions when using k-mers.
  -k, --kmer_size=N      K-mer size to use to create filter. [25]
```

### Notes:
It is suggested that at least 3 spaced seeds be used, and ideally, the references used to be similar in size to each other.

<a name="3"></a>
## 2. Classifying and Analyzing Sequences with biobloommicategorizer

Once you have the miBF created, you can use it with `biobloommicategorizer` to categorize sequences. Currently only fasta and fastq file are support (can be paired end). Native gzip file reading is support but for optimal performance it is recommended for the use to use `pigz -cd`.

Classification of reads in a tsv format (no sequences, only read names) to `stdout` can be generated with this command: 
```bash
./biobloommicategorizer -f filterName.bf inputReads.fa
```
`-f` is the path to the `prefix.bf` file generated in the previous step.

After execution, a `<prefix>_summary.tsv` (by default the prefix is empty) will be generated similar to normal BBT.

The options above are the bare minimum options you must use to run the program, but it likely that you will need to add additional options for your needs. For example, a fasta output (with the classification in the read header) can be outputed with `--fa` rather that outputting a tsv file. You can also use the `--help` parameter for a listing of the options. Indeed, the multiple space seeds used in this example may not be optimal for your purposes. By default, each header in the file will be indexed as a seperate entry.

### Help dialog of biobloommicategorizer:

```
Usage: biobloommicategorizer [OPTION]... -f [FILTER] [FILE]...
Usage: biobloommicategorizer [OPTION]... -f [FILTER1] -e [FILE1.fq] [FILE2.fq]
The input format may be FASTA, FASTQ, and compressed with gz. Output will be to
stdout with a summary file outputted to a prefix.
  -p, --prefix=N         Output prefix to use. Otherwise will output to current
                         directory.
  -f, --filter=N         Path of miBF '.bf' file
  -e, --paired_mode      Uses paired-end information. For BAM or SAM files, if
                         they are poorly ordered, the memory usage will be much
                         larger than normal. Sorting by read name may be needed.
  -s, --min_FPR=N        Minimum -10*log(FPR) threshold for a match. [100.0]
  -t, --threads=N        The number of threads to use. [1]
      --fa               Output categorized reads in Fasta output.
      --fq               Output categorized reads in Fastq output.
      --version          Display version information.
  -h, --help             Display this dialog.
  -v, --verbose          Display verbose output
  -I, --interval         the interval to report file processing status [10000000]
Advanced options:
  -a, --frameMatches=N   Min number seed matches needed in a frame to match [1]
                         Ignored if k-mers used when indexing.
  -m, --multi=N          Multi Match threshold.[2]
  -r, --streak=N         Number of additional hits needed to skip classification. [10]
  -c, --minNoSat         Minimum count of non saturated matches. Increasing this value
                         will filter out repetitive or low quality sequences.[0]
  -b, --bestHitAgree     Filters out all matches where best hit is ambiguous because
                         match count metrics do not agree.
      --debug            debug filter output mode.
      
Report bugs to <cjustin@bcgsc.ca>.
```
### Notes:
We suggested that sequences being classified be of fixed length in the same run. It should not be problem if the a few different sizes (e.g. 2 different sizes in barcode trimmed Chromium data) are but we do not recommend a set of sequences with all different lengths. This is due to limitations on our current implementation relating to the false postive rate calcuations (please email cjustin@bcgsc.ca for more details). It may be possible to optimize the code for variable lengths.

<a name="3"></a>
## 3. Output files
**biobloommimaker**

 * `<prefix>.sdsl`: This is the raw Bloom filter with rank information
   interleaved in. The sdsl extension come from the sdsl library used.
 * `<prefix>.bf`: Contains filter information in the header (size, k-mer
   size, spaced seeds, used etc.) followed by a vector of index IDs.
 * `<prefix>_ids.txt`: Plain text file containing the ID names
   used in the filter. The names are usually extracted from the
   reference fasta file headers or the filenames themselves. It is
   possible to alter them here if the default naming does not match the
   user's desired output names.

These files need to have the same prefix and be in the same directory for the filter to function.

**biobloommicategorizer**
* `<prefix>_summary.tsv`: A summary file of the counts of hits to each filter. Similar to the summary file provided by `biobloomcategorier`.
* STDOUT output: Seqeunce classificiation results are outputted to stdout
	* `--tsv`: The default output. Tab seperated file that contains the seqeunce name and classfication results in addition to the count metrics used to classify the seqeunce.
	* `--fa`:  FastA format output. Classification results are part of the header.
	* `--fq`:  FastQ format output. Classification results are part of the header.

<a name="4"></a>
## 4. Multiple Spaced seeds
The design of spaced seed use in BBT must be either palindromic or have another complementry seed. The design of seed depends on the degree of error tollerance needed, and the specificity of classification needed. For example, if similar sequences are indexed, then longer seed are recommended. If classifying to a divergent or highly erroreous sequence lower weight seeds are recommended.

Here is an example of some of randomly generates seed sets that can be used in our tool:
```
#seed 1: least specific, 3 seeds
001100000010111000010100000000110011000101011000100010100001 010010101100000001000011001100001100110000100000001101010010 100001010001000110101000110011000000001010000111010000001100

#seed 2: Shorter seeds, 4 seeds (more memory usage)
10101101011000001000010000110111100111110101101010 11100000000111010001110111010001111100001111011000 00011011110000111110001011101110001011100000000111 01010110101111100111101100001000010000011010110101

#seed 3: Long seeds (more specific), 4 seeds
00000101100001000001000000101000101000010000001010000010101100100001000001010000 00010000001000010010001010010110000101100001000000000000000000010100101000000111 11100000010100101000000000000000000010000110100001101001010001001000010000001000 00001010000010000100110101000001010000001000010100010100000010000010000110100000

#seed 4: Highly dense seeds (very specific, not as error tollerant) 3 seeds
11011110111011101001111011101111101110110010000110 10111111110110011110001100110001111001101111111101 01100001001101110111110111011110010111011101111011
```
The downsides to using longer seeds entail being unable to handle indels as well limitations on classifying shorter sequences (seeds must be smaller than the sequence).

Optimizing spaced seeds in an open area of research. However, in practice we have found that even randomly generated spaced seeds tend to have better sensitivity than k-mers. To generate some compatible random spaced seeds for experimentation they script `util/designSS.py` is worth looking into.

<a name="5"></a>
## 5. Discussion and FAQs
Have a question? Feel free to post an issue or email cjustin@bcgsc.ca.

### Do you output multimapping sequences?
We provided multimaps if a sequence is sufficently similar two multiple indexed references. This can be tweaked by changing the `-m` parameter, which is the difference between the number of hits in the best and secondary hit counts before we consider the match too close to tell apart. Thus, increasing it will allow more secondary hits to be outputted.

###  What if my genomes are very similar?
For every entry in the same spaced seed, we only allow for one reference sequence to take ownership of it. Thus, if a sequences is shared between multiple classes it is randomly assigned to be owned by one the references. In a similar vein, random collisions can occur which will also lead to one of the references to take ownership. We have measures in place to tollerate these events by attempting to fairly distributing collisions and marking entries as "saturated" if they occur in the same locus, which works well against random collisions. However, it is suggested to minimize using very similar reference sequence in the same miBF or to group them together. Another option is to increase the specificity of the seeds used (use longer, denser spaced seeds).

### What are the memory requirements?
The memory required to use this datastructure is dependent on the number of spaced seeds used and the size of the reference used in the file. Other lesser factor include the occupancy of the Bloom filter. Compared to normal BBT expect memory usage to be 3-4 times higher. The number of bits per basepair with default settings and 3 spaced seeds is around 5 bytes.

### Is there a limit to number of sequences I can index a time?
You can currently index a maximum of 32768 sequences at a time.

### Anything else I should know?
The more sequences you index, the better the specifity and sensivity of the method because adding more sequences helps to reduce the FPR of other sequences in the filter. For example, if filtering for viruses in a patient sample, it may be helpful to add human sequences as they will reduce the FPR of the viral sequences in the filter, thus requiring fewer seed hits to be considered a match (higher sensitivity).
