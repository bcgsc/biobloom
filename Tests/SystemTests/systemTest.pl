#!/usr/bin/perl -w

#Created Tuesday, May 23 2017
#Author: cjustin
#Assumes dwgsim is in your path

use warnings;
use strict;
use diagnostics;
use IO::File;

#TODO make hardcoded values into options
#my $timeCmd = "/usr/bin/time -v ";
my $timeCmd = "/usr/bin/time -v ";
my $srcDir  = "/home/cjustin/git/biobloom/Tests/SystemTests/bin/";
#my $srcDir  = "~/arch/centos6/bin/";
my $step    = 0;
if ( $ARGV[0] ) {
	$step = $ARGV[0];
}
my @files = glob('*.fa');    #put some fasta files in the path

#TODO add more automatic checks rather than manual evaluation of results

if ( $step == 0 ) {

	#for each header file found in directory:
	foreach my $filename (@files) {

		my $filterPrefix = $filename;
		$filterPrefix =~ s/\.fa$//;
		mkdir($filterPrefix);
		chdir($filterPrefix);

		#Generate ReadSets for evaluation (if they don't already exist)
		unless ( -e "$filterPrefix.bwa.read1.fastq" ) {
			my $readSetCmd =
			  "dwgsim -N 1000000 -1 150 -2 150 -y 0 ../$filename $filterPrefix";
			print STDERR `$readSetCmd`;
		}

		#Generate base bloom filters
		my $filterGenCmd =
		    $timeCmd
		  . $srcDir
		  . "biobloommaker -t 8 -p $filterPrefix ../$filename";
		print STDERR `$filterGenCmd`;

		#simplest test
		my $testCMD1 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -t 8 -p $filterPrefix -f $filterPrefix.bf $filterPrefix.bwa.read1.fastq";
		print STDERR `$testCMD1`;

		#test paired reads
		my $testCMD2 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -t 8 -p $filterPrefix -e -f $filterPrefix.bf $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq";
		print STDERR `$testCMD2`;

		#test integer scores
		my $testCMD3 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -t 8 -s 5 -p $filterPrefix -e -f $filterPrefix.bf $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq";
		print STDERR `$testCMD3`;

		#test best hit option
		#test paired reads
		my $testCMD4 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -t 8 -s 1 -p $filterPrefix -e -f $filterPrefix.bf $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq";
		print STDERR `$testCMD4`;

		#inclusive mode
		#test paired reads
		my $testCMD5 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -t 8 -s 1 -i -p $filterPrefix -e -f $filterPrefix.bf $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq";
		print STDERR `$testCMD5`;

		my $testCMD6 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -t 8 --fa -g -s 1 -i -p $filterPrefix -e -f $filterPrefix.bf $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq";
		print STDERR `$testCMD6`;

		#test print filter
		my $testCMD7 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -d -t 8 --fa -g -s 1 -i -p $filterPrefix -e -f $filterPrefix.bf $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq > temp.fq";
		print STDERR `$testCMD7`;

		#ordered mode
		my $testCMD8 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -c -d -t 8 -g -i -p $filterPrefix -e -f $filterPrefix.bf $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq > temp.fq";
		print STDERR `$testCMD8`;

		system("rm *.fa");

		chdir("../");
	}
}

#multibloom filter tests
if ( $step <= 1 ) {
	my $filterSets = "";
	mkdir("multibloom");
	chdir("multibloom");
	foreach my $file (@files) {
		my $filterPrefix = $file;
		$filterPrefix =~ s/\.fa$//;
		$filterPrefix = "../" . $filterPrefix . "/" . $filterPrefix;

		unless ( -e "$filterPrefix.bwa.read1.fastq" ) {
			my $readSetCmd =
			  "dwgsim -N 1000000 -1 150 -2 150 -y 0 ../$file $filterPrefix";
			print STDERR `$readSetCmd`;
		}

		unless ( -e "$filterPrefix.bf" ) {

			#Generate base bloom filters
			my $filterGenCmd =
			    $timeCmd
			  . $srcDir
			  . "biobloommaker -t 8 -p $filterPrefix ../$file";
			print STDERR `$filterGenCmd`;
		}

		$filterSets .= $filterPrefix . ".bf ";
	}

	foreach my $file (@files) {
		my $filterPrefix = $file;
		$filterPrefix =~ s/\.fa$//;
		$filterPrefix = "../" . $filterPrefix . "/" . $filterPrefix;
		my $testCMD1 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -t 8 -p multibloom "
		  . " -e -f \"$filterSets\" $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq";
		print STDERR `$testCMD1`;

		my $testCMD2 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -t 8 -p multibloom "
		  . " -e -s 1 -f \"$filterSets\" $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq";
		print STDERR `$testCMD2`;

		my $testCMD3 =
		    $timeCmd
		  . $srcDir
		  . "biobloomcategorizer -t 8 -p multibloom "
		  . " -e -c -f \"$filterSets\" $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq";
		print STDERR `$testCMD3`;
	}
	chdir("../");
}

#progressive bloomfilter tests @TODO
#multibloom filter tests
if ( $step <= 2 ) {

	#for each header file found in directory:
	foreach my $filename (@files) {

		my $filterPrefix = $filename;
		$filterPrefix =~ s/\.fa$//;
		mkdir($filterPrefix);
		chdir($filterPrefix);

		#Generate ReadSets for evaluation (if they don't already exist)
		unless ( -e "$filterPrefix.bwa.read1.fastq" ) {
			my $readSetCmd =
			  "dwgsim -N 1000000 -1 150 -2 150 -y 0 ../$filename $filterPrefix";
			print STDERR `$readSetCmd`;
		}

		#Generate base bloom filters
		my $filterGenCmd =
		    $timeCmd . "bash -c '"
		  . $srcDir
		  . "biobloommaker -r 0.15 -i -t 8 -n 50000 -p "
		  . $filterPrefix
		  . "_progressive <(seqtk trimfq -L 1000 ../$filename) $filterPrefix.bwa.read1.fastq $filterPrefix.bwa.read2.fastq'"
		  ;
		print STDERR `$filterGenCmd`;
		chdir("../");
	}
}
