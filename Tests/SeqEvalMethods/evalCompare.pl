#!/usr/bin/perl -w

#Created Tuesday, May 23 2017
#Author: cjustin
#Assumes dwgsim is in your path

use warnings;
use strict;
use diagnostics;
use IO::File;

#Run Make to compile code
#system("../../configure");
print STDERR `make`;

#@TODO: Download test data (if it doesn't already exist)

#Generate ReadSets for evaluation (if they don't already exist)
unless ( -e "human.bwa.read1.fastq" ) {
	my $readSetCmd = "dwgsim -N 1000000 -1 150 -2 150 -y 0 human.fa human";
	print STDERR `$readSetCmd`;
}
unless ( -e "mouse.bwa.read1.fastq" ) {
	my $readSetCmd = "dwgsim -N 1000000 -1 150 -2 150 -y 0 mouse.fa mouse";
	print STDERR `readSetCmd`;
}
unless ( -e "ecoli.bwa.read1.fastq" ) {
	my $readSetCmd = "dwgsim -N 1000000 -1 150 -2 150 -y 0 ecoli.fa ecoli";
	print STDERR `$readSetCmd`;
}

#Generate base bloom filter (if it doesn't already exist)
unless ( -e "human.bf" ) {
	my $filterGenCmd = "BioBloomMaker/biobloommaker -p human human.fa";
	print STDERR `filterGenCmd`;
}

##Rename current SeqEval.h to SeqEval_temp.h
system("mv ../../Common/SeqEval.h ../../Common/SeqEval_temp.h");

#for each header file found in directory:
my @files = glob('*SeqEval.h');
foreach my $filename (@files) {

	#copy as SeqEval.h
	print STDERR $filename . "\n";
	system( "cp -f " . $filename . " ../../Common/SeqEval.h" );

	#run Make again
	print STDERR `make`;

	#run evaluation with -w option
	my $humanCMD =
"BioBloomCategorizer/biobloomcategorizer -t 8 -s 1 -eiw --fa -p human -f human.bf human.bwa.read1.fastq human.bwa.read2.fastq";
	system($humanCMD);
	my $mouseCMD =
"BioBloomCategorizer/biobloomcategorizer -t 8 -s 1 -eiw --fa -p mouse -f human.bf mouse.bwa.read1.fastq mouse.bwa.read2.fastq";
	system($mouseCMD);
	my $ecoliCMD =
"BioBloomCategorizer/biobloomcategorizer -t 8 -s 1 -eiw --fa -p ecoli -f human.bf ecoli.bwa.read1.fastq ecoli.bwa.read2.fastq";
	system($ecoliCMD);

	#parse and generate summary tsv file
	parseResults( $filename, "human", "human_human_1.fa", "human_human_2.fa" );
	parseResults( $filename, "mouse", "mouse_human_1.fa", "mouse_human_2.fa" );
	parseResults( $filename, "ecoli", "ecoli_human_1.fa", "ecoli_human_2.fa" );
}

sub parseResults {
	my $method    = shift;
	my $type      = shift;
	my $file1     = shift;
	my $file2     = shift;
	my $threshold = 0;
	my $increment = 0.001;
	my $maxCount  = 1;

	my @values1;
	my @values2;

	my $fh1 = IO::File->new( $file1, "r" );
	my $fh2 = IO::File->new( $file2, "r" );
	my $header1 = $fh1->getline();
	my $header2 = $fh2->getline();
	my $seq1    = $fh1->getline();
	my $seq2    = $fh2->getline();

	while ($seq1) {
		my @tmp1 = split( /\s/, $header1 );
		my @tmp2 = split( /\s/, $header2 );

		my $value1 = $tmp1[1];
		my $value2 = $tmp2[1];

		push( @values1, $value1 );
		push( @values2, $value2 );

		$header1 = $fh1->getline();
		$header2 = $fh2->getline();
		$seq1    = $fh1->getline();
		$seq2    = $fh2->getline();
	}
	$fh1->close();
	$fh2->close();

	while ( $maxCount >= $threshold ) {
		my $matchCount = 0;

		for ( my $i = 0 ; $i < scalar(@values1) ; ++$i ) {
			if ( $values1[$i] > $threshold && $values2[$i] > $threshold ) {
				$matchCount++;
			}
		}
		print $method . "\t" . $type . "\t" . $matchCount . "\t" . $threshold . "\n";
		$threshold += $increment;
	}
}

#rename current SeqEval_temp.h back to SeqEval.h
system("mv -f ../../Common/SeqEval_temp.h ../../Common/SeqEval.h");

#clean up
