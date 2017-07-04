#!/usr/bin/perl -w
#Written by Justin Chu (cjustin@bcgsc.ca)

use warnings;
use strict;
use diagnostics;
use IO::File;
use Getopt::Std;
use vars qw($opt_f $opt_s $opt_1 $opt_2);

#f - filter name list, s - score threshold, 1 - file1, 2 - file2
getopts('f:s:1:2:');

my $score = 0.15;
if ($opt_s) {
	$score = $opt_s;
}

my @filterNames = split( /\s/, $opt_f );
my %summarizedCounts;

for ( my $i = 0 ; $i < scalar(@filterNames) ; $i++ ) {
	for ( my $j = $i + 1 ; $j < scalar(@filterNames) ; $j++ ) {
		$summarizedCounts{$i}{$j} = 0;
	}
}

#read in read pairs
my $fh1   = new IO::File( $opt_1, "r" );
my $line1 = $fh1->getline();
my $fh2   = new IO::File( $opt_2, "r" );
my $line2 = $fh2->getline();

while ($line1) {
	chomp($line1);
	if ( $line1 =~ /^@/ ) {
		my @tempArr1 = split( /\s/, $line1 );
		my @tempArr2 = split( /\s/, $line2 );
		for ( my $i = 1 ; $i < scalar(@tempArr1) ; $i++ ) {
			for ( my $j = $i + 1 ; $j < scalar(@tempArr1) ; $j++ ) {
				if (   $tempArr1[$i] > $score
					&& $tempArr2[$i] > $score
					&& $tempArr1[$j] > $score
					&& $tempArr2[$j] > $score )
				{
					#index offset by 1 due to read name in the front
					$summarizedCounts{ $i - 1 }{ $j - 1 }++;
				}
			}
		}
	}
	$line1 = $fh1->getline();
	$line2 = $fh2->getline();
}
$fh1->close();
$fh2->close();

for ( my $i = 0 ; $i < scalar(@filterNames) ; $i++ ) {
	for ( my $j = $i + 1 ; $j < scalar(@filterNames) ; $j++ ) {
		print $filterNames[$i] . "\t"
		  . $filterNames[$j] . "\t"
		  . $summarizedCounts{$i}{$j} . "\n";
	}
}

