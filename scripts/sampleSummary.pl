#!/usr/bin/perl -w

#########################################################################################################
# Program: sampleSummary.pl
# Purpose: read sample's _count_summary.csv, write to a summary file for all samples
# Output: .csv file
# Author: Wenjie Deng
# Date: 2019-04-15
# Modified: 2020-02-05
###########################################################################################################

use strict;

my $usage = "perl sampleSummary.pl inputCountSummayfile outSampleSummaryfile\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $indir = shift || '.'; # default current directory
my $samplecount = 0;

open IN, $infile or die "couldn't open $infile: $!\n";
open OUT, ">>", $outfile or die "couldn't open $outfile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^Sample ID,Consensus/) {
		if (-z $outfile) {
			print OUT "$line\n";
		}		
	}else {
		print OUT "$line\n";
	}
}
close IN;
close OUT;

