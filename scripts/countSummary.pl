#!/usr/bin/perl -w

#########################################################################################################
# Program: countSummary.pl
# Input: log file and AA variant fasta file
# Output: .csv file
# Author: Wenjie Deng
# Date: 2019-04-10
###########################################################################################################

use strict;

my $usage = "perl countSummary.pl sampleid input_log_file input_aa_variant_fasta_file output_count_summary_file\n";
my $sid = shift or die $usage;
my $logfile = shift or die $usage;
my $inFastaFile = shift or die $usage;
my $outfile = shift or die $usage;

my $totalSeq = my $variantcount = my $flag = my $cleancount = 0;

open LOG, $logfile or die "couldn't open $logfile: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
print OUT "Sample ID,Consensus sequences,Incomplete,Large deletion,Frame shift,Internal stop codon,Intact gp160,Fraction of intact gp160,NT_variants,AA_variants,Variants > 9.5%\n";
while (my $line = <LOG>) {
	chomp $line;
	if ($line =~ /Total (\d+) sequences in input fasta file/) {
		my $seqcount = $1;
		print OUT "Consensus sequences,$seqcount\n";
		if ($logfile =~ /_pe_log.txt/) {
			print OUT "After QC,$seqcount\n";
		}
	}
	if ($line =~ /\d+ sequences' family size fewer than (\d+) or min agreement fewer than (.*?), (\d+) clean sequences/) {
		print OUT "fs >= $1 and min agreement >= $2,$3\n";
	}elsif ($line =~ /total sequences: (\d+)/) {
		$cleancount = $1;
		print OUT "$sid,$cleancount";
	}elsif ($line =~ /hits cover the region of \d+ \- \d+: \d+, large deletion: (\d+), in frame: \d+, out frame: (\d+), not cover: (\d+)/) {
		print OUT ",$3,$1,$2";
		$cleancount -= $3;
	}elsif ($line =~ /complete AA sequences: (\d+), defective AA sequence: (\d+)/) {
		$totalSeq = $1;
		my $fraction = int($totalSeq / $cleancount * 100 + 0.5) / 100;
		print OUT ",$2,$totalSeq,$fraction";
	}elsif ($line =~ /total \d+ sequences. (\d+) variants/) {
		if (!$flag) {
			print OUT ",$1";
		}
		if ($flag == 1) {
			print OUT ",$1";
			last;
		}
		++$flag;
	}
}
close LOG;

open INFASTA, $inFastaFile or die "couldn't open $inFastaFile: $!\n";
while (my $line = <INFASTA>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(.*)_(\d+)$/) {
		my $count = $2;
		my $fraction = $count / $totalSeq;
		if ($fraction > 0.095) {
			++$variantcount;
		}else {
			last;
		}
	}
}
print OUT ",$variantcount\n";
close INFASTA;
close OUT;
