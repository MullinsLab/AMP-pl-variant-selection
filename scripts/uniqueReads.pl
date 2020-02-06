#!/usr/bin/perl -w

#######################################################################################
# Copyright © 2013 Mullins Lab, University of Washington
# SOFTWARE COPYRIGHT NOTICE 
# This program and its documentation are the copyright of the Mullins Lab, University
# of Washington. All rights are reserved.
# 
# This program is free software. It is supplied without any warranty. You can 
# redistribute it and/or modify it. Mullins Lab is not responsible for its use, 
# misuse, or functionality.
#
# This program was developed by Wenjie Deng <dengw@uw.edu>
#######################################################################################

use strict;
use Getopt::Long;

my %option = (
	'sid'  => '',
	'in'   => '',
	'type' => '',
);

my $usage = "\nusage: perl uniqueReads.pl [-option value]

options:  
-sid    sample ID
-in     input fasta file
-type   type (sequence type: NA or AA)
";

GetOptions (\%option, 'sid=s', 'in=s', 'type=s');

my $sid = $option{'sid'} or die $usage;
my $inFile = $option{'in'} or die $usage;
my $type = $option{'type'} or die $usage;
my $outFile = my $nameFile = my $variantFreqFile = $inFile;
$outFile =~ s/\.fasta/_variants.fasta/;
$nameFile =~ s/\.fasta/_variantList.csv/;
$variantFreqFile =~ s/\.fasta/_variantFreq.csv/;

my $count = my $uniqueCount = 0;
my $seq = my $seqName = '';

my (%seqCount, $uniqDup);
#my %seqLen;
#my %nameTag;
#my (@uniqSeqs, @fwdUniqSeqs, @revUniqSeqs, $uniqDup, $fwdUniqDup, $revUniqDup);
open INFASTA, $inFile or die "couldn't open $inFile: $!\n";
while (my $line = <INFASTA>) {
	chomp $line;
	next if ($line =~ /^\s*$/);
	if ($line =~ /^>(\S+)/) {
		if ($seq) {
			unless ($seq =~ /^\-+$/) {	# in case all gaps in sequence
				if (!$seqCount{$seq}) {
					$seqCount{$seq} = 0;
				}
				push @{$uniqDup->{$seq}}, $seqName;
				++$seqCount{$seq};							
				++$count;					
			}		
			$seqName = $seq = "";
		}
		$seqName = $1;	
	}else {
		$seq .= uc $line;
	}		
}
if ($seq) {
	unless ($seq =~ /^\-+$/) {	# in case all gaps in sequence		
		if (!$seqCount{$seq}) {
			$seqCount{$seq} = 0;
		}
		push @{$uniqDup->{$seq}}, $seqName;
		++$seqCount{$seq};			
		++$count;	
	}
	$seqName = $seq = "";
}
close INFASTA;

open OUT, ">",$outFile or die "couldn't open $outFile: $!\n";
open NAME, ">", $nameFile or die "couldn't open $nameFile: $!\n";
if ($type eq "NA" or $type eq "AA") {
	open FREQ, ">", $variantFreqFile or die "couldn't open $variantFreqFile: $!\n";
	print FREQ "Variant,Count,Frequency\n";
}
foreach my $seq (sort {$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
	$uniqueCount++;
	my $name = '';
	my $frequency = int($seqCount{$seq} / $count * 10000 + 0.5) / 10000;
	if ($uniqueCount =~ /^\d$/) {
		$name = $sid.'_env_'.$type.'_V00'.$uniqueCount."_".$frequency."_".$seqCount{$seq};
	}elsif ($uniqueCount =~ /^\d\d$/) {
		$name = $sid.'_env_'.$type.'_V0'.$uniqueCount."_".$frequency."_".$seqCount{$seq};
	}else {
		$name = $sid.'_env_'.$type.'_V'.$uniqueCount."_".$frequency."_".$seqCount{$seq};
	}
	print OUT ">",$name,"\n",$seq,"\n";
	print NAME $name, ",", join(';', @{$uniqDup->{$seq}}), "\n";
	if ($type eq "NA" or $type eq "AA") {
		my $freq = int($seqCount{$seq} / $count * 10000 + 0.5) / 10000;
		print FREQ "$name,$seqCount{$seq},$freq\n";
	}
}
print "total $count sequences. $uniqueCount variants\n";
close OUT;
close NAME;
if ($type eq "NA" or $type eq "AA") {
	close FREQ;
}
