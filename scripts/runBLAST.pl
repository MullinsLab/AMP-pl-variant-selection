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
use File::Basename;

my %option = (
	'in'     => '',
	'ref'    => '',
	'path'   => '',
	'mt'     => 1,
	'mm'     => -1,
	'go'     => 1,
	'ge'     => 2,
	'h'      => '',
);

my $usage = "\nUsage: perl runBLAST.pl [-option value]

options:
-in      input reads fasta file
-ref     reference fasta file
-path    blast path
-h       usage help

";

GetOptions (\%option, 'in=s', 'ref=s', 'blast=s', 'mt=i', 'mm=i', 'go=i', 'ge=i', 'h');

my $inFasta = $option{'in'} or die $usage;
my $refFasta = $option{'ref'} or die $usage;
my $reward = $option{'mt'};
my $penalty = $option{'mm'};
my $gapopen = $option{'go'};
my $gapextend = $option{'ge'};
my $help = $option{'h'};
die $usage if ($help);
my $blastPath = $option{'blast'};
my $outXML = $inFasta;
if ($inFasta =~ /\.mafft\.fa/) {
	$outXML =~ s/\.mafft\.fa//;
}
$outXML =~ s/\.fasta/.xml/;

print "Making BLAST db ... ";
my $nhrfile = $refFasta.".nhr";
my $ninfile = $refFasta.".nin";
my $nsqfile = $refFasta.".nsq";
my $rv = 0;
if (-e $nhrfile and -e $ninfile and -e $nsqfile) {
	print "blast db already exists.\n";
}else {
	my $formatdbLog = $refFasta . ".log";
	$rv = system("$blastPath/makeblastdb -in $refFasta -dbtype nucl -logfile $formatdbLog");
	unless ($rv == 0) {
		die "\nmakeblastdb failed: $rv\n";
	}
	print "done.\n";
}

print "BLASTing ... ";
$rv = system ("$blastPath/blastn -task blastn -db $refFasta -query $inFasta -out $outXML -reward $reward -penalty $penalty -gapopen $gapopen -gapextend $gapextend -dust no -outfmt 5");
unless ($rv == 0) {
	die "\nBLAST failed: $rv\n";
}
print "done.\n";
