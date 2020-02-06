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
use File::Basename;
use File::Path;
use Getopt::Long;

my %option = (
	'fas' => '',
	'xml' => '',
	'ref' => '',
	'rs'  => 1,
	're'  => 0,
	'dlx' => '',
	'h'   => '',
);

my $usage = "\nusage: perl retrieveRegion.pl [-option value]

options:
-fas    input blasted sequence fasta file
-xml    input blast output in xml format
-ref    input reference sequence used for blast in fasta format
-rs     region start position in reference sequence (default: 1, beginning of reference)
-re     region end position in reference sequence (default: 0, end of reference)
-dlx    flag for deleting xml file after running the script (default: false)
-h      usage help

";

GetOptions (\%option, 'fas=s', 'xml=s', 'ref=s', 'rs=i', 're=i', 'dlx', 'h');
my $fas = $option{'fas'} or die $usage;
my $xml = $option{'xml'} or die $usage;
my $refFile = $option{'ref'} or die $usage;
my $startPos = $option{'rs'};
my $endPos = $option{'re'};
my $dlx = $option{'dlx'};
my $help = $option{'h'};
die $usage if $help;

my $uncoverfile = $xml;
$uncoverfile =~ s/_REN_(\w+)\.xml/_env_$1.xml/;
$uncoverfile =~ s/_GP_(\w+)\.xml/_gag_$1.xml/;
$uncoverfile =~ s/\.xml/_NT_end_missing.fasta/;
my $ldeletionfile = $xml;
$ldeletionfile =~ s/_REN_(\w+)\.xml/_env_$1.xml/;
$ldeletionfile =~ s/_GP_(\w+)\.xml/_gag_$1.xml/;
$ldeletionfile =~ s/\.xml/_NT_largeDeletion.fasta/;
my $outframefile = $xml;
$outframefile =~ s/_REN_(\w+)\.xml/_env_$1.xml/;
$outframefile =~ s/_GP_(\w+)\.xml/_gag_$1.xml/;
$outframefile =~ s/\.xml/_NT_outframe.fasta/;
my $inframefile = $xml;
$inframefile =~ s/_REN_(\w+)\.xml/_env_$1.xml/;
$inframefile =~ s/_GP_(\w+)\.xml/_gag_$1.xml/;
$inframefile =~ s/\.xml/_NT_inframe.fasta/;
my $ntdefectivefile = $xml;
$ntdefectivefile =~ s/_REN_(\w+)\.xml/_env_$1.xml/;
$ntdefectivefile =~ s/_GP_(\w+)\.xml/_gag_$1.xml/;
$ntdefectivefile =~ s/\.xml/_NT_stop.fasta/;
my $nafile = $xml;
$nafile =~ s/_REN_(\w+)\.xml/_env_$1.xml/;
$nafile =~ s/_GP_(\w+)\.xml/_gag_$1.xml/;
$nafile =~ s/\.xml/_NT.fasta/;
my $aafile = $xml;
$aafile =~ s/_REN_(\w+)\.xml/_env_$1.xml/;
$aafile =~ s/_GP_(\w+)\.xml/_gag_$1.xml/;
$aafile =~ s/\.xml/_AA.fasta/;
my $aadefectivefile = $xml;
$aadefectivefile =~ s/_REN_(\w+)\.xml/_env_$1.xml/;
$aadefectivefile =~ s/_GP_(\w+)\.xml/_gag_$1.xml/;
$aadefectivefile =~ s/\.xml/_AA_stop.fasta/;

my ($readName, $seq, %readAlignStart, %readAlignEnd, $referenceSeq, $readSeq, %alignLen, %refAlignNaLen, @alignReads);
my %mismatchCount = my %nameFrame = my %nameSeq = ();
my $fileName = basename($xml);
$fileName =~ s/\.(.*?)$//;
my $queryCount = my $queryLen = my $hitCount = my $hitFlag = my $hitStart = my $hitEnd = 0;
my $qStart = my $qEnd = my $regionFlag = my $rCutoffHitCount = my $notPassRcutHitCount = my $frame = 0;
my $inframecount = my $outframecount = my $completecount = my $defectivecount = my $ldeletioncount = 0;
my $targetqstart = my $targetqend = my $seqcount = my $startflag = my $endflag = 0;
my $refSeq = my $alignedReadSeq = my $alignedRefSeq = '';
my ($refRegionSeq, $refRegionStatus, $readRegionSeq, $readRegionStatus, $readRegionDup);
open REF, $refFile or die "couldn't open $refFile: $!\n";
while (my $line = <REF>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$refSeq .= $line;
	}
}
close REF;
$refSeq =~ s/\-//g;
my $refLen = length $refSeq;
if ($endPos == 0) {
	$endPos = $refLen;
}
print "refLen: $refLen\n";
print "start $startPos, end $endPos\n";

my $seqname = '';
open FAS, $fas or die "couldn't open $fas: $!\n";
while (my $line = <FAS>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(.*)/) {
		$seqname = $1;
		++$seqcount;
	}else {
		$line =~ s/\-//g;
		$nameSeq{$seqname} .= uc $line;
	}
}

open XML, $xml or die "couldn't open $xml: $!\n";
open UNCOVER, ">", $uncoverfile or die "couldn't open $uncoverfile: $!\n";
open INFRAME, ">", $inframefile or die "couldn't open $inframefile: $!\n";
open OUTFRAME, ">", $outframefile or die "couldn't open $outframefile: $!\n";
open DELETION, ">", $ldeletionfile or die "couldn't open $ldeletionfile: $!\n";
open NA, ">", $nafile or die "couldn't open $nafile: $!\n";
open AA, ">", $aafile or die "couldn't open $aafile: $!\n";
open AADEFECTIVE, ">", $aadefectivefile or die "couldn't open $aadefectivefile: $!\n";
open NTDEFECTIVE, ">", $ntdefectivefile or die "couldn't open $ntdefectivefile: $!\n";
while (my $line = <XML>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /<Iteration_query-def>(.*)<\/Iteration_query-def>/) {
		$readName = $1;
		$queryCount++;
		$hitFlag = $targetqstart = $targetqend = $startflag = $endflag = 0;
	}elsif ($line =~ /<Iteration_query-len>(\d+)<\/Iteration_query-len>/) {
		$queryLen = $1;
	}elsif ($line =~ /<Hit_num>1<\/Hit_num>/) {	# sometimes there are more than one hit, we just get the first one.
		$hitFlag = 1;
		$hitCount++;
	}elsif ($hitFlag and $line =~ /<Hsp_query-from>(\d+)<\/Hsp_query-from>/) {
		$qStart = $1 - 1;
	}elsif ($hitFlag and $line =~ /<Hsp_query-to>(\d+)<\/Hsp_query-to>/) {
		$qEnd = $1;
	}elsif ($hitFlag and $line =~ /<Hsp_hit-from>(\d+)<\/Hsp_hit-from>/) {
		$hitStart = $1;
	}elsif ($hitFlag and $line =~ /<Hsp_hit-to>(\d+)<\/Hsp_hit-to>/) {
		$hitEnd = $1;
		unless ($hitEnd > $hitStart) {
			my $tmp = $hitStart;
			$hitStart = $hitEnd;
			$hitEnd = $tmp;
		}
	}elsif ($hitFlag and $line =~ /<Hsp_hit-frame>(.*)<\/Hsp_hit-frame>/) {
		$frame = $1;
	}elsif ($hitFlag and $line =~ /<Hsp_qseq>(.*)<\/Hsp_qseq>/) {
		$alignedReadSeq = uc $1;
	}elsif ($hitFlag and $line =~ /<Hsp_hseq>(.*)<\/Hsp_hseq>/) {
		$alignedRefSeq = uc $1;
		if ($hitStart <= $startPos and !$startflag) {
			if (!defined $nameFrame{$readName}) {
				$nameFrame{$readName} = $frame;
			}elsif ($nameFrame{$readName} != $frame) {
				die "oppsite orientation: $readName\n";
			}
			$startflag = 1;
			if ($frame == -1) {	# reverse complement
				$alignedReadSeq = ReverseComplement ($alignedReadSeq);
				$alignedRefSeq = ReverseComplement ($alignedRefSeq);
			}
			my $alignedRefLen = length $alignedRefSeq;
			my @alignedRefNas = split //, $alignedRefSeq;
			my @alignedReadNas = split //, $alignedReadSeq;
			my $idx = 0;
			if ($frame == -1) {
				$targetqstart = $queryLen - $qEnd;
			}else {
				$targetqstart = $qStart;
			}
			for (my $i = 0; $i < $alignedRefLen; $i++) {
				if ($alignedReadNas[$i] =~ /[A-Za-z]/) {
					++$targetqstart;
				}
				if ($alignedRefNas[$i] =~ /[A-Za-z]/) {
					$idx++;
					if ($idx == $startPos-$hitStart+1) {
						last;
					}
				}
			}
		}

		if ($hitEnd >= $endPos and !$endflag) {
			if (!defined $nameFrame{$readName}) {
				$nameFrame{$readName} = $frame;
			}elsif ($nameFrame{$readName} != $frame) {
				die "oppsite orientation: $readName\n";
			}
			$endflag = 1;
			if ($frame == -1) {	# reverse complement
				$alignedReadSeq = ReverseComplement ($alignedReadSeq);
				$alignedRefSeq = ReverseComplement ($alignedRefSeq);
			}
			my $alignedRefLen = length $alignedRefSeq;
			my @alignedRefNas = split //, $alignedRefSeq;
			my @alignedReadNas = split //, $alignedReadSeq;
			my $idx = 0;
			if ($frame == -1) {
				$targetqend = $queryLen - $qEnd;
			}else {
				$targetqend = $qStart;
			}
			for (my $i = 0; $i < $alignedRefLen; $i++) {
				if ($alignedReadNas[$i] =~ /[A-Za-z]/) {
					++$targetqend;
				}
				if ($alignedRefNas[$i] =~ /[A-Za-z]/) {
					$idx++;
					if ($idx == $endPos-$hitStart+1) {
						last;
					}
				}
			}
		}
	}elsif ($hitFlag and $line =~ /<\/Hit_hsps>/) {
		if ($startflag and $endflag) {
			$rCutoffHitCount++;
			my $seq = $nameSeq{$readName};
			if ($nameFrame{$readName} == -1) {
				$seq = ReverseComplement($seq);
			}
			my $retrievedRead = substr($seq, $targetqstart-1, $targetqend-$targetqstart+1);
			if ($retrievedRead =~ /AAA$/) {
				$retrievedRead =~ s/A$//;
			}
			my $retrievedReadLen = length $retrievedRead;
			my $retrievedRefLen = $endPos - $startPos + 1;
			my $envname = $readName;
			$envname =~ s/_REN_/_env_/;
			if ($retrievedReadLen / $retrievedRefLen < 0.8) {
				++$ldeletioncount;
				print DELETION ">$envname\n$retrievedRead\n";
			}else {
				my @retrievedReads = ();
				my $ambflag = 0;
				if ($retrievedRead =~ /^[ACGT]+$/) { # no ambiguity
					push @retrievedReads, $retrievedRead;
				}else {	# ambiguity
					$ambflag = 1;
					my @ambidxes = my @unambiguous = ();
					my @nas = split //, $retrievedRead;
					for (my $i = 0; $i < $retrievedReadLen; $i++) {
						unless ($nas[$i] =~ /[ACGT]/) {
							push @unambiguous, deambiguous($nas[$i]);
							push @ambidxes, $i;
						}
					}
					my @finals = ();
					my $flag = 0;
					foreach my $node (@unambiguous) {
						++$flag;
						my @splits = split //, $node;
						if ($flag == 1) {
							@finals = @splits;
						}else {
							my @temps = ();
							foreach my $final (@finals) {
								foreach my $split (@splits) {
									my $comb = $final.$split;
									push @temps, $comb;
								}
							}
							@finals = @temps;
						}
					}
					foreach my $final (@finals) {
						my @splits = split //, $final;
						my $splitsidx = 0;
						foreach my $ambidx (@ambidxes) {
							$nas[$ambidx] = $splits[$splitsidx];
							++$splitsidx;
						}
						my $splitseq = join('', @nas);
						push @retrievedReads, $splitseq;
					}
				}
				my $combinationidx = 0;
				foreach my $retrievedRead ( @retrievedReads) {
					my $name = $envname;
					if ($ambflag) {
						++$combinationidx;
						$name .= '_C'.$combinationidx;
					}
					if ($retrievedReadLen % 3 == 0) { # in frame
						++$inframecount;
						print INFRAME ">$name\n$retrievedRead\n";
						my @nas = split //, $retrievedRead;
						my $aaSeq = "";
						my $dflag = 0;
						for (my $i = 0; $i < $retrievedReadLen; $i += 3) {
							my $codon = $nas[$i].$nas[$i+1].$nas[$i+2];
							my $aa = translation($codon);
							if ($aa eq '*' and ($i != $retrievedReadLen-3)) {
								$dflag = 1;
							}
							$aaSeq .= $aa;
						}
						if ($dflag == 0) {
							++$completecount;
							if ($aaSeq =~ /\*$/) {
								$aaSeq =~ s/\*$//;
							}
							print AA ">$name\n$aaSeq\n";
							print NA ">$name\n$retrievedRead\n";
						}else {
							++$defectivecount;
							print AADEFECTIVE ">$name\n$aaSeq\n";
							print NTDEFECTIVE ">$name\n$retrievedRead\n";
						}
					}else {
						++$outframecount;
						print OUTFRAME ">$name\n$retrievedRead\n";
					}
				}
			}
		}else {
			++$notPassRcutHitCount;
			print UNCOVER ">$readName\n$nameSeq{$readName}\n";
		}
	}
}
close XML;
close INFRAME;
close OUTFRAME;
close UNCOVER;
close DELETION;
close AA;
close AADEFECTIVE;
close NTDEFECTIVE;

print "total sequences: $seqcount\n";
print "Total queries: $queryCount, hits: $hitCount\n";
print "hits cover the region of $startPos - $endPos: $rCutoffHitCount, large deletion: $ldeletioncount, in frame: $inframecount, out frame: $outframecount, not cover: $notPassRcutHitCount\n";
print "complete AA sequences: $completecount, defective AA sequence: $defectivecount\n";

if ($dlx) {
	print "removing $xml ...\n";
	unlink $xml;
}


sub ReverseComplement {
	my $seq = shift;
	my $rcSeq = reverse $seq;
	$rcSeq =~ tr/ACGTacgt/TGCAtgca/;
	return $rcSeq;
}

sub translation {
	my $codon = shift;
	my %codon2aa = (
		"ATT" => "I",
		"ATC" => "I",
		"ATA" => "I",
		"CTT" => "L",
		"CTC" => "L",
		"CTA" => "L",
		"CTG" => "L",
		"TTA" => "L",
		"TTG" => "L",
		"GTT" => "V",
		"GTC" => "V",
		"GTA" => "V",
		"GTG" => "V",
		"TTT" => "F",
		"TTC" => "F",
		"ATG" => "M",
		"TGT" => "C",
		"TGC" => "C",
		"GCT" => "A",
		"GCC" => "A",
		"GCA" => "A",
		"GCG" => "A",
		"GGT" => "G",
		"GGC" => "G",
		"GGA" => "G",
		"GGG" => "G",
		"CCT" => "P",
		"CCC" => "P",
		"CCA" => "P",
		"CCG" => "P",
		"ACT" => "T",
		"ACC" => "T",
		"ACA" => "T",
		"ACG" => "T",
		"TCT" => "S",
		"TCC" => "S",
		"TCA" => "S",
		"TCG" => "S",
		"AGT" => "S",
		"AGC" => "S",
		"TAT" => "Y",
		"TAC" => "Y",
		"TGG" => "W",
		"CAA" => "Q",
		"CAG" => "Q",
		"AAT" => "N",
		"AAC" => "N",
		"CAT" => "H",
		"CAC" => "H",
		"GAA" => "E",
		"GAG" => "E",
		"GAT" => "D",
		"GAC" => "D",
		"AAA" => "K",
		"AAG" => "K",
		"CGT" => "R",
		"CGC" => "R",
		"CGA" => "R",
		"CGG" => "R",
		"AGA" => "R",
		"AGG" => "R",
		"TAA" => "*",
		"TAG" => "*",
		"TGA" => "*",
	);
	return $codon2aa{$codon};
}

sub deambiguous {
	my $ambiguity = shift;
	my %deambi = (
		"R" => "AG",
		"Y" => "CT",
		"S" => "CG",
		"W" => "AT",
		"K" => "GT",
		"M" => "AC",
		"B" => "CGT",
		"V" => "ACG",
		"D" => "AGT",
		"H" => "ACT",
		"N" => "ACGT",
	);
	return $deambi{$ambiguity};
}
