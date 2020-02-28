#!/usr/bin/perl -w

#########################################################################################################
# Program: retrieveVariantRE.pl
# Purpose: retrieve RE nucleotide consensus sequence for the AA variants >= 9.5%
# Input: _REN_p.fasta, _AA_variant_freq.csv and _AA_variant_member.csv
# Output: fasta file
# Author: Wenjie Deng
# Date: 2019-08-15
###########################################################################################################

use strict;

my $usage = "perl retrieveVariantRE.pl input_REN_file input_env_NA_file input_aa_variant_freq_file input_aa_variant_member_file HXB2_RE_ref_file\n";
my $nafile = shift or die $usage;
my $envfile = shift or die $usage;
my $aafreqfile = shift or die $usage;
my $aamemberfile = shift or die $usage;
my $REreffile = shift or die $usage;

my $totalSeq = my $variantcount = my $flag = my $cleancount = 0;
my (@variants, %variantFreq, %variantUMIs);

open FREQ, $aafreqfile or die "couldn't open $aafreqfile: $!\n";
while (my $line = <FREQ>) {
	chomp $line;
	next if ($line =~ /^\s*$/ or $line =~ /Variant,Count,Frequency/);
	my ($name, $count, $freq) = split (/,/, $line);
	if ($freq > 0.095) {
		push @variants, $name;
		$variantFreq{$name} = $freq;
	}
}
close FREQ;

open MEM, $aamemberfile or die "couldn't open $aamemberfile: $!\n";
while (my $line = <MEM>) {
	chomp $line;
	next if ($line =~ /^\s*$/);
	my ($variant, $namestring) = split (/,/, $line);
	if ($variantFreq{$variant}) {
		my @names = split /\;/, $namestring;
		foreach my $name (@names) {
#			my @fields = split /_/, $name;
#			my $umi = $fields[$#fields];
			my $umi = "";
			if ($name =~ /_env_(\w+)/) {
				$umi = $1;
			}else {
				die "name not formatted: $name\n";
			}
			$variantUMIs{$variant}{$umi} = 1;
		}
	}
}
close MEM;

my $refname = my $refseq = "";
open REF, $REreffile or die "couldn't open $REreffile: $!\n";
while (my $line = <REF>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$refname = $1;
	}else {
		$refseq .= uc $line;
	}
}
close REF;

foreach my $variant (@variants) {
	print "variant >= 9.5%: $variant\n";
	my $umi = "";
	my $familycount = 0;
	my (%envUMI, %renUMI);
	open ENV, $envfile or die "couldn't open $envfile: $!\n";
	while (my $line = <ENV>) {
		chomp $line;
		if ($line =~ /^>(\S+)/) {
			++$familycount;
			$flag = 0;
			my $name = $1;
#			my @fields = split /_/, $name;
#			$umi = $fields[$#fields];
			if ($name =~ /_env_(\w+)/) {
				$umi = $1;
			}else {
				die "name not formatted: $name\n";
			}
			if ($variantUMIs{$variant}{$umi}) {
				$flag = 1;
			}
		}elsif ($flag) {
			$envUMI{$umi} .= $line;
		}
	}
	close ENV;
	open REN, $nafile or die "couldn't open $nafile: $!\n";
	while (my $line = <REN>) {
		chomp $line;
		if ($line =~ /^>(\S+)/) {
			$flag = 0;
			my $name = $1;
#			my @fields = split /_/, $name;
#			$umi = $fields[$#fields];
			if ($name =~ /_REN_(\w+)/) {
				$umi = $1;
			}else {
				die "name not formatted: $name\n";
			}
			if ($variantUMIs{$variant}{$umi}) {
				$flag = 1;
			}
		}elsif ($flag) {
			$line =~ s/\-//g;
			$renUMI{$umi} .= uc $line;
		}
	}
	close REN;
	my $REfile = $envfile;
	$REfile =~ s/\w+\.fasta/$variant\.fasta/;
	$REfile =~ s/_env_AA_/_RE_NT_/;
	$REfile =~ s/_\d+\.fasta/\.fasta/;
	open RE, ">", $REfile or die "couldn't open $REfile: $!\n";
	foreach my $umi (keys %renUMI) {
		if ($envUMI{$umi}) {
			my $idx = index($renUMI{$umi}, $envUMI{$umi});
			my $RElen = $idx + length($envUMI{$umi});
			my $re = substr($renUMI{$umi}, 0, $RElen);
			print RE ">$umi\n$re\n";
		}else {
			die "No UMI: $umi in envUMI\n";
		}
	}
	close RE;

	# muscle alignment
	my $REalignfile = my $REconsfile = $REfile;
	$REalignfile =~ s/.fasta/.aln/;
	system("muscle -quiet -in $REfile -out $REalignfile");
	# calculate consensus
	my $REcons = calculate_consensus($REalignfile, 0.5);
	my $REvariant = $REalignfile;
#	$REvariant =~ s/.aln//;
	$REvariant =~ /([\w\.]+)\.aln/;
	$REvariant = $1;
	print "$REvariant consensus: $REcons\n";
	$REconsfile =~ s/.fasta/_cons.fasta/;
	open RECONS, ">", $REconsfile or die "couldn't open $REconsfile: $!\n";
	print RECONS ">$refname\n$refseq\n";
	print RECONS ">consensus\n$REcons\n";
	close RECONS;
	# check rev in frame and no internal stop codon
	my $REconsalignfile = $REconsfile;
	$REconsalignfile =~ s/.fasta/.aln/;
	system("muscle -quiet -in $REconsfile -out $REconsalignfile");
	my $alignrefseq = my $alignconsseq = "";
	open ALN, $REconsalignfile or die "couldn't open $REconsalignfile: $!\n";
	while (my $line = <ALN>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>(\S+)/) {
			my $name = $1;
			if ($name =~ /HXB2/) {
				$flag = 1;
			}else {
				$flag = 0;
			}
		}elsif ($flag) {
			$alignrefseq .= $line;
		}else {
			$alignconsseq .= $line;
		}
	}
	close ALN;
	my $alignlen = length $alignrefseq;
	my @alignrefnas = split //, $alignrefseq;
	my @alignconsnas = split //, $alignconsseq;
	my $idx = 0;
	my $conscodon = "";
	for (my $i = 0; $i < $alignlen; $i++) {
		if ($alignrefnas[$i] =~ /[ACGT]/) {
			++$idx;
		}
		if ($idx <= 75) { # rev1
			$conscodon .= $alignconsnas[$i];
		}
		if ($idx == 75) {
			last;
		}
	}
	$conscodon =~ s/\-//g;
	my $conscodonlen = length $conscodon;
	if ($conscodonlen % 3 == 0) { # in frame
		my @nas = split //, $conscodon;
		for (my $i = 0; $i < length $conscodon; $i += 3) {
			my $codon = $nas[$i].$nas[$i+1].$nas[$i+2];
			my $aa = translation($codon);
			if ($aa eq "*") {
				die "rev1 has internal stop codon $codon\n";
			}
		}
	}else {
		die "rev1 out of frame: $conscodon, $conscodonlen\n";
	}
	my $variantSynthesisFile = $REfile;
	$variantSynthesisFile =~ s/_RE_NT_V/_RE_pb/;
	$variantSynthesisFile =~ s/_RE_(.*?)_(.*?)\.fasta/_RE_$1.fasta/;
	$variantSynthesisFile =~ s/\.fasta/_s.fasta/;
	my $synthseqname = $variantSynthesisFile;
	if ($synthseqname =~ /(\w+)\.fasta/) {
		$synthseqname = $1;
	}
	my $varidx = "";
	if ($synthseqname =~ /_RE_pb(.*?)_s/) {
		$varidx = $1;
	}
	$synthseqname .= " Variant=$varidx Frequency=$variantFreq{$variant} Total#families=$familycount";
	open SYNTH, ">", $variantSynthesisFile or die "couldn't open $variantSynthesisFile: $!\n";
	print SYNTH ">$synthseqname\n$REcons\n";
	close SYNTH;
}



sub calculate_consensus {
	my $file = shift;
	my $cutoff = shift;
	my $count = 0;
	my $cons = my $seq = "";
	my (@seqs, %naPosCount);
	my @nas = ("A", "C", "G", "T", "-");
	open IN, $file or die "couldn't open $file: $!\n";
	while (my $line = <IN>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>/) {
			if ($count) {
				push @seqs, $seq;
			}
			++$count;
			$seq = "";
		}else {
			$seq .= uc $line;
		}
	}
	# last seq
	my $alignlen = length $seq;
	push @seqs, $seq;
	close IN;
	foreach my $seq (@seqs) {
		my @seqnas = split //, $seq;
		for (my $i = 0; $i < $alignlen; $i++) {
			++$naPosCount{$i}{$seqnas[$i]};
		}
	}
	for (my $i = 0; $i < $alignlen; $i++) {
		my $nacons = "";
		foreach my $na (@nas) {
			if ($naPosCount{$i}{$na}) {
				my $nafreq = $naPosCount{$i}{$na} / $count;
				if ($nafreq >= $cutoff) {
					$nacons = $na;
					last;
				}
			}
		}
		if ($nacons) {
			$cons .= $nacons;
		}else {
			die "at positon $i, no na exceed cutoff of $cutoff";
		}
	}
	return $cons;
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
