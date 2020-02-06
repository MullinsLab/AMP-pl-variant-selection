#!/usr/bin/perl -w

#######################################################################################
# Calculate number of point mutations between each variant and the most common variant.
# Output distance list file and histogram file
#######################################################################################

use strict;
use Getopt::Long;

my %option = (
	'in'    => '',
	'sname' => '',
);

my $usage = "\nusage: perl calcDistanceFromMostCommonVariant.pl [-option value]

options:  
-in		input sequence alignment fasta file
-sname  input sanger variant member file
";

GetOptions (\%option, 'in=s', 'sname=s');

my $inFile = $option{'in'} or die $usage;
my $snamefile = $option{'sname'};
my $outDistFile = my $outHistFile = $inFile;
$outDistFile =~ s/_variants\.aln/_variantDist.csv/;
$outHistFile =~ s/_variants\.aln/_variantHist.csv/;
my (@names, %nameSeq, %nameSeqcount, @mostcommonNAs, %distVariantcount, %distSeqcount, %snameStatus, @snames, %snameNAs);
my $count = my $alignlen = my $seqcount = 0;
my $seq = my $seqName = my $mostcommonvariant = '';

if (-s $snamefile) {
	open SNAME, $snamefile or die "couldn't open $snamefile: $!\n";
	while (my $line = <SNAME>) {
		chomp $line;
		my ($svname, $members) = split /\,/, $line;
		$snameStatus{$svname} = 1;
		push @snames, $svname;
	}
	close SNAME;
}

open IN, $inFile or die "couldn't open $inFile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if ($line =~ /^\s*$/);
	if ($line =~ /^>(.*)$/) {
		if ($count == 1) {
			$alignlen = length $nameSeq{$seqName};
		}
		$seqName = $1;
		unless ($snameStatus{$seqName}) {
			push @names, $seqName;
			if ($seqName =~ /_(\d+)$/) {
				$nameSeqcount{$seqName} = $1;
				$seqcount += $1;
			}
			++$count;
		}				
	}else {
		$nameSeq{$seqName} .= uc $line;
	}		
}
close IN;

foreach my $name (@names) {
	if ($name =~ /_V001_/) {
		$mostcommonvariant = $name;
		@mostcommonNAs = split //, $nameSeq{$name};
		last;
	}
}

open DIST, ">", $outDistFile or die "couldn't open $outDistFile: $!\n";
print DIST "Reference,Variant,Mismatches\n";
foreach my $name (@names) {
	unless ($name =~ /_V001_/) {
		my $dist = 0;
		my @nas = split //, $nameSeq{$name};
		for (my $i = 0; $i < $alignlen; $i++) {
			if ($nas[$i] =~ /[A-Za-z]/ and $mostcommonNAs[$i] =~ /[A-Za-z]/ and $nas[$i] ne $mostcommonNAs[$i]) {
				++$dist;
			}
		}
		print DIST "$mostcommonvariant,$name,$dist\n";
		++$distVariantcount{$dist};
		$distSeqcount{$dist} += $nameSeqcount{$name};
	}	
}
print DIST "$mostcommonvariant,$mostcommonvariant,0\n";
++$distVariantcount{0};
$distSeqcount{0} += $nameSeqcount{$mostcommonvariant};
print DIST "\n";
close DIST;

open HIST, ">", $outHistFile or die "couldn't open $outHistFile: $!\n";
print HIST "Reference,Mismatches,Variant_count,Sequence_count,Sequence_fraction\n";
foreach my $dist (sort{$a <=> $b} keys %distVariantcount) {
	my $fraction = int($distSeqcount{$dist} / $seqcount * 10000 + 0.5) / 10000;
	print HIST "$mostcommonvariant,$dist,$distVariantcount{$dist},$distSeqcount{$dist},$fraction\n";
}
print HIST "\n";
close HIST;

if (@snames) {
	foreach my $sname (@snames) {
		@{$snameNAs{$sname}} = split //, $nameSeq{$sname};
	}
	open DIST, ">>", $outDistFile or die "couldn't open $outDistFile: $!\n";
	open HIST, ">>", $outHistFile or die "couldn't open $outHistFile: $!\n";
	foreach my $sname (@snames) {
		%distVariantcount = %distSeqcount = ();
		foreach my $name (@names) {
			my $dist = 0;
			my @nas = split //, $nameSeq{$name};
			for (my $i = 0; $i < $alignlen; $i++) {
				if ($nas[$i] =~ /[A-Za-z]/ and $snameNAs{$sname}[$i] =~ /[A-Za-z]/ and $nas[$i] ne $snameNAs{$sname}[$i]) {
					++$dist;
				}
			}
			print DIST "$sname,$name,$dist\n"; ###############
			++$distVariantcount{$dist};
			$distSeqcount{$dist} += $nameSeqcount{$name};
		}
		foreach my $dist (sort{$a <=> $b} keys %distVariantcount) {
			my $fraction = int($distSeqcount{$dist} / $seqcount * 10000 + 0.5) / 10000;
			print HIST "$sname,$dist,$distVariantcount{$dist},$distSeqcount{$dist},$fraction\n";
		}
		print DIST "\n";
		print HIST "\n";
	}
	close DIST;
	close HIST;
}

print "total $seqcount sequences, $count variants\n";
