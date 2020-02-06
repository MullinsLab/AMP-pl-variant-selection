#!/usr/bin/perl -w

#########################################################################################################
# Program: AaNtVariantRelation.pl
# Purpose: From _aa_variant_member.csv and _nt_variant_member.cvs, get aa and nt variant relationship
# Author: Wenjie Deng
# Date: 2019-04-01
###########################################################################################################

use strict;

my $usage = "perl AaNtVariantRelation.pl aa_member_file nt_member_file outfile\n";
my $inaafile = shift || die $usage;
my $inntfile = shift or die $usage;
my $outfile = shift or die $usage;
my (%nameVariant, %variantStatus);

open NT, $inntfile or die "couldn't open $inntfile: $!\n";
while (my $line = <NT>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my ($variant, $namelist) = split /\,/, $line;
	my @names = split /\;/, $namelist;
	foreach my $name (@names) {
		$nameVariant{$name} = $variant;
	}
}
close NT;

open AA, $inaafile or die "couldn't open $inaafile: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $line = <AA>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my @ntvariants = ();
	my ($aavariant, $namelist) = split /\,/, $line;
	my @names = split /\;/, $namelist;
	foreach my $name (@names) {
		my $ntvariant = $nameVariant{$name};
		if (!$variantStatus{$ntvariant}) {
			$variantStatus{$ntvariant} = 1;
			push @ntvariants, $ntvariant;
		}
	}
	print OUT $aavariant,',',join(';', @ntvariants),"\n";
}
close AA;
close OUT;
