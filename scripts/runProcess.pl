#!/usr/bin/perl -w


use strict;
use Getopt::Long;
use File::Basename;
use File::Copy;

my %option = (
	'sid'     => '',
	'outpath' => '.',
	'in'      => '',
	'sanger'  => '',
	'envref'  => 'panels/HXB2_5970-9012.fasta',
	'reref'   => 'panels/hxb2-re.fasta',
	'rs'      => 256,
	're'      => 2826,
	'scripts' => 'scripts/',
	'blast'   => '/usr/local/bin',
	'h'       => '',
);

my $usage = "\nUsage: perl runProcess.pl [-option value]

options:
-sid     sample ID
-outpath output data path
-in      input PacBio consensus fasta file
-sanger  input sanger sequences fasta file
-envref  HXB2 env reference fasta file
-reref   HXB2 RE reference fasta file
-rs      region start position in reference sequence (default: 1, beginning of reference)
-re      region end position in reference sequence (default: 0, end of reference)
-scripts where perl scripts are ()
-blast   where blastn and makeblastdb are
-h       usage help

";

GetOptions (\%option, 'sid=s', 'outpath=s', 'in=s', 'sanger=s', 'envref=s', 'reref=s', 'rs=i', 're=i', 'scripts=s', 'blast=s', 'h');
my $sid = $option{'sid'} or die $usage;
my $outpath = $option{'outpath'} or die $usage;
my $infasta = $option{'in'} or die $usage;
my $sangerfile = $option{'sanger'};
my $refFasta = $option{'envref'} or die $usage;
my $reref = $option{'reref'} or die $usage;
my $startPos = $option{'rs'};
my $endPos = $option{'re'};
my $scripts = $option{'scripts'};
my $blast = $option{'blast'};
my $help = $option{'h'};
die $usage if ($help);

unless (-e $outpath) {
	mkdir $outpath;
}

$sid =~ s/_REN_pe//; # sid is only the sample id
$sid =~ s/_REN_p//; # sid is only the sample id
my $infastabasename = basename($infasta);
$outpath =~ s/\/$//g;
my $outfasta = "$outpath/$infastabasename";
copy($infasta, $outfasta) or die "Copy failed: $!\n";
my $logfile =  my $outXML = $outfasta;
$logfile =~ s/.fasta/_log.txt/;
$outXML =~ s/\.fasta/.xml/;
my $rv = 0;
print "\n### input from: $infasta, output to: $outpath ###\n";
print "\n... running BLAST ... \n";
$rv = system("$scripts/runBLAST.pl -in $outfasta -ref $refFasta -blast $blast >> $logfile");
unless ($rv == 0) {
	die "runBLAST.pl failed: $rv\n";
}
print "\n... retrieving reads ...\n";
$rv = system("$scripts/retrieveRegion.pl -fas $outfasta -xml $outXML -ref $refFasta -rs $startPos -re $endPos >> $logfile");
unless ($rv == 0) {
	die "retrieveRegion.pl failed: $rv\n";
}

my $ntfile = my $aafile = $outfasta;
$ntfile =~ s/_REN_pe.fasta/_env_pe_NT.fasta/;
$ntfile =~ s/_REN_p.fasta/_env_p_NT.fasta/;
$aafile =~ s/_REN_pe.fasta/_env_pe_AA.fasta/;
$aafile =~ s/_REN_p.fasta/_env_p_AA.fasta/;
$ntfile =~ s/_GP_pe.fasta/_gag_pe_NT.fasta/;
$ntfile =~ s/_GP_p.fasta/_gag_p_NT.fasta/;
$aafile =~ s/_GP_pe.fasta/_gag_pe_AA.fasta/;
$aafile =~ s/_GP_p.fasta/_gag_p_AA.fasta/;
print "\n... calculating nucleotide sequence variants ...\n";
$rv = system("$scripts/uniqueReads.pl -sid $sid -in $ntfile -type NA >> $logfile");
unless ($rv == 0) {
	die "uniqueReads.pl failed: $rv\n";
}

print "\n... calculating amino acid sequence variants ...\n";
$rv = system("$scripts/uniqueReads.pl -sid $sid -in $aafile -type AA >> $logfile");
unless ($rv == 0) {
	die "uniqueReads.pl failed: $rv\n";
}

my $ntvariantfile = my $ntvariantalignfile = $ntfile;
$ntvariantfile =~ s/\.fasta/_variants.fasta/;
$ntvariantalignfile =~ s/\.fasta/_variants.aln/;
my $aavariantfile = my $aavariantalignfile = $aafile;
$aavariantfile =~ s/\.fasta/_variants.fasta/;
$aavariantalignfile =~ s/\.fasta/_variants.aln/;
my $sangerntvariantfile = my $sangeraavariantfile = '';
if ($sangerfile) {
	print "\n... running BLAST for sanger sequences ... \n";
	$rv = system("$scripts/runBLAST.pl -in $sangerfile -ref $refFasta -blast $blast >> $logfile");
	unless ($rv == 0) {
		die "runBLAST.pl failed: $rv\n";
	}

	my $outXML = $sangerfile;
	$outXML =~ s/\.fasta/.xml/;
	print "\n... retrieving sanger sequences ...\n";
	$rv = system("$scripts/retrieveRegion.pl -fas $sangerfile -xml $outXML -ref $refFasta -rs $startPos -re $endPos >> $logfile");
	unless ($rv == 0) {
		die "retrieveRegion.pl failed: $rv\n";
	}

	my $sangerntfile = $sangerfile;
	$sangerntfile =~ s/_REN_sanger\.fasta/_env_sanger_NT.fasta/;
	$sangerntfile =~ s/_GP_sanger\.fasta/_gag_sanger_NT.fasta/;
	my $sangeraafile = $sangerfile;
	$sangeraafile =~ s/_REN_sanger\.fasta/_env_sanger_AA.fasta/;
	$sangeraafile =~ s/_GP_sanger\.fasta/_gag_sanger_AA.fasta/;
	print "\n... calculating sanger nucleotide sequence variants ...\n";
	$rv = system("$scripts/uniqueReads.pl -sid $sid -in $sangerntfile -type SNA >> $logfile");
	unless ($rv == 0) {
		die "uniqueReads.pl failed: $rv\n";
	}

	print "\n... calculating sanger amino acid sequence variants ...\n";
	$rv = system("$scripts/uniqueReads.pl -sid $sid -in $sangeraafile -type SAA >> $logfile");
	unless ($rv == 0) {
		die "uniqueReads.pl failed: $rv\n";
	}

	$sangerntvariantfile = $sangerntfile;
	$sangerntvariantfile =~ s/\.fasta/_variants.fasta/;
	$sangeraavariantfile = $sangeraafile;
	$sangeraavariantfile =~ s/\.fasta/_variants.fasta/;
	if (-s $sangerntfile) {
		print "\n... merging sanger and consensus variants ...\n";
		open SANGERNT, $sangerntvariantfile or die "couldn't open $sangerntvariantfile: $!\n";
		open OUT, ">>", $ntvariantfile or die "couldn't open $ntvariantfile: $!\n";
		while (my $line = <SANGERNT>) {
			chomp $line;
			print OUT $line,"\n";
		}
		close SANGERNT;
		close OUT;

		open SANGERAA, $sangeraavariantfile or die "couldn't open $sangeraavariantfile: $!\n";
		open OUT, ">>", $aavariantfile or die "couldn't open $aavariantfile: $!\n";
		while (my $line = <SANGERAA>) {
			chomp $line;
			print OUT $line,"\n";
		}
		close SANGERAA;
		close OUT;
	}
}

print "\n... aligning nucleotide sequences ...\n";
$rv = system("muscle -quiet -in $ntvariantfile -out $ntvariantalignfile");
unless ($rv == 0) {
	die "muscle alignment failed: $rv\n";
}

print "\n... calculating nucleotide sequence mismatches ...\n";
if (-s $sangerntvariantfile) {
	my $sangerntvariantmemberfile = $sangerntvariantfile;
	$sangerntvariantmemberfile =~ s/_variants\.fasta/_variantList.csv/;
	$rv = system("$scripts/calcDistanceFromMostCommonVariant.pl -in $ntvariantalignfile -sname $sangerntvariantmemberfile >> $logfile");
}else {
	$rv = system("$scripts/calcDistanceFromMostCommonVariant.pl -in $ntvariantalignfile >> $logfile");
}
unless ($rv == 0) {
	die "calcDistanceFromMostCommonVariant.pl failed: $rv\n";
}

print "\n... aligning amino acid sequences ...\n";
$rv = system("muscle -quiet -in $aavariantfile -out $aavariantalignfile");
unless ($rv == 0) {
	die "muscle alignment failed: $rv\n";
}

print "\n... calculating amino acid sequence mismatches ...\n";
if (-s $sangeraavariantfile) {
	my $sangeraavariantmemberfile = $sangeraavariantfile;
	$sangeraavariantmemberfile =~ s/_variants\.fasta/_variantList.csv/;
	$rv = system("$scripts/calcDistanceFromMostCommonVariant.pl -in $aavariantalignfile -sname $sangeraavariantmemberfile >> $logfile");
}else {
	$rv = system("$scripts/calcDistanceFromMostCommonVariant.pl -in $aavariantalignfile >> $logfile");
}
unless ($rv == 0) {
	die "calcDistanceFromMostCommonVariant.pl failed: $rv\n";
}

my $aavariantmemberfile = my $aantvariantmemberfile = $aafile;
$aavariantmemberfile =~ s/\.fasta/_variantList.csv/;
$aantvariantmemberfile =~ s/\.fasta/_NT_variantList.csv/;
my $ntvariantmemberfile = $ntfile;
$ntvariantmemberfile =~ s/\.fasta/_variantList.csv/;
$rv = system("$scripts/AaNtVariantRelation.pl $aavariantmemberfile $ntvariantmemberfile $aantvariantmemberfile >> $logfile");
unless ($rv == 0) {
	die "AaNtVariantRelation.pl failed: $rv\n";
}

print "\n... retrieve > 9.5% variant RE sequences ...\n";
my $aavariantfreqfile = $aavariantmemberfile;
$aavariantfreqfile =~ s/List.csv/Freq.csv/;
$rv = system("$scripts/retrieveVariantRE.pl $outfasta $ntfile $aavariantfreqfile $aavariantmemberfile $reref >> $logfile");
unless ($rv == 0) {
	die "retrieveVariantRE.pl failed: $rv\n";
}

print "\n... create summary file ...\n";
my $countsummaryfile = $logfile;
$countsummaryfile =~ s/_log.txt/_count_summary.csv/;
$rv = system("$scripts/countSummary.pl $sid $logfile $aavariantfile $countsummaryfile");
unless ($rv == 0) {
	die "countSummary.pl failed: $rv\n";
}

print "\n... write sample summary file ...\n";
my @fields = split /\//, $outfasta;
pop @fields;
pop @fields;
my $samplesummaryfile = join ('/', @fields);
$samplesummaryfile .= "/sample_summary.csv";
$rv = system("$scripts/sampleSummary.pl $countsummaryfile $samplesummaryfile");
unless ($rv == 0) {
	die "sampleSummary.pl failed: $rv\n";
}

my $opidir = "$outpath/ofPossibleInterest/";
my $defectdir = "$outpath/Defectives/";
unless (-e $opidir) {
	mkdir $opidir;
}
unless (-e $defectdir) {
	mkdir $defectdir;
}
my $inframefile = my $largedeletionfile = my $ntdistfile = my $nthistfile = my $ntstopfile = my $ntoutframefile = $ntfile;
my $aadistfile = my $aahistfile = my $aastopfile = $aafile;
$inframefile =~ s/\.fasta/_inframe.fasta/;
$largedeletionfile =~ s/\.fasta/_largeDeletion.fasta/;
$ntdistfile =~ s/\.fasta/_variantDist.csv/;
$nthistfile =~ s/\.fasta/_variantHist.csv/;
$ntstopfile =~ s/\.fasta/_stop.fasta/;
$ntoutframefile =~ s/\.fasta/_outframe.fasta/;
$aadistfile =~ s/\.fasta/_variantDist.csv/;
$aahistfile =~ s/\.fasta/_variantHist.csv/;
$aastopfile =~ s/\.fasta/_stop.fasta/;
move ($ntoutframefile, $defectdir);
unlink $ntoutframefile;
move ($aastopfile, $defectdir);
unlink $aastopfile;
move ($ntstopfile, $defectdir);
unlink $ntstopfile;
move ($outXML, $opidir);
unlink $outXML;
move ($ntfile, $opidir);
unlink $ntfile;
move ($aafile, $opidir);
unlink $aafile;
move ($ntvariantfile, $opidir);
unlink $ntvariantfile;
move ($aavariantfile, $opidir);
unlink $aavariantfile;
move ($inframefile, $opidir);
unlink $inframefile;
move ($largedeletionfile, $opidir);
unlink $largedeletionfile;
move ($ntdistfile, $opidir);
unlink $ntdistfile;
move ($nthistfile, $opidir);
unlink $nthistfile;
move ($aadistfile, $opidir);
unlink $aadistfile;
move ($aahistfile, $opidir);
unlink $aahistfile;
my $RENTfiles = $ntfile;
$RENTfiles =~ s/_env_pe_NT\.fasta/_RE_NT_V\*/;
$RENTfiles =~ s/_env_p_NT\.fasta/_RE_NT_V\*/;
for my $file (glob $RENTfiles) {
    move ($file, $opidir) or die $!;
    unlink $file;
}

print "\n... done ...\n\n";
