#!/usr/bin/perl
#./calcGC.pl /g/bork1/hildebra/SNP/GNMassSimu2/AssmblGrp_1/metag/scaffolds.fasta.filt  /g/bork1/hildebra/SNP/GNMassSimu2/AssmblGrp_1/metag/ContigStats/sgc.2
use warnings;
use strict;

my $inF = $ARGV[0];
my $outF = $ARGV[1];

open I,"<$inF" or die "Can't open $inF";
open O,">$outF" or die "Can't open $outF";
print O "contig\tGC\n";

my $curTag = ""; my $GC=0; my $AT=0;
my $cnt=0;
while (my $line = <I>){
	if ($line =~ m/^>(.*)/){
		if ($cnt > 0){ print O "$curTag\t". sprintf('%.3f', ($GC/($GC+$AT)*100)) ."\n";}
		$curTag = $1; $GC=0;$AT=0;
		$cnt++; next;
	}
	$GC += ($line =~ tr/G//);	$GC += ($line =~ tr/C//);
	$AT += ($line =~ tr/A//);	$AT += ($line =~ tr/T//);
}
print O "$curTag\t". sprintf('%.3f', ($GC/($GC+$AT)*100)) ."\n";
close O; close I;