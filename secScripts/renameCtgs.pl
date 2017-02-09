#!/usr/bin/perl
#./renameCtgs.pl /g/scb/bork/hildebra/SNP/GNMass/alien-11-2-0/assemblies/metag/scaffolds.fasta
use warnings;
use strict;

sub newHD;

my $inF = $ARGV[0];
my $tag = $ARGV[1];
my $tmpOut = $inF.".tmp";
my $transOut = $inF.".lnk";

open I,"<$inF";
open O,">$tmpOut";
open O2,">$transOut"; my $ohd = ""; my $ntag = ""; my $seq = "";
my $cnt = 0; 
my $line = <I>;
$ntag = newHD($line,$seq);
$cnt++;
while ($line = <I>){
	chomp $line;
	if ($line =~ m/>/){
		chomp $line; $ohd = $line;
		#next if (length($seq) ==0);
		$ntag = newHD($line,$seq);
		$seq =~ s/(.{1,80})/$1\n/gs;
		print O "$ntag\n$seq"; 
		print O2 "$ntag\t$ohd\n";
		$cnt++;
		$seq = "";
	} else {
		$seq .= $line;
	}
}
$seq =~ s/(.{1,80})/$1\n/gs;
$ntag = newHD("",$seq);
print O "$ntag\n$seq"; $cnt++;
print O2 "$ntag\t$ohd\n";
close I; close O; close O2;

system("rm $inF; mv $tmpOut $inF");


sub newHD($ $){
	my ($line, $seq) = @_;
	my $LendTag = "="; #;
	if ($line =~ m/^>.*_length_(\d+).*/){
		$ntag = ">$tag"."__C$cnt"."_L=$1$LendTag";
	} else {
		$ntag = ">$tag"."__C$cnt"."_L=".length($seq).$LendTag;
	}

}
