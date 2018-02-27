#!/usr/bin/perl
use warnings;
use strict;
use Mods::GenoMetaAss qw(gzipwrite gzipopen);

my $inF = $ARGV[0];
my $outF = $ARGV[1];
my $mode = $ARGV[2];


my ($I,$OK) = gzipopen($inF,"mp3 output file [in]",1); 
my($O,$OK2) = gzipwrite($outF,"mp3 [out]",1);

while (my $line = <$I>){
	next unless ($line =~ m/\d+\./);
	chomp $line;
	my @spl = split /\t/,$line;
	next unless (@spl > 1); chomp $spl[1];
	#my $classi = "";
	#if ($spl[3] eq "Pathogenic" || $spl
	print $O "$spl[1]\t$spl[6]\n";
}

close $I; close $O;
