#!/usr/bin/env perl

use strict; use warnings;use threads;
use Mods::IO_Tamoc_progs qw(getProgPaths buildMapperIdx);
my $smtBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";


my ($inBam,$smplTag,$doBody) = @ARGV;
die "undefined contig marker" if ($smplTag eq "");
#die if (system "rm $inBam.bai; $smtBin index $inBam" );
return "$inBam" if (-e "$inBam.hd");
my $tmp= `$smtBin view -H $inBam `;
#die "$tmp\n";
my @XX = split '\n',$tmp;
#open I,">$inBam.sam";
foreach my $l (@XX){
	if ($l =~ m/^\@SQ/){
		next unless ($l =~ m/SN:$smplTag/);
	}
	print $l."\n";
}
#close I;
#system "$smtBin reheader $inBam.sam $inBam > $inBam.tmp"; 
if ($doBody){
system "$smtBin view $inBam >> $inBam.sam;";
system "$smtBin  view -b $inBam.sam > $inBam; $smtBin  index $inBam";
system "rm $inBam.sam;touch $inBam.hd";
}
#die "$inBam\n";
#return "$inBam";

exit (0);