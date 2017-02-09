#!/usr/bin/perl -w
use strict;


foreach my $file (@ARGV) {
open (FILE, $file);

my $output_file = "$file.gt1000genes.txt";
my $output_file2 = "$file.leftovers.txt";

open (OUTPUT, ">$output_file");
open (OUTPUT2, ">$output_file2");

	foreach (<FILE>) {	
 	chomp;

#let's filter on length

	my ($genecall, $seqid) = split(/:/);
	my($node, $nodenum, $length, $lengthnum, $cov, $covnum) = split(/_/, $seqid);
 	my $add = $lengthnum + 36;
	
# print "$add\n";
	
	if ($add > 1000) {
	print OUTPUT "$_\n";
	}

	else {
	print OUTPUT2 "$_\n";
	}

	}}


exit;
