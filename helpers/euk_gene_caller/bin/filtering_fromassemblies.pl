#!/usr/bin/perl -w
use strict;


foreach my $file (@ARGV) {
open (FILE, $file);

my $output_file = "$file.sub1000.fa";
my $output_file2 = "$file.gt1000.af";

open (OUTPUT, ">$output_file");
open (OUTPUT2, ">$output_file2");

	while (my $header = <FILE>) {	
		chomp($header);
		my $sequence = <FILE>;
		chomp($sequence);


#let's filter on length
	my($node, $nodenum, $length, $lengthnum, $cov, $covnum) = split(/\t/);
	$length = $length +36;
	
	if ($length <= 1000) {
	print OUTPUT "$header\n$sequence\n";
	}

	else {
	print OUTPUT2 "$header\n$sequence\n";
	}

	}}


exit;
