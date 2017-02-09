#!/usr/bin/perl -w
use strict;
	
			
	foreach my $file (@ARGV) {
	open (FILE, $file);

	my $output_file = "$file.output.txt";
	open (OUTPUT, ">$output_file");

	while (<FILE>) {
 	chomp;
		$_ =~ s/^\s+//;
		$_ =~ s/\s+$//;

		my ($seqname, $caller, $codingtype, $start, $stop, $score, $strand, $frame, $geneid) = split (/\t/);
 		$geneid =~ tr/ /_/;
 		$geneid =~ s/^/>/;
  		my $new = join(':', $geneid, $seqname);

  	  print OUTPUT "$geneid\t$new\n";
	}
}

exit;




			
