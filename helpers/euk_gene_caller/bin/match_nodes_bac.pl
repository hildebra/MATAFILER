#!/usr/bin/perl -w
use strict;

#this gets your reference table
my $ref = $ARGV[0];
my $fasta = $ARGV[1];
my %reffile = &build_seq_info_hash($ref);
#print_hash(%reffile); 
#exit;


	my $output_file = "$fasta.renamed.fa";
	open (OUTPUT, ">$output_file");
	
	open (FASTA, $fasta);
			
	while (my $header = <FASTA>) {
		chomp($header);
		my $sequence = <FASTA>;
		chomp($sequence);
		
		my ($seqid) = split (/\t/, $header);
		my $match = $reffile{$seqid} || "Broken";
		print OUTPUT "$match\n$sequence\n";

}

exit;


####################

sub build_seq_info_hash{
	my $file = shift;
	my %hash;

	open(IN, $file);
	foreach (<IN>) {
		chomp;
		my($a, $b) = split(/\t/);
		
	$a =~ s/^\s+//;
	$a =~ s/\s+$//;
	
	$b =~ s/^\s+//;
	$b =~ s/\s+$//;
	
	$hash{$a} = $b;
	
					}
	close(IN);
	return %hash;
}

#subroutine of general use; prints a hash out
sub print_hash{
  my(%hash) = @_;
  foreach my $key(sort keys %hash) {
    print "'$key'\t'$hash{$key}'\n";
  }
}

