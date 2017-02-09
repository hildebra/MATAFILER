#!/usr/bin/perl -w
use strict;

#this gets your reference table
my $file = $ARGV[0];
my $ref1 = $ARGV[1];
my $ref2 = $ARGV[2];

my $ref3 = "";
$ref3 = $ARGV[3] if (@ARGV > 3);
my %reffile1 = &build_seq_info_hash($ref1);
my %reffile2 = &build_seq_info_hash($ref2);
my %reffile3;
%reffile3 = &build_seq_info_hash3($ref3) if ($ref3 ne "");
# print_hash(%reffile2); 
# exit;
 
	my $output_file = "$file.ALLmatched.txt";
	open (OUTPUT, ">$output_file");
	open (FILE, $file);

	foreach (<FILE>) {
 		chomp;
		my $matchmgm = $reffile1{$_} || "nocall";
		my $matchaug = $reffile2{$_} || "nocall";
		my $taxonomy = $reffile3{$_} || "nocall";
		
# 		if ($matchmgm eq "nocall" || $matchaug eq "nocall"){		
		print OUTPUT "$_\t$matchmgm\t$matchaug\t$taxonomy\n";
# }

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

	my ($geneid, $seqname) = split (/:/, $b);
	$seqname = join("", ">", $seqname);
	$hash{$seqname} = $b;
	
					}
	close(IN);
	return %hash;
}

sub build_seq_info_hash3{
	my $file = shift;
	my %hash;

	open(IN, $file) or die "Can't open $file (build_seq_info_hash3)\n";
	foreach (<IN>) {
		chomp;
		my($a, $b, $c) = split(/\t/);
		
	$a =~ s/^\s+//;
	$a =~ s/\s+$//;
	
	$b =~ s/^\s+//;
	$b =~ s/\s+$//;

	$c =~ s/^\s+//;
	$c =~ s/\s+$//;

	my $seqname = join("", ">", $a);
	$hash{$seqname} = $c;
	
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

