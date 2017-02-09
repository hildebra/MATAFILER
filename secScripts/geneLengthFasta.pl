#!/usr/bin/env perl
use Mods::GenoMetaAss qw(gzipwrite gzipopen readFasta);

my ($inF,$outF) = @ARGV;

my $fr = readFasta($inF,1);
my %FA = %{$fr};

open O,">$outF" or die "Can;t open output DB length file $outF\n";
foreach my $k (keys %FA){
	print O "$k\t".length($FA{$k})."\n";
}
close O;

print "Done\n";