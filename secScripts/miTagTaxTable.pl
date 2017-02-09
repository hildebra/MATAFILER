#! /usr/bin/perl -w
#perl site_taxon_table.pl Family * > family.ITS.tab ### usage

use strict;

my $tax_level_a= lc shift @ARGV;
my $outF = shift @ARGV;
my @tlvls = split(/,/,$tax_level_a);
my %sites;
my %taxa;

foreach my $file (@ARGV) {
        open my $FHANDLE, "<$file" or die "Could not open $file: $!\n";
		my $tag = $file;
		$tag =~ s/.*\///;
		$tag =~ s/\.hiera\.txt$//;
        my %column ; my $cset = 0;
        while (my $row=<$FHANDLE>) {
                chomp $row;
				my @temp = split /\t/, $row;
				if ($cset == 0) {
					foreach my $l (@tlvls){
						for (my $i=0; $i<scalar @temp; ++$i) {
								if (lc $temp[$i] eq $l) { $column{$l} = ($i-1); last; } #new LCA: get rid of first entry
								if ($i+1 == scalar @temp) { die "Could not find given taxon level $l in $file\n"; }
						}
						$sites{$l}{$tag} = {};
					}
					$cset=1;
				} else {
					shift @temp; #new LCA: get rid of first entry
					foreach my $l (@tlvls){
						#die "\n".$column+1 ."\n";
						my $k = "";
						if (@temp <= $column{$l}){
							$k = join (';',@temp[0 .. $#temp]);
							$k .= join(";", "?" x ($column{$l} - $#temp));
						} else {
							$k = join (';',@temp[0 .. $column{$l}]);
						}
						#print "$column{$l} @temp\n$k\n$file\n";
						#die $k."\n";
						++$sites{$l}{$tag}->{$k};
						++$taxa{$l}{$k};
					}
				}
        }
        close $FHANDLE;
}
foreach my $l (@tlvls){
	my @taxa_keys = sort keys %{$taxa{$l}};
	my @sites_keys = sort keys %{$sites{$l}};
	open O,">$outF.$l.txt";
	foreach my $site (@sites_keys) { print O "\t$site"; }
	print O "\n";
	foreach my $key (@taxa_keys) {
		print O $key;
		foreach my $site (@sites_keys) { 
			if ($sites{$l}{$site}->{$key}) { print O "\t$sites{$l}{$site}->{$key}"; }
			else { print O "\t0"; }
		}
		print O "\n";
	}
	close O;
}