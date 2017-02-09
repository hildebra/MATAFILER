#!/usr/bin/env perl
#reads in NCBI tax and uses LCA to derrive taxonomy of read / contig
use warnings;
use strict;
use Mods::GenoMetaAss qw(readNCBItax);


my $taxFile = "/g/bork3/home/hildebra/DB/eggNOG10/all_species_data.txt";
my $hr = readNCBItax($taxFile);