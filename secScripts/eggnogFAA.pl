#!/usr/bin/env perl
use warnings;
use strict;


#/g/bork1/huerta/eggnog_hmm

my $hmmBin = "/g/bork1/huerta/eggnog_hmm/hmmer-3.1b1-linux-intel-x86_64_PATCHED/src/hmmscan";
my $inProt = "/g/scb/bork/hildebra/SNP/GNMass/alien-11-0-0/assemblies/metag/genePred/proteins.shrt.shrtHD.faa";
my $hmmDB = "/g/bork1/huerta/eggnog_hmm/bactNOG/bactNOG.hmm";
my $hitOut = "/g/scb/bork/hildebra/SNP/GNMass/alien-11-0-0/assemblies/metag/genePred/protEggNOG.txt";

my $cmd = "$hmmBin --cpu 10 -o /dev/null --tblout $hitOut  $hmmDB    $inProt";
die $cmd ."\n";

