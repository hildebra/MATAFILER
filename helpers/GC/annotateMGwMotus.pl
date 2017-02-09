#!/usr/bin/env perl
#gets based on luis motu clustering script the genomes
#  ./annotateMGwMotus.pl /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR 12
use warnings;
use strict;

use Mods::GenoMetaAss qw(readClstrRev );
use Mods::IO_Tamoc_progs qw(getProgPaths);

sub readMotuTax;
sub readGene2mlinkage;
sub lambdaBl;


my $SpecID="/g/bork3/home/hildebra/DB/MarkerG/specI/";
#progenomes.specIv2_2


my $rarBin = getProgPaths("rare");#"/g/bork5/hildebra/dev/C++/rare/rare";
my $lambdaBin = getProgPaths("lambda");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda";
my $lambdaIdxBin = $lambdaBin."_indexer";#getProgPaths("");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda_indexer";

my $GCd = $ARGV[0];
my $BlastCores = $ARGV[1];
my $MGdir = "$GCd/FMG/";
system "mkdir -p $MGdir/tax" unless (-d "$MGdir/tax");

my $motuDir = "/g/bork3/home/hildebra/DB/MarkerG/mOTU";
#load motu DBs...
my ($hr1,$hr2) = readMotuTax("$motuDir/mOTU.v1.1.padded.motu.linkage.map");
my %LG2motu = %{$hr1}; my %motu2tax = %{$hr2};
$hr1 = readGene2mlinkage("$motuDir/mOTU.v1.1.padded.motu.map");
my %gene2LG = %{$hr1};
#annotate against DB using lambda

if (0){ #too general
	my $tar = "$GCd/compl.incompl.95.fna"; my $DB = "$motuDir/263MetaRef10MGv9.cal.v2.nr.padded.fna";	my $taxblastf = "$GCd/compl.incompl.95.motuAss.tmp.m8";
	lambdaBl($tar,$DB,$taxblastf);
}

#assign each COG separately
my @catsPre = split/\n/,`cat $GCd/FMG.subset.cats`; my %cats;
system "mkdir -p $MGdir" unless (-d $MGdir);
foreach (@catsPre){
	my @spl  = split /\t/;
	#$cats{$spl[0]} = $spl[2];
	my $samBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";
	my $cmd = "$samBin faidx $GCd/compl.incompl.95.fna ". join (" ", split(/,/,$spl[2]) ) . " > $MGdir/$spl[0].fna\n";
	system $cmd unless (-e "$MGdir/$spl[0].fna");
	$cmd = "$samBin faidx $GCd/compl.incompl.95.prot.faa ". join (" ", split(/,/,$spl[2]) ) . " > $MGdir/$spl[0].faa";
	system $cmd unless (-e "$MGdir/$spl[0].faa");
	
	lambdaBl("$MGdir/$spl[0].fna","$SpecID/$spl[0].rep.fna","$MGdir/tax/$spl[0].tmp.m8");
	

}

print "Done initial blast\n$MGdir\n";
exit(0);
#TOGO




#####################################################

sub lambdaBl($ $ $){
	my ($tar,$DB, $taxblastf) = @_;

	my $cmd="";
	if (!-f $DB.".dna5.fm.sa.val"  ) {
		print "Building LAMBDA index anew (may take up to an hour)..\n";
		my $cmdIdx = "$lambdaIdxBin -p blastn -t $BlastCores -d $DB";
		#system "touch $DB.dna5.fm.lf.drv.wtc.24";
		if (system($cmdIdx)){die ("Lamdba ref DB build failed\n$cmdIdx\n");}
	}
	$cmd .= "$lambdaBin -t $BlastCores -id 93 -nm 100 -p blastn -e 1e-40 -so 7 -sl 16 -sd 1 -b 5 --output-columns \"std qlen slen\" -pd on -q $tar -d $DB -o $taxblastf\n";
		#die $cmd."\n";
	system $cmd unless (-e $taxblastf);
}


sub readMotuTax($){
	my ($inF) = @_;
	my %gene2motu;
	my %motu2tax;
	open I,"<$inF";
	while (my $l = <I>){
		chomp $l; my @spl=split/\t/,$l;
		$gene2motu{$spl[0]} = $spl[8];
		$motu2tax{$spl[8]} = join(";",@spl[1,2,3,4,5,6,7]) if (!exists $motu2tax{$spl[8]} );
	}
	close I;
	return (\%gene2motu,\%motu2tax);
}

sub readGene2mlinkage($){
	my ($inF) = @_;
	my %gene2LG;
	open I,"<$inF";
	while (my $l = <I>){
		chomp $l; my @spl=split/\t/,$l;
		$gene2LG{$spl[0]} = $spl[2];
	}
	close I;
	return (\%gene2LG);
}














