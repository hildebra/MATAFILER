#!/usr/bin/env perl
#gets based on luis motu clustering script the genomes
#  ./annotateMGwMotus.pl /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR 12
use warnings;
use strict;

use Statistics::RankCorrelation;


use Mods::GenoMetaAss qw(readClstrRev );
use Mods::IO_Tamoc_progs qw(getProgPaths);

sub readMotuTax;
sub readGene2mlinkage;
sub lambdaBl;
sub readNCBI;
sub calculate_spearman_correlation;
sub read_matrix;


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
$hr1 = readNCBI("/g/bork3/home/hildebra/DB/NCBI/ncbi_tax_table_synonyms_2014-01-23.txt");
my %NTax = %{$hr1};
$hr1 = read_matrix("$GCd/FMG.subset.mat");
my %FMGmatrix= %{$hr1};

#annotate against DB using lambda

if (0){ #too general
	my $tar = "$GCd/compl.incompl.95.fna"; my $DB = "$motuDir/263MetaRef10MGv9.cal.v2.nr.padded.fna";	my $taxblastf = "$GCd/compl.incompl.95.motuAss.tmp.m8";
	lambdaBl($tar,$DB,$taxblastf);
}

my %specIid;my %specItax;
open I,"<$SpecID/progenomes.specIv2_2";
while (<I>){next if (m/^#/);chomp; my @xx = split /\t/;$specIid{$xx[1]} = $xx[0];$xx[1]=~m/^(\d+)\./; $specItax{$xx[0]}=$1;}
close I;

#assign each COG separately
my @catsPre = split/\n/,`cat $GCd/FMG.subset.cats`; my %cats;
system "mkdir -p $MGdir" unless (-d $MGdir);
my %FMGcutoffs = (COG0012=>94.8,COG0016=>95.8,COG0018=>94.2,COG0172=>94.4,COG0215=>95.4,COG0495=>96.4,COG0525=>95.3,COG0533=>93.1,COG0541=>96.1,
COG0552=>94.5,COG0048=>98.4,COG0049=>98.7,COG0052=>97.2,COG0080=>98.6,COG0081=>98,COG0085=>97,COG0087=>99,COG0088=>99,COG0090=>98.8,COG0091=>99,
COG0092=>99,COG0093=>99,COG0094=>99,COG0096=>98.6,COG0097=>98.4,COG0098=>98.7,COG0099=>98.9,COG0100=>99,COG0102=>99,COG0103=>98.4,
COG0124=>94.5,COG0184=>98.2,COG0185=>99,COG0186=>99,COG0197=>99,COG0200=>98.4,COG0201=>97.2,COG0202=>98.4,COG0256=>99,COG0522=>98.6);
my %gene2specI; my %specItaxname; my %SpecIgenes;
foreach (@catsPre){
	my %Q2S;
	my %specIcnt;
	my @spl  = split /\t/;
	#$cats{$spl[0]} = $spl[2];
	my @genes = split(/,/,$spl[2]);
	my $COG = $spl[0];
	my $samBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";
	my $cmd = "$samBin faidx $GCd/compl.incompl.95.fna ". join (" ", @genes ) . " > $MGdir/$COG.fna\n";
	system $cmd unless (-e "$MGdir/$COG.fna");
	$cmd = "$samBin faidx $GCd/compl.incompl.95.prot.faa ". join (" ", @genes ) . " > $MGdir/$COG.faa";
	system $cmd unless (-e "$MGdir/$COG.faa");	
	lambdaBl("$MGdir/$COG.fna","$SpecID/$COG.rep.fna","$MGdir/tax/$COG.tmp.m8");
	die "no cutoff for $COG\n" unless (exists $FMGcutoffs{$COG});
	my $reqID = $FMGcutoffs{$spl[0]};
	open I,"<$MGdir/tax/$spl[0].tmp.m8" or die $!;
	while (my $line = <I>){
		chomp $line; @spl = split /\t/,$line;
		next if ($spl[2]<($reqID*1) || $spl[3] < ($spl[13]*0.9));
		#next if (exists $Q2S{$spl[0]});
		$spl[1] =~ m/^(.*\..*)\./;
		die "can't find specI $1\n" unless (exists $specIid{$1});
		#get specI assignments
		my $speci= $specIid{$1};
		if (exists($Q2S{$spl[0]} )){
			if ($Q2S{$spl[0]} =~ m/$speci/){next;}
			$Q2S{$spl[0]} .= ",".$speci; #print "m";
		} else {
			$Q2S{$spl[0]} = $speci;
		}
		$specIcnt{$speci}++;
		if (exists($SpecIgenes{$speci}{$COG} )){
			$SpecIgenes{$speci}{$COG} .= ",".$spl[0];
		} else {
			$SpecIgenes{$speci}{$COG} = $spl[0];
		}
		
	}
	close I;
	#print "Found ".keys(%Q2S)." assignments (".@genes.")\n";
	
	foreach my $k (keys %Q2S){
		last; #deactivated, need to distribute genes rather by correlations to mean
		my @spl = split /,/,$Q2S{$k};
		my $finalSpeci = $spl[0];
		my $set=0;
		#routine to distribute genes as good as possible between specIs
		for (my $i=(@spl-1);$i>=0;$i--){
		last;
			my $s = $spl[$i];
			if ( ($specIcnt{$s}>1 && $i!=0 ) || $set){#delete entries where too many specIs are available
				$specIcnt{$s}--;
			} else {
				$specIcnt{$s}--;
				$finalSpeci = $spl[$i]; $set=1;
			}
		}
		
		#$gene2specI{$k} = $finalSpeci;
	}
}
#die;
print "Done initial blast\n$MGdir\n";
my %meanFMGcors = getCorrs(\%FMGmatrix,\%SpecIgenes);
calculate_spearman_correlation($FMGmatrix{},$FMGmatrix{})



#write specI assignments for markerG
my %specIcnts;
open O,">$MGdir/gene2specI.txt";
foreach my $k (keys %gene2specI){
	print O "$k\t$gene2specI{$k}\n";
	$specIcnts{$gene2specI{$k}}++;
}
close O;

my %histo;
for my $k (sort {$specIcnts{$a} <=> $specIcnts{$b}} keys %specIcnts) {
	$histo{int $specIcnts{$k}/10}++;
    print "$k $specIcnts{$k} $NTax{$specItax{$k}}\n" ;#if ($specIcnts{$k}>=40);   # bbb c aaaa
}
foreach (sort (keys %histo)){
	print "$_\t$histo{$_}\n";
}

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

sub readNCBI($){
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






# subroutine definitions

sub convert_values_to_ranks {
	my $values = shift;

	# code below is slightly unintuitive, but we have to compute
	# average rank in cases where there are ties.
	my $idx = 0;
	my %sorted; # $sorted{$val} = [ $rank1, $rank2, $rank3 ];
	foreach my $val ( sort { $a <=> $b } @$values ) {
		push @{$sorted{$val}}, $idx;
		$idx++
	}

	# compute the average rank for a given value
	my %average_ranks =
		map { $_ => List::Util::sum(@{$sorted{$_}}) / scalar(@{$sorted{$_}}) }
		keys %sorted;

	# encode the values using average rank of the value in the list
	my @ranks = map { $average_ranks{$_} } @$values;

	if ( wantarray ) { return @ranks; }
	return \@ranks;
}

sub calculate_spearman_correlation {
	my $n1 = shift;
	my $n2 = shift;

	if ( scalar(@$n1) != scalar(@$n2) ) {
		die "Error: spearman correlation given two lists of unequal size!\n";
	}

	my $ranked_n1 = convert_values_to_ranks($n1);
	my $ranked_n2 = convert_values_to_ranks($n2);

	my $sum_diff_squared = 0;
	foreach my $idx ( 0 .. scalar(@$ranked_n1)-1 ) {
		my $diff = $ranked_n1->[$idx] - $ranked_n2->[$idx];
		$sum_diff_squared += $diff * $diff;
	}

	my $N   = scalar(@$ranked_n1);
	my $rho = 1 - ( 6 * $sum_diff_squared / ($N * ($N*$N-1)) );
	return $rho;
}

sub read_matrix($){
	my ($mF) = @_;
	my %oM;
	open I,"<$mF" or die "Can't open $mF\n";
	my $cnt=0;
	while (<I>){
		chomp; $cnt++;
		my @row = split /\t/;
		if ($cnt==1){
			shift @row;
			$oM{header} = \@row;
		} else {
			my $ID = shift @row;
			$oM{$ID} = \@row;
		}
	}
	close I;
	return \%oM;
}

#check which mean abundance profile single MG have, and how multi MGs correlate to this 
#then selects best correlating MG to be "the one" that just fits
sub getCorrs(){
	
	foreach my $k(keys %SpecIgenes){
		my @tarGenes; my $gcnt =0;
		foreach my $c(keys %{$SpecIgenes{$k}}){
			if ($SpecIgenes{$k}{$c} !~ m/,/){#single copy, use this gene
				if ($gcnt == 0){
					@tarGenes = @{$FMGmatrix{ $SpecIgenes{$k}{$c} }};
					die "@tarGenes\n".@tarGenes."\n";
				} else {
					for (my $j=0;$j<@tarGenes;$j++){
						$tarGenes[$j] += ${$FMGmatrix{ $SpecIgenes{$k}{$c} }}[$j];
					}
				}
				$gcnt++;
			}
		}
		
		#doesn't need norm, since we do spearman correlation
		#but now corr and see which genes just fit best of the multi choices..
		
		
	}
}











