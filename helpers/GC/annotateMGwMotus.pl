#!/usr/bin/env perl
#annotates specI's in the dataset
#also creates specI abundance tables
#  ./annotateMGwMotus.pl /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR 12


use warnings;
use strict;

use Mods::GenoMetaAss qw(readClstrRev median);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use List::Util;


sub readMotuTax;
sub readGene2mlinkage;
sub lambdaBl;
sub readNCBI;sub read_speci_tax;

sub calculate_spearman_correlation;
sub read_matrix; sub getCorrs;sub passBlast;
sub add2geneList; sub rm4geneList;
sub sanityCheckCorr;
sub specImatrix;
sub createAreadSpecItax;

my $SpecID="/g/bork3/home/hildebra/DB/MarkerG/specI/"; my $freeze11=1;
#my $SpecID="/g/bork3/home/hildebra/DB/MarkerG/specI_2017";my $freeze11=0;
#progenomes.specIv2_2
my $globalCorrThreshold = 0.8; # determines cutoff, when still to accept correlating genes into specI
my $reblast=0;#do blast again?

my $rarBin = getProgPaths("rare");#"/g/bork5/hildebra/dev/C++/rare/rare";
my $lambdaBin = getProgPaths("lambda");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda";
my $lambdaIdxBin = $lambdaBin."_indexer";#getProgPaths("");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda_indexer";
my $samBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";

if (@ARGV == 0){
	die "Not enough input args: use ./annotateMGwMotus.pl [path to GC] [# Cores]\n";
}

my $GCd = $ARGV[0];
my $BlastCores = $ARGV[1];
my $MGdir = "$GCd/FMG/";
system "mkdir -p $MGdir/tax" unless (-d "$MGdir/tax");

my $motuDir = "";#"/g/bork3/home/hildebra/DB/MarkerG/mOTU";
#load motu DBs...
#my ($hr1,$hr2) = readMotuTax("$motuDir/mOTU.v1.1.padded.motu.linkage.map");
#my %LG2motu = %{$hr1}; my %motu2tax = %{$hr2};
#$hr1 = readGene2mlinkage("$motuDir/mOTU.v1.1.padded.motu.map");
#my %gene2LG = %{$hr1};
#$hr1 = readNCBI("/g/bork3/home/hildebra/DB/NCBI/ncbi_tax_table_synonyms_2014-01-23.txt");
#my %NTax = %{$hr1};

my $hr1 = read_matrix("$GCd/FMG.subset.mat");
my %FMGmatrix= %{$hr1};
#$hr1 = read_speci_tax("$SpecID/specI.tax");
#my %specItax = %{$hr1};
#annotate against DB using lambda

if (0){ #too general
	my $tar = "$GCd/compl.incompl.95.fna"; my $DB = "$motuDir/263MetaRef10MGv9.cal.v2.nr.padded.fna";	my $taxblastf = "$GCd/compl.incompl.95.motuAss.tmp.m8";
	lambdaBl($tar,$DB,$taxblastf);
}

my %specIid;my %specItax;
open I,"<$SpecID/progenomes.specIv2_2";
while (<I>){next if (m/^#/);chomp; my @xx = split /\t/;$specIid{$xx[1]} = $xx[0];$xx[1]=~m/^(\d+)\./; $specItax{$xx[0]}=$1;}
close I;
my %specItaxM; #real matrix with tax levels
#die "$specItax{specI_v2_Cluster1309}\n";
my $specIfullTax = createAreadSpecItax(\%specItax,"$SpecID/specI.tax2");
my $xtrLab= "";$xtrLab= ".rep" if ($freeze11);

#assign each COG separately
my @catsPre = split/\n/,`cat $GCd/FMG.subset.cats`;
system "mkdir -p $MGdir" unless (-d $MGdir);
my %FMGcutoffs = (COG0012=>94.8,COG0016=>95.8,COG0018=>94.2,COG0172=>94.4,COG0215=>95.4,COG0495=>96.4,COG0525=>95.3,COG0533=>93.1,COG0541=>96.1,
COG0552=>94.5,COG0048=>98.4,COG0049=>98.7,COG0052=>97.2,COG0080=>98.6,COG0081=>98,COG0085=>97,COG0087=>99,COG0088=>99,COG0090=>98.8,COG0091=>99,
COG0092=>99,COG0093=>99,COG0094=>99,COG0096=>98.6,COG0097=>98.4,COG0098=>98.7,COG0099=>98.9,COG0100=>99,COG0102=>99,COG0103=>98.4,
COG0124=>94.5,COG0184=>98.2,COG0185=>99,COG0186=>99,COG0197=>99,COG0200=>98.4,COG0201=>97.2,COG0202=>98.4,COG0256=>99,COG0522=>98.6);
my %specItaxname; my %SpecIgenes; 
my %COGDBLass; my %COGass;
my %Q2S;
foreach (@catsPre){
	my %specIcnt;
	my @spl  = split /\t/;
	#$cats{$spl[0]} = $spl[2];
	my @genes = split(/,/,$spl[2]);
	my $COG = $spl[0];
	
	#actual blast (heavy)
	my $cmd = "$samBin faidx $GCd/compl.incompl.95.fna ". join (" ", @genes ) . " > $MGdir/$COG.fna\n";
	system $cmd unless (-e "$MGdir/$COG.fna");
	$cmd = "$samBin faidx $GCd/compl.incompl.95.prot.faa ". join (" ", @genes ) . " > $MGdir/$COG.faa";
	system $cmd unless (-e "$MGdir/$COG.faa");	
	
	#print "$MGdir/tax/$COG${xtrLab}.tmp.m8\n";
	lambdaBl("$MGdir/$COG.fna","$SpecID/$COG${xtrLab}.fna","$MGdir/tax/$COG${xtrLab}.tmp.m8"); #.rep
	#next;
	
	die "no cutoff for $COG\n" unless (exists $FMGcutoffs{$COG});
	my $reqID = $FMGcutoffs{$spl[0]};
	open I,"<$MGdir/tax/$spl[0]${xtrLab}.tmp.m8" or die $!;
	while (my $line = <I>){
		chomp $line; @spl = split /\t/,$line;
		next if (passBlast(\@spl,$reqID)==0);
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
		
		if (exists($SpecIgenes{$speci}{$COG} )){
			$SpecIgenes{$speci}{$COG} .= ",".$spl[0];
			$COGDBLass{$COG}++;
		} else {
			$SpecIgenes{$speci}{$COG} = $spl[0];
			$COGass{$COG}++;
		}		
	}
	close I;
}
#die;
#print "Done initial blast\n$MGdir\n";
#foreach (sort {$COGass{$a} cmp $COGass{$b}} keys %COGass){print "$_ $FMGcutoffs{$_}:  $COGass{$_}($COGDBLass{$_})\n";}

my %gene2specI; my %specIprofiles; my %SpecIgenes2;
#sort out best multi hit by correlation analysis
getCorrs(\%FMGmatrix,\%SpecIgenes);

#check that all corrs here check out well..
sanityCheckCorr();

#second round.. go through blast again and reassign everything, that hasn't got a hit yet

undef %Q2S; my $xtraEntry=0;
@catsPre = split/\n/,`cat $GCd/FMG.subset.cats`;
foreach (@catsPre){
#last;
	my %specIcnt;
	my @spl  = split /\t/;
	#$cats{$spl[0]} = $spl[2];
	my @genes = split(/,/,$spl[2]);
	my $COG = $spl[0];
	
	my $reqID = $FMGcutoffs{$spl[0]};
	#secondary scanning of assignments...
	open I,"<$MGdir/tax/$spl[0]${xtrLab}.tmp.m8" or die "Can't open $MGdir/tax/$spl[0]${xtrLab}.tmp.m8\n";
	while (my $line = <I>){
		chomp $line; @spl = split /\t/,$line;
		my $gid = $spl[0];
		#1 blast good enough
		next if (passBlast(\@spl,$reqID)==0);
		
		#2 gene not already assigned in first high confidence pass
		next if (exists $gene2specI{$gid});
		$spl[1] =~ m/^(.*\..*)\./;
		die "can't find specI $1\n" unless (exists $specIid{$1});
		#get specI assignments
		my $speci= $specIid{$1};
		
		#3 this MG has not been assigned in the high confidence initial assignments
		if (exists($SpecIgenes2{$speci}{$COG} )){next;}
		
		#4 correlate to species core, to make sure the gene kind of makes sense..
		if (exists($specIprofiles{$speci})){
			my $corr = correlation($specIprofiles{$speci},$FMGmatrix{ $gid }) ;
			next if ($corr < $globalCorrThreshold);
			#print $corr."\t";
		}


		if (exists($Q2S{$gid} )){
			if ($Q2S{$gid} =~ m/$speci/){next;}
			$Q2S{$gid} .= ",".$speci; #print "m";
		} else {
			$Q2S{$gid} = $speci;
		}
		$specIcnt{$speci}++;
		$xtraEntry++;
		#and block gene slot
		add2geneList($speci,$COG,$gid);
		#$gene2specI{$spl[0]} = $speci;
	}
	close I;
	#print "Found ".keys(%Q2S)." assignments (".@genes.")\n";
	
	foreach my $k (keys %Q2S){
		last;
		my @spl = split /,/,$Q2S{$k};
		my $finalSpeci = $spl[0];
		my $set=0;
		#routine to distribute genes as good as possible between specIs
		for (my $i=(@spl-1);$i>=0;$i--){
			my $s = $spl[$i];
			if ( ($specIcnt{$s}>1 && $i!=0 ) || $set){#delete entries where too many specIs are available
				$specIcnt{$s}--;
			} else {
				$specIcnt{$s}--;
				$finalSpeci = $spl[$i]; $set=1;
			}
		}
		
		$gene2specI{$k} = $finalSpeci;
	}
}

#die;
sanityCheckCorr();

my %specIcnts;
foreach my $k (keys %gene2specI){$specIcnts{$gene2specI{$k}}++;}

#write specI assignments for markerG
open O,">$MGdir/gene2specI.txt";
foreach my $k (keys %gene2specI){
	print O "$k\t$gene2specI{$k}\n";
}
close O;

#create abundance profile
specImatrix("$MGdir/specI.mat",$specIfullTax);





my %histo;
for my $k (sort {$specIcnts{$a} <=> $specIcnts{$b}} keys %specIcnts) {
	$histo{int $specIcnts{$k}/10}++;
   # print "$k $specIcnts{$k} $NTax{$specItax{$k}}\n" ;#if ($specIcnts{$k}>=40);   # bbb c aaaa
}
foreach (sort (keys %histo)){
	print "$_\t$histo{$_}\n";
}
print "Extra: $xtraEntry\n";
exit(0);
#TOGO




#####################################################


sub rebase($){ #calculates the profile for each SI based on the 40 marker genes
	my ($hr1)=@_; my %specIs = %{$hr1};
	my $medianC=0;
	foreach my $sp (keys %specIs){
	#last;
		next if ($sp =~ m/,/); #from where do these come?
		my @tar;my $MGn=0;
		foreach my $gid (@{$specIs{$sp}}){
			$MGn++;
			#now get genes from matrix
			if ($medianC){
				for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){push(@{$tar[$j]}, ${$FMGmatrix{ $gid }}[$j]);}
			} else {
				for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$tar[$j] +=  ${$FMGmatrix{ $gid }}[$j];}
			}
		}
		if ($medianC){
			my @tar2;for (my $j=0;$j<scalar(@tar);$j++){$tar2[$j] = median(@{$tar[$j]});} $specIprofiles{$sp} = \@tar2;
		} else {
			for (my $j=0;$j<scalar(@tar);$j++){$tar[$j] /= $MGn;}$specIprofiles{$sp} = \@tar;
		}
		# if ($MGn>2);
	}
}

sub sanityCheckCorr(){
	#final sanity check, that marker genes are correlating
	my $wrongGene=0; my $corrGene=0;
	my %specIGcorrs;my %specIcnts; my %specIs;
	foreach my $k (keys %gene2specI){$specIcnts{$gene2specI{$k}}++;push(@{$specIs{$gene2specI{$k}}},$k);}
	
	#rebase..
	rebase(\%specIs);
	
	
	#actual correlation check
	foreach my $k (keys %gene2specI){
		#$specIGset{$gene2specI{$k}}
		my $sIS = $gene2specI{$k};
		foreach my $sI (split(/,/,$sIS)){ #in case several SIs have been assigned to gene..
			next unless ($specIcnts{$sI} > 2);
			next unless (exists($specIprofiles{$sI}));
			#if (nonZero($specIprofiles{$sI}) < 3){ next;}
			my $corr = correlation($specIprofiles{$sI},$FMGmatrix{ $k }) ;
			#print $corr." ";
			if ($corr < $globalCorrThreshold){$wrongGene++;rm4geneList($sI,$k);$specIcnts{$sI}--;}
			$corrGene++;
			push (@{$specIGcorrs{$sI}},$corr);
		}
	}
	print "correct(corr): $corrGene, rem(corr): $wrongGene. \n";
	
	undef %specIs;
	foreach my $k (keys %gene2specI){push(@{$specIs{$gene2specI{$k}}},$k);}
	rebase(\%specIs);
}

sub passBlast($ $){
	my $spl = shift;my $reqID = shift;
	my $lengthGood = 0;
	if (@{$spl} < 12 || !defined($spl->[13])|| !defined($spl->[12])){die "not enough entries:\n@{$spl}\n";}
	$lengthGood=1 if ($spl->[3] >= ($spl->[13]*0.1) && $spl->[3] >= ($spl->[12]*0.7)  ); #subject query
	if ($spl->[2] < ($reqID*0.985) || !$lengthGood)  {return 0;}
	return 1;
}
sub lambdaBl($ $ $){
	my ($tar,$DB, $taxblastf) = @_;
	my $cmd="";
	
	if (!-d $DB.".lambda/"  ) {
		print "Building LAMBDA index anew (may take up to an hour)..\n";
		my $cmdIdx = "$lambdaIdxBin -p blastn -t ".int($BlastCores)." -d $DB";
		if (system($cmdIdx)){die ("Lamdba ref DB build failed\n$cmdIdx\n");}
	}
	$cmd .= "$lambdaBin -t $BlastCores -id 93  -nm 100 -p blastn -e 1e-40 -q $tar -oc \"std qlen slen\" -i $DB.lambda -o $taxblastf\n";
	#die $cmd."\n";
	system $cmd if (!-e $taxblastf || $reblast);
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
sub read_speci_tax($){
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
#correlation routines
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

	if ( scalar(@{$n1}) != scalar(@{$n2}) ) {
		die "Error: spearman correlation given two lists of unequal size!\n@{$n2}\n@{$n1}\n";
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
 
sub mean {
   my ($x)=@_;
   my $num = scalar(@{$x}) - 1;
   my $sum_x = '0';
   for (my $i = 1; $i < scalar(@{$x}); ++$i){
      $sum_x += $x->[$i];
   }
   my $mu_x = $sum_x / $num;
   return($mu_x);
}
 
### ss = sum of squared deviations to the mean
sub ss {
   my ($x,$mean_x,$y,$mean_y)=@_;
   my $sum = '0';
   for (my $i=1;$i<scalar(@{$x});++$i){
     $sum += ($x->[$i]-$mean_x)*($y->[$i]-$mean_y);
   }
   return $sum;
}
 sub correl {
   my($ssxx,$ssyy,$ssxy)=@_;
   if ($ssyy==0 || $ssxy==0){return 0;}
   my $sign=$ssxy/abs($ssxy);
   my $correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
   return $correl;
}
sub correlation {
   my ($x,$y) = @_;
   my ($mean_x) = mean($x);
   my $mean_y = mean($y);
   my $ssxx=ss($x,$mean_x,$x,$mean_y);
   my $ssyy=ss($y,$mean_x,$y,$mean_y);
   my $ssxy=ss($x,$mean_x,$y,$mean_y);
   my $correl=correl($ssxx,$ssyy,$ssxy);
   my $xcorrel=sprintf("%.4f",$correl);
   return($xcorrel);
 
}
 


sub read_matrix($){
	my ($mF) = @_;
	print "Reading matrix\n";
	my %oM;
	open I,"<$mF" or die "Can't open $mF\n";
	my $cnt=0;
	while (<I>){
		chomp; $cnt++;
		my @row = split /\t/;
		my $ID = shift @row;
		if ($cnt==1){
			$oM{header} = \@row;
		} else {
			#for (my $i=0;$i<@row;$i++){if ($row[$i] < 2){$row[$i]=0;} } 
			#for (my $i=0;$i<@row;$i++){$row[$i] = sqrt ($row[$i]);} 
			$oM{$ID} = \@row;
		}
	}
	close I;
	return \%oM;
}
sub specImatrix($$){
	my ($oF,$hr) = @_;
	my %sTax = %{$hr};
	#print "@{$sTax{specI_v2_Cluster34}}\n";
	open O,">$oF";
	my @bkgrnd; my @dblCh;
	foreach my $gid (keys %FMGmatrix){
		next if ($gid eq "header");
		for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$dblCh[$j] +=  ${$FMGmatrix{ $gid }}[$j]/40;}
		next if (exists ($gene2specI{$gid}));
		for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$bkgrnd[$j] +=  ${$FMGmatrix{ $gid }}[$j] / 40;}
	}
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";
	
	print O "SpecI\t".join ("\t",@{$FMGmatrix{ header }})."\n";
	#print O "SUM\t\t".join ("\t",@dblCh)."\n";
	print O "?\t".join ("\t",@bkgrnd)."\n";
	foreach my $si (keys %specIprofiles){
		if (exists($specItax{$si})){
#			print O "$si\t$specItax{$si}\t@{$sTax{$si}}\t".join ("\t",@{$specIprofiles{ $si }})."\n";
			print O "$si\t".join ("\t",@{$specIprofiles{ $si }})."\n";
		} else {
			print "Can't find $si\n";
		}
		
		#for (my $j=0;$j<scalar(@{$specIprofiles{ $gid }});$j++){$bkgrnd[$j] +=  ${$specIprofiles{ $gid }}[$j];}
	}
	close O;
	#calculating higher level abundance matrix
	my $oFx = $oF; $oFx =~ s/\.[^\.]*$//;
	my @taxLs = ("superkingdom","phylum","class","order","family","genus","species");
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";
	for (my $t=0;$t<@taxLs;$t++){
		#sum up to hi lvl
		my %thisMap;
		foreach my $si (keys %specIprofiles){
			die "doesnt exist: $si\n" if (!exists($sTax{$si}));
			print "ERR: $si    @{$sTax{$si}}\n" if (@{$sTax{$si}} <= $t);
			my $clvl = join (";",@{$sTax{$si}}[0 .. $t]);
			
			#if ($clvl eq "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia albertii"){
				#print "$si\n@{$specIprofiles{ $si }}\n";
			#}
			
			if (exists($thisMap{$clvl})){
				for (my $kl=0;$kl<scalar(@{$specIprofiles{ $si }});$kl++){
					${$thisMap{$clvl}}[$kl] += ${$specIprofiles{ $si }}[$kl];
				}
			} else {
				$thisMap{$clvl} = [@{$specIprofiles{ $si }}];
			}
		}
		open O,">$oFx.$taxLs[$t]";
		print O "$taxLs[$t]\t".join ("\t",@{$FMGmatrix{ header }})."\n";
		print O "?\t".join ("\t",@bkgrnd)."\n";
		foreach my $kk (keys %thisMap){
			print O "$kk\t".join("\t",@{$thisMap{$kk}}) . "\n";
		}
		close O;
	}
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";

}
sub nonZero($){
	my ($ar) = @_;
	my @a = @{$ar};
	my $nc =0 ;
	foreach my $x (@a){ $nc ++ if ($x > 0);}
	return $nc;
}

sub add2geneList($ $ $){
	my ($k,$c,$sg) = @_;
	my $ret=0;
	if (exists($gene2specI{$sg})){
		$gene2specI{$sg} .= ",".$k; $ret=1;#print "should not happen\n";#$gene2specI{$spl[$i]}   $k\n";
	} else {$gene2specI{$sg} = $k; } 
	$SpecIgenes2{$k}{$c}=$sg;#set mark to block this MG in this specI...
	return $ret;
}

sub createAreadSpecItax{
	my ($hr1,$file) = @_;#\$specIid,"$SpecID/specI.tax2");
	my %sNTID = %{$hr1};
	my %ret; my $cnt = -1;
	if (-e $file){
		open I,"<$file";
		while (my $l = <I>){
			$cnt++;
			chomp $l; my @spl = split /\t/,$l;
			#if ($cnt == 0){}
			my $id = shift @spl;
			$ret{$id} = \@spl;
		}
		close I;
	}
	#print "$ret{specI_v2_Cluster1309}\n";
	my @taxids; my @spids;
	foreach my $k (keys %sNTID){
		next if (exists($ret{$k}));
		#jaimes ete prog
		#print "$k\n";
		push (@taxids,$sNTID{$k});
		push (@spids,$k);

	}
	#die "python /g/bork3/home/hildebra/dev/python/get_ranks.py ".join(" ",@taxids);
	if (@taxids>0){
		print "Detecting ".@taxids." new taxids\n";
		my $cmd = "python /g/bork3/home/hildebra/dev/python/get_ranks.py ".join(" ",@taxids);
		my $tret= `$cmd`;
		my @newT = split /\n/,$tret;
		if (@newT > 0){
			open O,">>$file";
			for (my $i=0;$i<@newT;$i++){
				my @spl = split /\t/,$newT[$i]; shift @spl;
				print O "$spids[$i]\t".join("\t",@spl)."\n";
			}
			close O;
		}
	}
	#die "@{$ret{specI_v2_Cluster34}}\n";
	return (\%ret);
}

sub rm4geneList($ $){
	my ($k,$sg) = @_;
	#die "ASDA";
	if ($gene2specI{$sg} !~ m/,/){delete $gene2specI{$sg};
	} else {#more complicated, remove from comma list
		my @spl = split /,/,$gene2specI{$sg};
		#print "@spl\n";
		@spl = grep (!/$k/i, @spl);
		#die "TODO: @spl  XX $k";
		$gene2specI{$sg} = join ",", @spl;
	}
	
	return unless (exists($SpecIgenes2{$k}));
	#print "A";
	#$c is unknown..
	my $entrFound=0;
	foreach my $c (keys %{$SpecIgenes2{$k}}){
		#print "L";
		if ($SpecIgenes2{$k}{$c} eq $sg){
			delete $SpecIgenes2{$k}{$c}; $entrFound=1;last;
		}
	}
	#die "not deleted from SpecIgenes2 $k $sg\n" unless ($entrFound==1);
	#return $ret;
}


#check which mean abundance profile single MG have, and how multi MGs correlate to this 
#then selects best correlating MG to be "the one" that just fits
sub getCorrs(){
	my $dblAssi=0;
	my $dblA=0;my $singlA=0;my $singlMultA=0;my $skippedSIs=0; my $newAssigns=0;
	foreach my $k(keys %SpecIgenes){
		my @tarGenes; my $gcnt =0;
		foreach my $c(keys %{$SpecIgenes{$k}}){
			if ($SpecIgenes{$k}{$c} =~ m/,/){next;}#single copy, use this gene
			my $gid = $SpecIgenes{$k}{$c};
			if ($Q2S{$gid} =~ m/,/){$dblA++;;next;}
			$singlA++; 
			#if ($gcnt == 0){@tarGenes = @{$FMGmatrix{ $gid }};
			#} else {
			for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$tarGenes[$j] += ${$FMGmatrix{ $gid }}[$j];}
			#} 
			$gcnt++;				
			#$gene2specI{$gid} = $k;   #can still be wrong assignment.. 
			$dblAssi++ if (add2geneList($k,$c,$gid));
			#if (exists($gene2specI{$gid})){
			#	$gene2specI{$gid} .= ",".$k; $dblAssi++;#print "should not happen\n";#$gene2specI{$spl[$i]}   $k\n";
			#} else {$gene2specI{$gid} = $k; } 
		}
		
		if ($gcnt < 5 && scalar(keys %{$SpecIgenes{$k}}) >= 5){ #no use, calc everything together up
			#These genes will later be checked again for the correlation to the grand mean
			foreach my $c(keys %{$SpecIgenes{$k}}){
			#last;
				if ($SpecIgenes{$k}{$c} !~ m/,/){next;}
				my @spl = split /,/,$SpecIgenes{$k}{$c};
				my @subTarG; my $sgcnt=0;
				foreach my $gid (@spl){#single copy, use this gene
					if ($Q2S{$gid} =~ m/,/){next;}#don't want multi species assignments
					#if ($sgcnt == 0){@subTarG = @{$FMGmatrix{ $gid }}; 
					#} else {
						for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$subTarG[$j] += ${$FMGmatrix{ $gid }}[$j];}
					#}
					$sgcnt++;
				} 
				if ($sgcnt>0){
					$singlMultA++;
					$gcnt++ ;
					for (my $j=0;$j<@subTarG;$j++){$tarGenes[$j] += $subTarG[$j] / $sgcnt;}
				}
			}

		}
		if ($gcnt ==0){$skippedSIs++;next;}	
		if (nonZero(\@tarGenes) < 3){ next;}#print "XX" ;
		
		#doesn't need norm, since we do spearman correlation
		#but now corr and see which genes just fit best of the multi choices..
		#print "@tarGenes\n";
		foreach my $c(keys %{$SpecIgenes{$k}}){
			if ($SpecIgenes{$k}{$c} =~ m/,/){#only look at multi assigned genes
				#die "$SpecIgenes{$k}{$c}\n";
				my @spl = split /,/,$SpecIgenes{$k}{$c};
				my @subCors; my $max=0;
				foreach my $sg (@spl){
					die "$sg doesn't exist in gene list \n" if (!exists($FMGmatrix{ $sg }));
					
					#my $corr = calculate_spearman_correlation(\@tarGenes,$FMGmatrix{ $sg }) ;
					#print "@{$FMGmatrix{ $sg }}\n";
					my $corr = correlation(\@tarGenes,$FMGmatrix{ $sg }) ;
					if ($corr > $max){$max = $corr;}
					push (@subCors,  $corr  );
				}
				
				
				if ($max < $globalCorrThreshold){next;}
				$newAssigns++;
				for (my $i=0;$i<@spl;$i++){
					unless ($subCors[$i]>$max-0.03){next;}
					my $sg = $spl[$i];
					#this assignment is what I need, now I know this gene is blocked for assignment to other SpecI's
					$dblAssi++ if (add2geneList($k,$c,$sg));

					
					for (my $j=0;$j<@tarGenes;$j++){$tarGenes[$j] += ${$FMGmatrix{ $sg }}[$j] ;}
				}

				#print  "@subCors\n";
			}
		}
		
		#norm vector
		for (my $j=0;$j<@tarGenes;$j++){$tarGenes[$j] /= $gcnt;}
		#and save the final specI profile..
		$specIprofiles{$k} = \@tarGenes;
	}
	print "double assignment $dblAssi; assigned: $newAssigns   Stats in Run: $dblA $singlA $singlMultA $skippedSIs\n";
}











