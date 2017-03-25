#!/usr/bin/perl
#builds from a contig file several stats to seperate contigs into single species
#$subparts =~ a=gene abundance g=GC e=essential k=kmer m=metabat s=microsats
#./separateContigs.pl /g/scb/bork/hildebra/SNP/GNMass/alien-11-374-0/
use warnings;
use strict;
#use Scalar::Util qw(looks_like_number);


sub readFasta;
sub readGFF;
sub geneAbundance; sub runMaxBin;
sub findMicrSat;
use Mods::GenoMetaAss qw(systemW is_integer reverse_complement_IUPAC);
use Mods::IO_Tamoc_progs qw(getProgPaths );

my $inD = $ARGV[0];
my $assD = $ARGV[1];
my $subparts = $ARGV[2];
my $readLength = $ARGV[3];


#my $rdCovBin = "/g/bork3/home/hildebra/dev/C++/rdCov/./rdCover";
my $rdCovBin =getProgPaths("readCov");
my $growthBin =getProgPaths("growthP");
my $FMGd = getProgPaths("FMGdir");#"/g/bork5/hildebra/bin/fetchMG/";
my $compoundBinningScr = getProgPaths("cmpBinScr");#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/compoundBinning.pl";
my $mrepsB = getProgPaths("mreps");#"/g/bork3/home/hildebra/bin/mreps/./mreps";


my $inScaffs = "$assD/scaffolds.fasta.filt";
my $proteins = "$assD/genePred/proteins.shrtHD.faa";
my $genesNT = "$assD/genePred/genes.shrtHD.fna";
my $genesAA = "$assD/genePred/proteins.shrtHD.faa";
my $rawINPUTf = $inD."input.txt";
my $outDab = $inD."/assemblies/metag/ContigStats/";
my $outD = $assD."/ContigStats/";
my $cmd = "";
my $cleanUp = 0;#1=do overwrite
if ($cleanUp){
	system "rm -r $outD";
}
system "mkdir -p $outD";
system "mkdir -p $outDab";

#figure out mapping names
my $SmplNm = `cat $inD/mapping/done.sto`;
$SmplNm =~ s/-smd.bam\n?//;
#my $inBAM = $inD."mapping/$SmplNm-smd.bam";
my $coverage = $inD."mapping/$SmplNm-smd.bam.coverage";
#die "$coverage\n";
#system "cp $inScaffs $outD";
######################### ini calcs
#$cmd = "";
#$cmd .= "cut -f1 -d \" \" $proteins > $proteins.shrtHD.faa\n" unless (-e "$proteins.shrtHD.faa");
#$cmd .= "cut -f1 -d \" \" $genesNT > $genesNT.shrtHD.fna\n" unless (-e "$genesNT.shrtHD.fna");
#systemW $cmd unless ($cmd eq "");
###############################   Coverage per gene  ###############################


geneAbundance($coverage) if ($subparts =~ m/a/);

###############################   GC content  ###############################
my $GCP = getProgPaths("calcGC_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/calcGC.pl";
if ((!-s "$outD/scaff.GC" || !-s "$outD/scaff.pergene.GC") && $subparts =~ m/g/){
	print "Analyzing GC content\n";
	systemW "$GCP $inScaffs $outD/scaff.GC";
	systemW "$GCP $genesNT $outD/scaff.pergene.GC";
	#systemW $cmd."\n";
} else {
	print "GC content already calculated\n";
}

###############################   essential proteins  ###############################
my $oDess = "$outD/ess100genes/";my $outDFMG = "$outD/FMG/";
if ($subparts =~ m/e/){
	print "Searching for core proteins in predicted genes..    ";
	###############################   fetchMG  ###############################
	my $FMGrwkScr = getProgPaths("FMGrwk_scr");
	system("mkdir -p $outDFMG");
	$cmd = "perl $FMGd/fetchMG.pl -m extraction -o $outDFMG -l $FMGd/lib -t 1 -d $genesNT $proteins"; # -x $FMGd/bin  -b $FMGd/lib/MG_BitScoreCutoffs.uncalibrated.txt
	$cmd .= "; $FMGrwkScr $outDFMG";
	if ( !-s "$outDFMG/COG0012.faa" && !-s "$outDFMG/COG0016.faa"){
		systemW $cmd;
	}
	system "rm  -rf $inD/assemblies/metag/genePred/*.cidx $outDFMG/temp $outDFMG/hmmResults";
	############################### essential 100 proteins ###############################
	my $hmmbin = getProgPaths("hmmsearch");#"/g/bork5/hildebra/bin/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmsearch";
	my $essDB = getProgPaths("essentialHMM");#"/g/bork5/hildebra/bin/multi-metagenome-master/R.data.generation/essential.hmm";
	my $essEukDB = "/g/bork3/home/hildebra/DB/HMMs/eukCore/eukCore.hmm"; #TODO
	
	systemW("mkdir -p $oDess");

	$cmd = "";
	$cmd .= "$hmmbin --domtblout $oDess/assembly.hmm.orfs.txt -o $oDess/assembly.hmmsout.txt --cut_tc --cpu 1 --notextw $essDB $proteins\n";
	$cmd .= "tail -n+4 $oDess/assembly.hmm.orfs.txt | sed \'s/\\s\\s*/ /g\' | sed \'s/^#.*//g\' | cut -f1,4 -d \" \" > $oDess/ess100.id.txt\n";

	if (!-s "$oDess/ess100.id.txt" && !-e "$oDess/e100split.sto"){
		print "\nDetecting 100 essential proteins\n";
		systemW $cmd ."\n";
	}
	my $protIDss = `cut -f1 -d " " $oDess/ess100.id.txt `;
	my $protClsTmp = `cut -f2 -d " " $oDess/ess100.id.txt `;
	my @protIDs = split("\n",$protIDss);
	my @protCls = split("\n",$protClsTmp);
	my $phr = readFasta("$proteins");
	my $ghr = readFasta("$genesNT");
	my %prots = %{$phr}; my %essProt;
	my %genes = %{$ghr};

	my %seen;
	#split 100 essentials into separate files for each gene class
	#my @unique = grep { ! $seen{$_}++ } @faculty;
	foreach ( grep { ! $seen{$_}++ } @protCls){
		open O,">$oDess/pe100_".$_.".faa";	close O;
		open O,">$oDess/ge100_".$_.".fna";	close O;
	}
	for (my$i=0;$i<@protIDs;$i++){
		my $id = $protIDs[$i];
		next if ($id eq "");
		if (!exists($prots{$id})){die "Can't find $id protein in file $proteins\n";}
		if (!exists($genes{$id})){die "Can't find $id protein in file $proteins\n";}
		open O,">>$oDess/pe100_".$protCls[$i].".faa" or die "Can't open $oDess/pe100_$protCls[$i].faa";
		print O ">$id\n$prots{$id}\n";
		close O;
		open O,">>$oDess/ge100_".$protCls[$i].".fna"or die "Can't open $oDess/ge100_$protCls[$i].faa";
		print O ">$id\n$genes{$id}\n";
		close O;
	}
	sleep (3);
	systemW("touch $oDess/e100split.sto;";
	#required for maxbin
	#rm $oDess/assembly.hmm*");
	#die ($protIDs[0]."\n");

	print "Done\n";
} else {
	print "No essential proteins requested\n";
}

###############################   growth prediction of sample  ###############################
if ($subparts =~ m/g/){
	systemW("cat $oDess/*.fna >$oDess/alle100.fna") unless (-e "$oDess/alle100.fna");
	$cmd = "$growthBin -f $oDess/alle100.fna -g $genesNT -c 0 -T 37 -m";
}

###############################   kmer content  ###############################
#print $cmd."\n";
if (!-s "$outD/scaff.4kmer.gz" && $subparts =~ m/k/){
	print "Analyzing kmer content\n";
	my $KCP = getProgPaths("kmerFreqs");#"perl /g/bork5/hildebra/bin/multi-metagenome-master/R.data.generation/calc.kmerfreq.pl";
	systemW "$KCP -i $inScaffs -m 100 -o $outD/scaff.4kmer";
	systemW "gzip  $outD/scaff.4kmer";
	#system "gzip $outD/scaff.4kmer";
	#systemW $cmd."\n";
}

sleep(2);
###############################   binning  ###############################
if ( !-s $inD."Binning/MaxBin/MB.summary" && $subparts =~ m/m/){
	runMaxBin();
}
if ( !-s $inD."Binning/MetaBat/MeBa.sto" &&  $subparts =~ m/m/){
	system "perl $compoundBinningScr $inD";
}

if (int(-s "$outD/microsat.txt") == 0 && $subparts =~ m/s/){
	findMicrSat($inScaffs,"$outD/microsat.txt");
}

sub findMicrSat{
	my ($scaffs,$rep) = @_;
	my $cmd = "$mrepsB -fasta -minperiod 2 -maxperiod 13 -maxsize 200 $scaffs > $rep\n";
	#die $cmd;
	print "Searching for microsattelites..";
	systemW $cmd;
	print "Done\n";
}



sub runMaxBin {
	print "Running MaxBin..\n";
	my $maxBin = getProgPaths("maxBin");#"perl /g/bork5/hildebra/bin/MaxBin-1.4.2/run_MaxBin.pl";
	my $rwkMB = getProgPaths("maxBinFrmt_scr");# "perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/maxBin_rewrk.pl";
	my $covCtgs = "$outDab/Coverage.percontig";
	my $essHMM = "$outD/ess100genes/assembly.hmm.orfs.txt";
	my $outF = "MB";
	my $outD2 = $inD."Binning/MaxBin/";

	print "Binning contigs with MaxBin\n";
	system ("mkdir -p $outD2");

	my $cmd  = "$maxBin -contig $inScaffs -out $outD2/$outF -abund $covCtgs -thread 10 -HMM $essHMM -AA $genesAA\n";
	$cmd .= "$rwkMB $outD2 $outF";
	system $cmd;
}

sub calcGeneCov($ $ $){
	my ($arr,$hr,$curChr) = @_;
	#print "geneCov\n";
	my %gff = %{$hr}; my $retStr = "";
	my @cov = @{$arr};
	my $cnt=1;
	unless (exists($gff{$curChr})){
		print "$curChr chromosome doesn't exist\n" ;
		return "";
	}
	foreach my $ge (keys %{$gff{$curChr}} ){
		my $st = $gff{$curChr}{$ge}{start}; my $en = $gff{$curChr}{$ge}{stop};
		if ($en > @cov){print "$en in $curChr ".@cov."\n";}
		my $gcov = 0; for ($st .. $en){ $gcov += $cov[$_];}
		$gcov /= ($en-$st);
		#die $gcov."\n";
		$retStr .= $curChr."_$cnt\t$gcov\n";
		$cnt++;
	}
	#die();
	return $retStr;
}

sub geneAbundance{
	my ($inF) = @_;
	#my $hr = readGFF($inD."assemblies/metag/genePred/genes.gff");
	my $outF = $inF . ".pergene"; 	my $outF2 = $inF . ".percontig";	
	my $outF3 = $inF . ".window";	my $outF4 = $inF . ".geneStats";
	my $outF5 = $inF . ".count_pergene";
	my $fileEnd2 = "";
	my $outFfin = $outDab . "Coverage.pergene$fileEnd2"; 	my $outF2fin = $outDab . "Coverage.percontig$fileEnd2";	
	my $outF3fin = $outDab . "Coverage.window$fileEnd2";	my $outF4fin = $outDab . "GeneStats.txt";	
	my $outF5fin = $outDab . "Coverage.count_pergene$fileEnd2";	
	if (-s $outFfin && -s $outF2fin&& -s $outF3fin&& -s $outF4fin && -s $outF5fin){print "Gene abundance was already calculated\n";return;}
	my $clnCmd = "";
	my $inFG = $inF.".gz";
	if (-e $inFG && !-s $inF){system "rm $inF";}
	if (-e $inFG && !-e $inF){
		systemW("gunzip -c $inFG > $inF"); $clnCmd .= "rm $inF\n";
	} elsif (-e $inF && !-e $inFG){
		$clnCmd .= "gzip $inF\n";
	}
	#open I,"<$inF" or die "Can't open $inF\n";	open O,">$outF" or die "Can't open output $outF\n";	open O2,">$outF2" or die "Can't open output $outF2\n";
	#1st delinearize	my $curChr = ""; my $cont=1; my $fresh = 1; my $expLength = 0;	my @chrCov=();	my $start = time; my $ctgCnt = 0; my $ccov = 0; #contig coverage
	#while (my $line = <I>){		chomp $line; 		my @spl = split(/\t/,$line);		my $chr = $spl[0];		if ($curChr ne $chr ){			unless  ($curChr eq ""){				#print "1";				$curChr =~ m/_length_(\d+)_/;				#if (length(@chrCov)<= $1){@chrCov = (@chrCov,(0)x($1-@chrCov+1)); }								#for (0 .. (@chrCov-1)){ $ccov += $chrCov[$_];} $ccov /= (@chrCov-1);				print O2 "$curChr\t$ccov\n"; $ccov = 0;				#print @chrCov." @chrCov $1\n";				#next if ($1 < 500); #because no gene calls on these				print O calcGeneCov(\@chrCov,$hr,$curChr);				$ctgCnt++;				if ($ctgCnt == 500){					my $duration2 = time - $start;					die "Time for 500 : ". $duration2." s\n";				}			}			#die $outF."\n";			$curChr = $chr;  $fresh = 1;			$curChr =~ m/_length_(\d+)_/;			$expLength = $1; @chrCov = ((0)x($expLength+1));			#my @test = (int($spl[1])..(int($spl[2])-1)); print "@test X $spl[3]\n";		}		#print $spl[3]."\n";		#if ($fresh){$fresh = 0;	if ($spl[1] != 0){@chrCov = (@chrCov,(0)x($spl[1]));}}
		#@chrCov = push(@chrCov,($spl[3])x($spl[2]-$spl[1]));		my $fill = $spl[3];		#$ccov += $spl[3] * ($spl[2] - $spl[1]-1);		foreach ($spl[1] .. ($spl[2]-1)){ $chrCov[$_] = $fill }		#if (!is_integer($spl[2])){die $spl[2]."  2 \n";}		#if (!is_integer($spl[1])){die $spl[1]."  1 \n";}		#if (!is_integer($spl[3])){die $spl[3]."  3 \n";}		#$chrCov[$spl[1]..$spl[2]-1] = $spl[3];		#$chrCov[11..22] = int($spl[3]);				#die @chrCov." @chrCov\n" if ($spl[3] == 6);	}
	#close I; close O; close O2;
	#my $duration = time - $start;
	#print "Gene abundances calculated in $duration s\n";
	#die ("$clnCmd");
	
	print "$rdCovBin $inF $assD/genePred/genes.gff $readLength";
	systemW "$rdCovBin $inF $assD/genePred/genes.gff $readLength";
	#$cmd .= "gzip -f $outF $outF2 $outF3\nmv $outF.gz $outFfin\n mv $outF2.gz $outF2fin\nmv $outF3.gz $outF3fin\n";
	systemW "mv $outF$fileEnd2 $outFfin\nmv $outF2$fileEnd2 $outF2fin\nmv $outF3$fileEnd2 $outF3fin\nmv $outF4$fileEnd2 $outF4fin\nmv $outF5$fileEnd2 $outF5fin\n";
	
	#die $cmd."\n";
	#systemW $cmd;
	#print "cont";
	system "echo \"$assD\" > $inD/assemblies/metag/assembly.txt";
	systemW($clnCmd);
}

sub readGFF($){
	my ($inF) = @_;
	my %gff;
	open I,"<",$inF or die "Can't find ".$inF."\n";
	my $curChrCnt = 0;
	while (my $line = <I>){
		next if ($line =~ m/^#/);
		chomp $line; 
		#die $line."\n";
		my @spl = split (/\t/,$line);
		my $chr = $spl[0];
		if (exists($gff{$chr})){
			$curChrCnt ++;
		} else {
			$curChrCnt = 0;
		}
		#die $spl[3]."\n";
		$gff{$chr}{$curChrCnt}{start} = $spl[3];
		$gff{$chr}{$curChrCnt}{stop} = $spl[4];
	}
	close I;
	print "Read gff\n";
	return \%gff;
}

sub readFasta($){
  my ($fil) = @_;
  my %Hseq;
  if (-z $fil){ return \%Hseq;}
  open(FAS,"<","$fil") || die("Couldn't open FASTA file $fil.");
    
     my $temp; 
     my $line; my $hea=<FAS>; chomp ($hea);
      my $trHe = ($hea);
      my @tmp = split(" ",$trHe);
      $trHe = substr($tmp[0],1);
      # get sequence
    while($line = <FAS>)
    {
      #next if($line =~ m/^;/);
      if ($line =~ m/^>/){
        chomp($line);
        $Hseq{$trHe} = $temp;
        $trHe = ($line);
        @tmp = split(" ",$trHe);
	$trHe = substr($tmp[0],1);
	$trHe =~ s/\|//g;
	
#        print $trHe."\n";
        $temp = "";
        next;
      }
    chomp($line);
    $line =~ s/\s//g;
    $temp .= ($line);
    }
    $Hseq{$trHe} = $temp;
  close (FAS);
    return \%Hseq;
}


