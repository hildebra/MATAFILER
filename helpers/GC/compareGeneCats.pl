#!/usr/bin/env perl
#script to blast a GC against other GCs and calculate the overlap between these
use strict;
use warnings;
#use Mods::GenoMetaAss qw(readMap qsubSystem emptyQsubOpt readFastHD prefix_find);
use Mods::IO_Tamoc_progs qw(getProgPaths );
my $lambdaBin = getProgPaths("lambda");
my $lambdaIdxBin = getProgPaths("lambdaIdx"); #lambda ref DB search
my $sortSepScr = getProgPaths("sortSepReLen_scr");#"perl $thisDir/secScripts/sepReadLength.pl";
my $bwt2Bin = getProgPaths("bwt2");#"/g/bork5/hildebra/bin/bowtie2-2.2.9/bowtie2";

my $DB = "/g/bork3/home/hildebra/data/SNP/GCs/SoilCatv2b1/compl.incompl.95.fna";
my @fnaRefs = ("/g/bork3/home/hildebra/DB/freeze11/freeze11.genes.representatives.fa","/g/bork3/home/hildebra/DB/GeneCats/Tara/Tara.fna","/g/bork3/home/hildebra/DB/GeneCats/IGC/1000RefGeneCat.fna");
my @fnRNms = ("frz11","tara","gut1000");
my $odir = "/g/bork3/home/hildebra/data/SNP/GCs/SoilCatv2b1/cmp2otherGC/";
my $tmpDir = "/local/hildebra/cmpGC/";
my @tdirs = ("9305","9350","20250","20205");
my $BlastCores = 40; my $doBwt2 = 0;
my $bwtCore = 40;
system "mkdir -p $odir" unless (-d $odir);
system "mkdir -p $tmpDir" unless (-d $tmpDir);
my $cnt = -1;
foreach my $query (@fnaRefs){
	$cnt++;my $cmd="";
	my $tmptaxblastf = "$odir/$fnRNms[$cnt].m8";
	my $tmptaxsamf = "$odir/$fnRNms[$cnt].sam";
	my $bwt35Log = "$odir/$fnRNms[$cnt].bwt2.log";
	if (!$doBwt2){
		if (!-f $DB.".dna5.fm.sa.val"  ) {
			print "Building LAMBDA index anew (may take up to an hour)..\n";
			my $cmdIdx = "$lambdaIdxBin -p blastn -t $BlastCores -d $DB";
			if (system($cmdIdx)){print("Lamdba ref DB build failed\n$cmdIdx\n",3);}
		}
		$cmd .= "$lambdaBin -t $BlastCores -id 95 -nm 200 -p blastn -e 1e-5 -q $query -d $DB -o $tmptaxblastf\n";
	} else {
		#bowtie2 way
		my $bwtIdx = $DB.".bw2";
		unless (-e $bwtIdx.".rev.2.bt2" && -e $bwtIdx.".4.bt2"&& -e $bwtIdx.".2.bt2"){
			print "Building bwt2 catalog..\n";
			my $cmdIdx = $bwt2Bin."-build --threads $bwtCore -q $DB $bwtIdx\n" ;
			if (system($cmdIdx)){print("Bowtie2 ref DB build failed\n$cmdIdx\n",3);}
			print "Done\n";
		}
		#take care of long reads
		my $tmpFNA = "$tmpDir/bwt24GC.fna";
		system "rm -f $tmpFNA;cp $query $tmpFNA";
		$cmd .= "$sortSepScr 7000 $tmpFNA\n";
		$cmd .= $bwt2Bin." --sensitive --local --norc --no-unal --no-hd --no-sq -p 1 ";
		$cmd .= "-x $bwtIdx -f -U  $tmpFNA.long > $tmptaxsamf 2> $bwt35Log\n";
		#and bulk of reads
		$cmd .= $bwt2Bin." --sensitive --local --norc --no-unal --no-hd --no-sq -p $bwtCore ";
		$cmd .= " -x $bwtIdx -f -U  $tmpFNA >> $tmptaxsamf 2>> $bwt35Log\n";
		
		#fix missing newlines
		$cmd .= "rm -f $tmpFNA\n";
	}
	print "Aligning sample $fnRNms[$cnt]...\n";
	system $cmd;
	print "Done with alignment\n";

}






