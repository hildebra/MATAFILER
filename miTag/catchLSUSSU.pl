#!/usr/bin/env perl
#./catchLSUSSU.pl /g/scb/bork/hildebra/data2/Soil_finland/Project_150108/Sample_S141_170/S141_170_TCCGCGAA-ATAGAGGC_L002_R1_001.fastq.gz /g/scb/bork/hildebra/data2/Soil_finland/Project_150108/Sample_S141_170/S141_170_TCCGCGAA-ATAGAGGC_L002_R2_001.fastq.gz /g/scb/bork/hildebra/data2/Soil_finland/tmp_16s/ /g/scb/bork/hildebra/data2/Soil_finland/tmp_16s/aligned/ 10 1
# ./catchLSUSSU.pl '/scratch/bork/hildebra/SimuB/simulated_metaG3SA_0//tmp/seqClean/filtered.1.fastq' '/scratch/bork/hildebra/SimuB/simulated_metaG3SA_0//tmp/seqClean/filtered.2.fastq' /tmp/hildebra/simulated_metaG3SA_0ITS//SMRNA/ /g/scb/bork/hildebra/SNP///SimuB/simulated_metaG3SA_0/ribos/ 12 Sb1 0 /scratch/bork/hildebra/SimuB//rnaDB/

use strict;
use warnings;
use Mods::GenoMetaAss qw(renameFastqCnts systemW);
use Mods::IO_Tamoc_progs qw(inputFmtSpades getProgPaths);

sub smrnaRunCmd; sub outfileCpy;


#my $tmpP = "/g/scb/bork/hildebra/data2/Soil_finland/tmp_16s/";
my $tmpP = $ARGV[2];
my $alignPath = $ARGV[3];
my $threads = $ARGV[4];
my $smpN = $ARGV[5];
my $doRiboAssembl = $ARGV[6];
my $path2DB = $ARGV[7];

if (@ARGV<7){die "Not enough input arguments!!\n";}

#my $smrPath = "/g/bork5/hildebra/bin/sortmerna-2.0/";
my $smrPath = getProgPaths("srtMRNA_path");

#my $path2DB = "$smrPath/rRNA_databases/";
my $smrnaBin = "$smrPath/./sortmerna";
my $mergeScript = getProgPaths("mergeRdScr");#"$smrPath/scripts/merge-paired-reads.sh";
my $unmergeScript = getProgPaths("unmergeRdScr");#"$smrPath/scripts/unmerge-paired-reads.sh";
my $spadesBin = getProgPaths("spades");#"/g/bork5/hildebra/bin/SPAdes-3.7.0-dev-Linux/bin/spades.py";
my $ltslcaP = "$alignPath/ltsLCA/";
my $singlMode = 0;


if (-e "$alignPath/ITS_pull.sto" && -e "$alignPath/SSU_pull.sto" && -e "$alignPath/LSU_pull.sto" ){
	print "All riboFind sortmerna targets seems to be complete\n";
	if (-e $alignPath."/Ass/allAss.sto"){
		print ", as well as assemblies\n";
		exit(0);
	}
} else  { #unpack reads 
	system "rm -fr $tmpP" if (-d $tmpP);
	system "mkdir -p $tmpP";
	system "mkdir -p $alignPath";
	#preparation of reads
	my @r1i = split(";",$ARGV[0]); my $r1="$tmpP/read1.tmp.fq";
	my @r2i = split(";",$ARGV[1]); my $r2="$tmpP/read2.tmp.fq";
	if ($r2i[0] eq "-1"){$singlMode=1;}
	#die $singlMode."\n";
	#die "$r1i[0]\n";
	my $interLeave = "$tmpP/interleave.fq";
	system "rm -f $r1 $r2 $interLeave";
	if (!$singlMode){
		print "File prep stage 1\n";
		for (my $i=0;$i<@r1i;$i++){
			# merge the files & prepare
			print "File pair $i\n";
			if ($r1i[$i] =~ m/\.gz$/){system "gunzip -c $r1i[$i] | perl -pe 's/\\n/\\t/ if \$. %4' >> $r1"; }#$r1[$i] = "$tmpP/read1.tmp.fq";
			else { system "perl -pe 's/\\n/\\t/ if \$. %4'  $r1i[$i] >> $r1";}
			print "perl -pe 's/\\n/\\t/ if \$. %4'  $r1i[$i] >> $r1\n";
			if ($r2i[$i] =~ m/\.gz$/){system "gunzip -c $r2i[$i] | perl -pe 's/\\n/\\t/ if \$. %4' >> $r2"; } #$r2[$i] = "$tmpP/read2.tmp.fq";
			else { system "perl -pe 's/\\n/\\t/ if \$. %4'  $r2i[$i] >> $r2";}
		}
		print "File prep stage 2\n";
		#print $r1."\n";
		#interleaving
		system "paste -d '\n' $r1 $r2 | tr '\t' '\n' > $interLeave";# bash $mergeScript $r1 $r2 $interLeave";
		system "rm -f $r1" ;#if ($r1i ne $r1);
		system "rm -f $r2" ;#if ($r2i ne $r2);
	} else {
		for (my $i=0;$i<@r1i;$i++){
			if ($r1i[$i] =~ m/\.gz$/){system "gunzip -c $r1i[$i] >>$interLeave";}
			else {system "cat $r1i[$i] >>$interLeave";}
		}
	}
	#die "File prep complete\n";
	#my $refDBits = "/g/bork3/home/hildebra/DB/MarkerG/ITS_fungi/sh_general_release_30.12.2014.fasta,/g/bork3/home/hildebra/DB/MarkerG/ITS_fungi/sh_general_release_30.12.2014.idx";
	
	my $ITSDBfa = getProgPaths("ITSdbFA"); $ITSDBfa =~ m/([^\/]+)$/; $ITSDBfa = $1; my $ITSDBidx = $ITSDBfa; $ITSDBidx =~ s/\.fa*$/\.idx/;

	my $refDBits = "$path2DB/$ITSDBfa,$path2DB/$ITSDBidx";
	
	my $refDBssu = "$path2DB/silva-euk-18s-id95.fasta,$path2DB/silva-euk-18s-id95.idx:$path2DB/silva-bac-16s-id90.fasta,$path2DB/silva-bac-16s-id90.idx:$path2DB/silva-arc-16s-id95.fasta,$path2DB/silva-arc-16s-id95.idx";
	my $refDBlsu = "$path2DB/silva-euk-28s-id98.fasta,$path2DB/silva-euk-28s-id98.idx:$path2DB/silva-bac-23s-id98.fasta,$path2DB/silva-bac-23s-id98.idx:$path2DB/silva-arc-23s-id98.fasta,$path2DB/silva-arc-23s-id98.idx";
	#memory dependent on input file
	#my $outFile = $alignPath."/reds_16S";
	#unless (system $cmd. "--ref $refDBits") {die "Failed\n$cmd\n";}

	my $curStone = "$alignPath/ITS_pull.sto";
	unless (-e $curStone){
		
		if (-e "$path2DB/$ITSDBfa" && -e "$path2DB/$ITSDBidx"){
			my $ltag = "reads_ITS";
			my $runner = smrnaRunCmd($tmpP."/$ltag",$refDBits,$interLeave,$alignPath);print $runner."\n";
			if (system $runner) {print "Error in $runner\n"; exit(90);}#unless (system $runner) {die "Failed\n$runner\n";}
			#renameFastqCnts($alignPath."/reads_ITS.r1.fq",$smpN."__ITS"); renameFastqCnts($alignPath."/reads_ITS.r2.fq",$smpN."__ITS");
			outfileCpy($tmpP."/$ltag",$alignPath);
			system "rm -f $ltslcaP/ITS_ass.sto" if (-e " $ltslcaP/ITS_ass.sto"); #fwd destruction of assignments
			system "touch $curStone" unless (`wc -l $alignPath/$ltag.r1.fq | cut -f1 -d' '` != `wc -l $alignPath/$ltag.r1.fq | cut -f1 -d' '` );
		} else {
			print "Skipping ITS since DB could not be found\n" 
		}
	}

	$curStone = "$alignPath/SSU_pull.sto";
	unless (-e $curStone){
		my $ltag = "reads_SSU";
		my $runner = smrnaRunCmd($tmpP."/$ltag",$refDBssu,$interLeave,$alignPath);#print $runner."\n";
		
		if (system $runner) {print "Error in $runner\n"; exit(88);}#unless (system $runner) {die "Failed\n$runner\n";}
		outfileCpy($tmpP."/$ltag",$alignPath);
		system "rm -f $ltslcaP/SSU_ass.sto" if (-e " $ltslcaP/SSU_ass.sto");
		
		#renameFastqCnts($alignPath."/reads_SSU.r1.fq",$smpN."__SSU"); renameFastqCnts($alignPath."/reads_SSU.r2.fq",$smpN."__SSU");
		system "touch $curStone";
	}
	$curStone = "$alignPath/LSU_pull.sto";
	unless (-e $curStone){
		my $ltag = "reads_LSU";
		my $runner = smrnaRunCmd($tmpP."/$ltag",$refDBlsu,$interLeave,$alignPath);
		print $runner."\n";
		system "rm -f $ltslcaP/LSU_ass.sto" if (-e " $ltslcaP/LSU_ass.sto");
		if (system $runner) {print "Error in $runner\n"; exit(89);}#unless (system $runner) {die "Failed\n$runner\n";}
		outfileCpy($tmpP."/$ltag",$alignPath);
		#renameFastqCnts($alignPath."/reads_LSU.r1.fq",$smpN."__LSU"); renameFastqCnts($alignPath."/reads_LSU.r2.fq",$smpN."__LSU");
		system "touch $curStone";
	}

	system "rm -f $interLeave";

}

system "rm -fr $tmpP";

#ribo assemblies.. difficult as high chance for chimeras.. maybe switch assembler later?
my $outP = $alignPath;
if ($doRiboAssembl && (!-d $outP."/Ass_ITS" || !-e $outP."/Ass_ITS/scaffolds.fasta" || !-e "$outP/Ass/allAss.sto") ){#assembly of ITS, LSU, SSU seperately
	my $K =  "-k 27,33,55,71,125";
	my $logDir = "$outP/Ass/logs/";
	my $Scmd = "mkdir -p $logDir\n\n";
	my $nodeTmp = $outP."Ass/SSU/";
	if (!-z "$outP/reads_SSU.r1.fq"){ #check if input is empty
		$Scmd .= "\nmkdir -p $nodeTmp\n$spadesBin $K ".inputFmtSpades([$outP."/reads_SSU.r1.fq"],[$outP."/reads_SSU.r2.fq"],[],$logDir)." -t $threads --meta -m 60 ";#SSU --mismatch-correction 
		$Scmd .= "-o $nodeTmp\nrm -f -r $nodeTmp/K* $nodeTmp/tmp $nodeTmp/mismatch_corrector $nodeTmp/corrected $nodeTmp/misc \n";
	}
	$nodeTmp = $outP."Ass/LSU/";
	if (!-z "$outP/reads_LSU.r1.fq"){
		$Scmd .= "\nmkdir -p $nodeTmp\n$spadesBin $K ".inputFmtSpades([$outP."/reads_LSU.r1.fq"],[$outP."/reads_LSU.r2.fq"],[],$logDir)." -t $threads --meta -m 60 ";#LSU
		$Scmd .= "-o $nodeTmp\nrm -f -r $nodeTmp/K* $nodeTmp/tmp $nodeTmp/mismatch_corrector $nodeTmp/corrected $nodeTmp/misc \n";
	}
	$nodeTmp = $outP."Ass/ITS/";
	if (!-z "$outP/reads_ITS.r1.fq"){
		$Scmd .= "\nmkdir -p $nodeTmp\n$spadesBin $K ".inputFmtSpades([$outP."/reads_ITS.r1.fq"],[$outP."/reads_ITS.r2.fq"],[],$logDir)." -t $threads --meta -m 60 ";#ITS
		$Scmd .= "-o $nodeTmp\nrm -f -r $nodeTmp/K* $nodeTmp/tmp $nodeTmp/mismatch_corrector $nodeTmp/corrected $nodeTmp/misc \n";
	}
	$Scmd .= "touch $outP/Ass/allAss.sto\n";
	systemW $Scmd;
}

 print "Finished pull out (and eventually assembly)\n";
 exit(0);


sub smrnaRunCmd( $ $ $ $){
	my ($outFile,$refDB,$interLeave,$finD) = @_;
	my $cmd = "set -e\n$smrnaBin --best 1 --reads $interLeave ";
	my $pairOpt = ""; if (!$singlMode){$pairOpt = "--paired_in ";}
	$cmd .= "--blast 1 -a $threads -e 0.1 --log $pairOpt  --fastx --aligned '$outFile' "; #-m 9000 
	$cmd .= "--ref '$refDB'\n";
	$cmd .= "rm -f $outFile.r*\n";
	if (!$singlMode){
		$cmd .= "$unmergeScript $outFile.fq $outFile.r1.fq $outFile.r2.fq\n";#rm -f $outFile.fq";
	}
	#die $cmd;
	return $cmd;
}

sub outfileCpy($ $){
	my ($outF,$finD) = @_;
	
	if ($singlMode){
		#die "$outF , $finD\n";
		if (-e "$outF.fq"){system "mv $outF.fq $finD";} #else {system "touch $outF.fq";}
	} else {
		system "rm -f $outF.fq";
		if (-e "$outF.r1.fq"){system "mv $outF.r1.fq $finD";} #else {system "touch $outF.r1.fq";}
		if (-e "$outF.r2.fq"){system "mv $outF.r2.fq $finD";} #else {system "touch $outF.r2.fq";}
	}
	print "Moved file to $finD\n";
}

