#!/usr/bin/env perl
#assigns LCA tax to LSU / SSU / ITS seq fragments (extracted from metag)
#ex ./lotus_LCA_blast.pl /g/scb/bork/hildebra/Tamoc/FinSoil/Sample_G2677_115/ribos FS31 20
#  ./lotus_LCA_blast.pl /g/scb/bork/hildebra/SNP/SimuG/simulated_FungiSA_5/ribos S5 20 /tmp/hildebra/simu/S5/
use strict;
use warnings;
use threads;
use Mods::GenoMetaAss qw(reverse_complement_IUPAC readFasta systemW);
use Mods::IO_Tamoc_progs qw(getProgPaths);



sub getTaxForOTUfromRefBlast;
sub splitBlastTax;
sub writeBlastHiera; sub runBlastLCA;
sub merge;
sub fastq2fna; sub flashed;
sub rework_tmpLines;

#binaries
my $flashBin = getProgPaths("flash");#"/g/bork3/home/hildebra/bin/FLASH-1.2.10/flash";
my $lambdaBin = getProgPaths("lambda");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda";
my $lambdaIdxBin = $lambdaBin."_indexer";#getProgPaths("");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda_indexer";
my $LCAbin = getProgPaths("LCA");#"/g/bork3/home/hildebra/dev/C++/LCA/./LCA";
my $srtMRNA_path = getProgPaths("srtMRNA_path");


#some more pars for LCA	
my @idThr = (97,95,93,91,88,78);
my $lengthTolerance = 0.85;
my $maxHitOnly = 0;
my $DoPar = 0;

if (@ARGV ==0 ){die"no Input arguments given\n";}

my $inD = $ARGV[0]; #expects TAMOC output files in this dir
$inD .= "/" unless ($inD =~ m/\/$/);
#dir to store LCA
my $outdir = $inD."ltsLCA/";
#dir to store merged reads in
my $mergeD = $inD."flashMerge/";
my $SmplID = $ARGV[1];#Name of the Sample
my $BlastCores = $ARGV[2]; #parrallel execution
my $tmpD = ""; #tmp dir for faster I/O
if (@ARGV > 4){
	$tmpD = $ARGV[4] ;
	system "rm -rf $tmpD;mkdir -p $tmpD";
}
my $singlReads=0;
if (@ARGV > 5){
	$singlReads = !$ARGV[5];
}
#die "$singlReads\n";
my $DBdir = $ARGV[3];

#datbases
#datbases
my $LSUdbFA = "$DBdir/SLV_128_LSU.fasta";
my $LSUtax = "$DBdir/SLV_128_LSU.tax";
my $SSUdbFA= "$DBdir/SLV_128_SSU.fasta";
my $SSUtax = "$DBdir/SLV_128_SSU.tax";
#my $ITSdbFA = "$DBdir/sh_refs_qiime_ver7_99_02.03.2015.fasta";my $ITStax = "$DBdir/sh_taxonomy_qiime_ver7_99_02.03.2015.txt";
my $ITSdbFA = "$DBdir/ITS_comb.fa";my $ITStax = "$DBdir/ITS_comb.tax";
my $PR2dbFA = "$DBdir/gb203_pr2_all_10_28_99p.fasta";
my $PR2tax = "$DBdir/PR2_taxonomy.txt";

#some cleanups (includes prev runs)
system "rm -f  $inD/*.blast $outdir/*riboRun_bl"; #$inD/reads_SSU.fq $inD/reads_LSU.fq $inD/reads_ITS.fq

system "mkdir -p $outdir" unless (-e $outdir);

my $blMode = 2; #2=lambda,1=Blast,3=sortmeRNA
my $inputOK=1; #flag
my $curStone = ""; #flag file for each subpart (SSU,LSU,ITS)
my @allInFa;
if ($singlReads){
	@allInFa = ($inD."reads_LSU.fq",$inD."reads_ITS.fq",$inD."reads_SSU.fq");
} else {
	@allInFa = ($inD."reads_LSU.r1.fq",$inD."reads_LSU.r2.fq",$inD."reads_ITS.r1.fq",$inD."reads_ITS.r2.fq",$inD."reads_SSU.r1.fq",$inD."reads_SSU.r2.fq");
}

if (!$singlReads){
	#pretty circular safety catch.. maybe remove later?
	if (!-z $allInFa[2] && (-z $allInFa[0] ||  -z $allInFa[1])){#LSU obviously wrong
		system "rm -f $inD/LSU_pull.sto $outdir/LSU_ass.sto";
	}
	if (!-z $allInFa[2] && (-z $allInFa[4] || -z $allInFa[5])){#SSU obviously wrong
		system "rm -f $inD/SSU_pull.sto $outdir/SSU_ass.sto";
	}

	#all are empty (also wrong)
	if (-z $inD."reads_LSU.r1.fq" && -z $inD."reads_LSU.r2.fq" && -z $inD."reads_ITS.r1.fq" && -z $inD."reads_ITS.r2.fq"
				&& -z $inD."reads_SSU.r1.fq" && -z $inD."reads_SSU.r2.fq" &&
				-z $inD."reads_LSU.r1.fq.gz" && -z $inD."reads_LSU.r2.fq.gz" && -z $inD."reads_ITS.r1.fq.gz" && -z $inD."reads_ITS.r2.fq.gz"
				&& -z $inD."reads_SSU.r1.fq.gz" && -z $inD."reads_SSU.r2.fq.gz"){
		print "input seems completely wrong.. run has to be repeated\n";
		system "rm -f $mergeD/*_pull.sto";
		exit(11);
	}
}

#all done already
if (-e "$outdir/Assigned.sto" && -e "$outdir/LSU_ass.sto" && -e "$outdir/ITS_ass.sto" && -e "$outdir/SSU_ass.sto"){
	print "All assigned already\n";
	exit (0);
}



#------------- LSU - SSU - ITS ------------------
my @tags=("ITS","SSU","LSU");
for (my $i=0;$i<3;$i++){
	$curStone = "$outdir/$tags[$i]_ass.sto";
	unless (-e $curStone && -e "$outdir/$tags[$i]riboRun_bl.hiera.txt"){
		my @dbfa; my @dbtax;
		
		if ($tags[$i] eq "LSU") {@dbfa = ($LSUdbFA); @dbtax = ($LSUtax) ;
		}elsif ($tags[$i] eq "SSU"){ @dbfa = ($PR2dbFA,$SSUdbFA); @dbtax = ($PR2tax,$SSUtax) ;
		}elsif ($tags[$i] eq "ITS"){ @dbfa = ($ITSdbFA); @dbtax = ($ITStax) ;
		}else {die "Unrecognized tag $tags[$i]\n";}
		
		#merge
		my $go=1;
		if (!$singlReads){
			$go = flashed($inD."reads_$tags[$i].r1.fq",$inD."reads_$tags[$i].r2.fq",$mergeD,$tags[$i],$inD); 
			if ($go==1){
				#sim search & LCA
				runBlastLCA($mergeD.$tags[$i],\@dbfa,\@dbtax,"$tags[$i]riboRun_bl",$blMode,$SmplID,$tags[$i],$go);
				
			} elsif($go==2) {$inputOK=0;print"Problem with $tags[$i] primary input; run needs to be repeated";}
		} else {#single reads, also have specific input fmt
			my $inFilX = $inD."reads_$tags[$i].fq";
			$inFilX = $inD."reads_$tags[$i].fq.gz" if (-e $inD."reads_$tags[$i].fq.gz" && !-e $inFilX);
			runBlastLCA($inFilX,\@dbfa,\@dbtax,"$tags[$i]riboRun_bl",$blMode,$SmplID,$tags[$i],$go);
			#$inD."reads_$tags[$i].fq"
		}
		system "touch $curStone";
	}
#	die;
}

#die();


#die();
#cleanup
system "rm -r $tmpD" if ($tmpD ne "");
if ($inputOK){
	print "$outdir/Assigned.sto";
	system "touch $outdir/Assigned.sto";
	system "gzip -q ".join (" ",@allInFa);
	exit(0);
}
print "Finished\n";
exit(0);


#flash merge & some file checks
sub flashed($ $ $ $ $){
	my ($r1,$r2,$outD,$outT,$primD) = @_;
	if ( -e $r2.".gz"){system "rm -f $r2; gunzip $r2.gz";}
	if ( -e $r1.".gz"){system "rm -f $r1; gunzip $r1.gz";}
	if (    !-e $r2 || !-e $r1 || (-z $r1 && !-z $r2) || (!-z $r1 && -z $r2) 
	||  (`wc -l $r1 | cut -f1 -d' '` != `wc -l $r2 | cut -f1 -d' '` ) ){
		print "Empty input $r1\n";
		system "rm -f $primD/$outT"."_pull.sto";
		#exit 33;
		return 2;
	}
	if (-z $r1 && -z $r2){return 0;}
	print "running flash..\n";
	my $mergCmd = "$flashBin -M 200 -o $outT -d $outD -t $BlastCores $r1 $r2";
	if (system $mergCmd){
		print "\n$mergCmd\nfailed\n";
		system "rm -f $primD/$outT"."_pull.sto";
		return 2;
	}# or 
	return 1;
}

sub fastq2fna($){
	my ($in) = @_;
	#deactivate, since lambda can just read fq...
	return $in;
	print $in."\n";
	my $out = $in;
	$out =~ s/f[ast]?q$/fna/;
	open I,"<$in" or die "Input fastq file not available";
	my $l = <I>; close I; if ($l =~ m/^>/){return $out;}
	system "cat $in | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\\t' '\\n' > $out";
	#die $out;
	return $out;
}

#find reads that could not match any LCA and redo with different DB
sub findUnassigned($ $ $ ){
	my ($BTr,$Fr,$outF) = @_;
	my %BT = %{$BTr}; my %Fas = %{$Fr};
	#my @t =keys %Fas; print "$t[0]\n";
	my $cnt =0; my $dcn =0 ;
	my @kk = keys %Fas;
	my $totFas = @kk;
	
	if ($outF eq ""){
		foreach my $k (keys %BT){
			my @curT = @{$BT{$k}};
			if ( @curT==0 || $curT[0] eq "?" ){$dcn++;} #|| $curT[1] eq "?"
			$cnt++;
		}
		if (@kk == 0){
			print "Total of ". ($cnt-$dcn)." / $cnt reads have LCA assignments\n";
		} else {
			print "$dcn / $cnt reads failed LCA assignments, checked $totFas reads.\n";
		}
		return;
	}
	foreach my $k (keys %BT){
		my @curT = @{$BT{$k}};
		#print $k."\t${$BT{$k}}[0]   ${$BT{$k}}[2]\n" 
		if ( @curT==0 || $curT[0] eq "?" ){#|| $curT[1] eq "?"){
			delete $BT{$k};
			$dcn++;
			
			#print ">".$k."\n".$Fas{$k}."\n";
		} #else {print $k."\t${$BT{$k}}[0]   ${$BT{$k}}[2]\n" ;}
		else {
			die "Can't find fasta entry for $k\n" unless (exists $Fas{$k});
			delete $Fas{$k};
		}
		$cnt ++;
		#die if ($cnt ==100);
	}
	@kk = keys %Fas;
	print "$dcn / $cnt reads failed LCA assignments\nWriting ".@kk." of previous $totFas reads for next iteration.\n";
	open O,">$outF" or die "can;t open unassigned fasta file $outF\n";
	foreach my $k(@kk){
		print O ">".$k."\n".$Fas{$k}."\n";
	}
	close O;
	return ($outF,\%BT,$dcn,\%Fas);
}

#main routine that does sim search & starts the LCA
sub runBlastLCA(){

	my ($queryO,$DBar,$DBtaxar,$id,$doblast,$SmplID,$MKname,$go) =@_;
	my $taxblastf_base = $outdir."$id";
	if ($tmpD ne ""){
		$taxblastf_base = $tmpD."$id";
	}
	my $doInter=1; my $doQuery=1; #fine control of what subparts to do..
	my $r1 = $queryO; my $r2 = $queryO;
	my $interLeaveO = "";
	my $hof1 = "$outdir/$id.hiera.txt"; my $hof2 = "$outdir/$id"."_inter.hiera.txt";
	#check for empty input
	if ($go == 2){ return ;}
	if (!$go){
		print "$id has empty input files\n";
		system "touch $hof1";
		return;
	}

	if ($singlReads){
		die "can't find single read input file: $queryO\n" unless (-e $queryO);
		$queryO = fastq2fna($queryO);
		#die "$queryO\n";
	}elsif ($queryO !~ m/\.fna$/){#merged fastq that need to be processed
		$queryO .= ".extendedFrags.fastq";
		$r1 .= ".notCombined_1.fastq";
		$r2 .= ".notCombined_2.fastq";
		$interLeaveO = $outdir."inter$id.fna";
		#print $r1." V\n";
		$r1 = fastq2fna($r1); $r2 = fastq2fna($r2); $queryO = fastq2fna($queryO);

		merge($r1,$r2,$interLeaveO) if ($doInter);
	} else {
		print "no interLeaveO\n";
		$doInter=0;
	}
	if (-z $interLeaveO){print "No interleaved files $id\n";$doInter=0;}
	if (-z $queryO){print "No merged files $id\n";$doQuery=0;}
	#my $BlastCores = 20;
#die "$queryO\n";
	#my $fasrA = readFasta( $queryO,1); 
	#my %fas=%{$fasrA}; my @t = keys %fas; die "$queryO\n@t\n$t[0]\n";
	#my $fasrI = readFasta($interLeaveO,1); 
	#print $DBar."\n";
	my @DBa = @{$DBar}; my @DBtaxa = @{$DBtaxar};
	my $query= $queryO; my $interLeave = $interLeaveO;
	my $BlastTaxRi = {}; my $BlastTaxR = {};
	my $fullBlastTaxRi = {}; my $fullBlastTaxR = {};
	my $simName = ""; my $leftover = 111;
	my $contInter = 1; my $contQuery = 1;
	my @taxouts=(); my @taxouts_inter = ();
	for (my $DBi=0;$DBi<@DBa; $DBi++){
		my $DB = $DBa[$DBi]; my $DBtax = $DBtaxa[$DBi];
		print "Running sim search $doblast..\n";
		my $taxblastf = $taxblastf_base.".$DBi.m8.gz";
		my $taxblastf2=$taxblastf_base.".$DBi.i.m8.gz";
		push @taxouts, $taxblastf;
		push @taxouts_inter, $taxblastf2;
		
		if ($interLeave ne "" && -e $interLeave){$doInter=1;} else {$doInter=0;}
		if (!-e $query){$doQuery=0;} else {$doQuery=1;}
		#die "$doQuery $doInter\n$interLeave\n";
		if ($doblast==1){
			die "using blast is deprecated and no longer supported!\n";
			$taxblastf.=".blast"; $simName="blast";
			print "Running Blast\n";
			my $mkBldbBin = "/g/bork5/hildebra/dev/lotus//bin//ncbi-blast-2.2.29+/bin/makeblastdb";
			my $blastBin = "/g/bork5/hildebra/dev/lotus//bin//ncbi-blast-2.2.29+/bin/blastn";
			my $cmd = "$mkBldbBin -in $DB -dbtype 'nucl'\n";
			unless (-f $DB.".nhr"){	system($cmd);}
			my $strand = "both";
			#-perc_identity 75
			$cmd = "";
			if ($doQuery && $contQuery){
				$cmd .= "$blastBin -query $query -db $DB -out $taxblastf -outfmt 6 -max_target_seqs 50 -evalue 0.1 -num_threads $BlastCores -strand $strand \n"; #-strand plus both minus
			}
			unless ( -e $taxblastf){	print "Running blast on combined files\n";
				system($cmd);
			} else {	print "Blast output $taxblastf does exist\n";}
			if ($doInter){
				$cmd .= "$blastBin -query $interLeave -db $DB -out $taxblastf2 -outfmt 6 -max_target_seqs 50 -evalue 0.1 -num_threads $BlastCores -strand $strand \n"; #-strand plus both minus
				unless ( -e $taxblastf2){	print "Running blast on interleaved files\n";system($cmd);
				} else {	print "Blast output $taxblastf2 does exist\n";}
			}
		} elsif ($doblast==2){
			$simName = "lambda";
			if (!-f $DB.".dna5.fm.sa.val"  ) {
				print "Building LAMBDA index anew (may take up to an hour)..\n";
				my $cmdIdx = "$lambdaIdxBin -p blastn -t $BlastCores -d $DB";
				#system "touch $DB.dna5.fm.lf.drv.wtc.24";
				if (system($cmdIdx)){die ("Lamdba ref DB build failed\n$cmdIdx\n");}
			}
			print "Starting LAMBDA similarity search..\n";
			my $tmptaxblastf = "$outdir/tax.m8";
			$tmptaxblastf = "$tmpD/tax.m8" unless ($tmpD eq "");
			my $cmd = "";
			$cmd = "";
			if ($doQuery && $contQuery){
				$cmd .= "$lambdaBin -t $BlastCores -id 75 -nm 200 -p blastn -e 1e-5 -so 7 -sl 16 -sd 1 -b 5 -pd on -q $query -d $DB -o $taxblastf\n";
			}
			#$cmd .= "\nmv $tmptaxblastf $taxblastf\n";
			if ($doInter && $contInter){
				$cmd .= "$lambdaBin -t $BlastCores -id 75 -nm 200 -p blastn -e 1e-5 -so 7 -sl 16 -sd 1 -b 5 -pd on -q $interLeave -d $DB -o $taxblastf2\n";
				#$cmd .= "\nmv $tmptaxblastf $taxblastf2\n";
				#unless ( -e $taxblastf2){	
				print "Running blast on interleaved files\n";
			}
			#die "$doQuery && $contQuery\n$cmd\n";
			system($cmd);
			#} else {	print "Blast output $taxblastf2 does exist\n";}}
			#system $cmd ;#or die "\n$cmd\n failed\n";
		} elsif ($doblast==3) {
			$simName = "smRNA";
			my $smrnaBin = "$srtMRNA_path/sortmerna";
			my $smrnaMkDB = "$srtMRNA_path/indexdb_rna";
			my $refDB = "$DB,$DB.sidx";
			unless (-e "$DB.sidx.stats"){	print "making DB..";system "$smrnaMkDB --ref $refDB";	print " done\n";}
			#die $outdir."inter$id.fna\n";
			print "Running sortmerna merge\n";
			my $cmd = "$smrnaBin --best 50 --reads $interLeave ";
			$cmd .= "--blast 1 -a 20 -e 0.1 -m 10000 --paired_in --fastx --aligned $taxblastf.i ";
			$cmd .= "--ref $refDB\n";
			print "$taxblastf.i\n";
			system $cmd;
			print "Running sortmerna\n";
			$cmd = "";
			if ($doQuery && $contQuery){
				$cmd .= "$smrnaBin --best 50 --reads $query ";
				$cmd .= "--blast 1 -a 20 -e 0.1 -m 10000 --aligned $taxblastf ";
				$cmd .= "--ref $refDB\n";
			}
			system $cmd;
			$taxblastf.=".blast";
			print "done sortmeRNA assignment\n";
		}

	}

	print "Running LCA..\n";
	my $LCcmd = "";
	
	

	#$hof1 = writeBlastHiera($fullBlastTaxR,$id,$simName);  undef $BlastTaxR ;
	$LCcmd .= "$LCAbin -i ". join(",",@taxouts) . " -r ".join(",",@DBtaxa)." -o $hof1 -LCAfrac 0.8 -showHitRead\n";
	#merging of interleave results
	if ($doInter){
		$LCcmd .= "$LCAbin -i ". join(",",@taxouts_inter) . " -r ".join(",",@DBtaxa)." -o $hof2 -LCAfrac 0.8 -showHitRead\n";
		$LCcmd.= "cat $hof2 >> $hof1; rm $hof2;";
	}
	
	#die $LCcmd."\n";
	systemW $LCcmd;
	system "rm -f @taxouts @taxouts_inter" ; #die();
	system "rm -f $queryO"."__U* $interLeaveO"."__U*";
}

#file operations on paired end files
sub merge($ $ $){
	my ($r1,$r2,$interleave) = @_;
	print "Interleaving 2 files..\n";
	open I1,"<$r1" or die "1   $r1\n".$!;open I2,"<$r2"or die "2   $r2\n".$!;open O,">$interleave"or die "interl   ".$!;
	my $line1=<I1>;my $line2=<I2>; $line2="";#read header & get rid of it..
	while(1){
		my $tmp="";
		while ($tmp !~ m/^>/){chomp $tmp;$line1.=$tmp;$tmp=<I1>; last unless defined $tmp;}
		#print $line1."\n";
		print O $line1; if (defined $tmp ){chomp $tmp; $line1 = $tmp."\n"; }
		$tmp=""; 
		while ($tmp !~ m/^>/){chomp $tmp; $line2.=$tmp;$tmp=<I2>;last unless defined $tmp;}
		#get one more line to get rid of header:
		print O reverse_complement_IUPAC($line2)."\n"; $line2 = "";#$tmp;
		last unless defined $tmp;
	}
	close O; close I1; close I2;
#	my $mergeScript = "/g/bork5/hildebra/bin/sortmerna-2.0/scripts/merge-paired-reads.sh";
#	system "bash $mergeScript $r1 $r2 $interLeave";

}

#function that goes line by line through Blast m8 table
sub getTaxForOTUfromRefBlast($ $ $ $){
	my ($blastout,$GGref,$interLMode,$BlastCores) = @_;
	#sp,ge,fa,or,cl,ph
	my %GG = %{$GGref};
	my @ggk = keys(%GG);
	my $maxGGdep = scalar(@{$GG{$ggk[0]}});
	@ggk = ();
	open B,"<",$blastout or die "Could not read $blastout\n";
	my $sotu = "";my $sID=0 ; my $sLength=0; 
	#my @sTax=(); my @sMaxTaxNum = ();
	my %retRef ;# my $retDRef ={};
	my $cnt=0;
	my @tmpLines = (); #stores Blast lines
	my @spl; #temp line delim
	my $minBit = 120; my $minEval = 1e-14;
	my %prevQueries = ();
	my $line = "";
	my @thrs; my $thrsCnt = 0;
	#for (my $i=0;$i<$BlastCores;$i++){$thrs[$i] = threads->create(\&empty_proc);}
	while ($line = <B>){
		$cnt++; chomp $line; #$line2 = $line;
		my @spl = split("\t",$line);
		my $totu = $spl[0]; #line otu
		$totu =~ s/^>//;
		if ($cnt == 1) {$sotu = $totu;}
		#print $line." XX $spl[11] $spl[10]\n"; die () if ($cnt > 50);
		
		#check if this is a 2 read-hit type of match (interleaved mode) & merge subsequently
		if ($interLMode){
			if (@tmpLines>0 && exists($prevQueries{ $spl[1]}) ){
				my @prevHit = @{$prevQueries{$spl[1]}};
				die "something went wrong with the inter matching: $prevHit[1] - $spl[1]\n" unless ( $prevHit[1] eq $spl[1] );
				$prevHit[11] += $spl[11];#bit score
				$prevHit[3] += $spl[3];$prevHit[5] += $spl[5];$prevHit[4] += $spl[4];#alignment length,mistmatches,gap openings
				$prevHit[2] = ($prevHit[2] + $spl[2]) / 2;
				#$tmpLines[-1] = \@prevHit;
				@spl = @prevHit;
				#next;
			} else {$prevQueries{ $spl[1] } = \@spl;}
		}
		if (($spl[11] < $minBit) || ($spl[10] > $minEval) ){ #just filter out..
			#print "ss\n"; 
			#next;
			$spl[2] =0; #simply deactivate this way...
		}
		
		if ($sotu eq $totu){
			push(@tmpLines,\@spl);
			if ($spl[2] > $sID && $spl[3] > $sLength){$sID = $spl[2]; $sLength = $spl[3];}
			if ($spl[3] > ($sLength*1.4) && $spl[2] > ($sID*0.9)) {$sID = $spl[2]; $sLength = $spl[3];} #longer alignment is worth it..
			#print $sID."\n";
		} else {
			#print "DD";
			#print "Maybe\n";
			
#			if (0 && $sotu ne ""){				if ($thrsCnt>= $BlastCores){$thrsCnt=0;}				#get result of old job
#				my $ret = $thrs[$thrsCnt]->join();				foreach (keys %{$ret}){$retRef{$_} = ${$ret}{$_};}				$thrs[$thrsCnt] = threads->create(\&rework_tmpLines,\@tmpLines,$sotu,$sID,$sLength,\%GG,$maxGGdep,$sMaxTax) ;
#				$thrsCnt++; print $thrsCnt." ";			}
			if ($sotu ne ""){
				my ($AR) = rework_tmpLines(\@tmpLines,$sID,$sLength,\%GG,$maxGGdep,); 
				$retRef{$sotu} = $AR;
			}
			$sotu = $totu; undef @tmpLines ; undef %prevQueries ;
			push(@tmpLines,\@spl);$prevQueries{ $spl[1] } = \@spl;
			$sID =  $spl[2]; $sLength = $spl[3];
		}
	}
	#last OTU in extra request
	if ($sotu ne ""){
		#my @spl = split("\t",$line);
		
		my ($AR) = rework_tmpLines(\@tmpLines,$sID,$sLength,\%GG,$maxGGdep);
		$retRef{$sotu} = $AR;
#		my($sTaxX,$sMaxTaxX) = LCA(\@sTax,\@sMaxTaxNum,$maxGGdep);
#		$ret{$sotu} = $sTaxX;
#		$retD{$sotu} = $sMaxTaxX;
	}
	close B;
	#debug 
	#my %ret = %{$retRef};	my @tmp = @{$ret{$sotu}};print "\n@tmp  $sotu\n";
	undef  @tmpLines;undef %prevQueries ;
	#print "Assigned $refDBname Taxonomy to OTU's\n",0;
	
	return (\%retRef);
}

