#!/usr/bin/env perl
#perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/buildTree.pl /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//allFNAs.fna /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//allFAAs.faa /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//categories4ete.txt /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2/tesssst/ 12 0 0.8
#ARGS: ./buildTree.pl [FNA] [FAA] [categoryFile] [outDir] [CPUs] [useEte? [1=ETE,0=this script]] [filter, ignore]

use warnings;
use strict;
use threads ('yield',
                 'stack_size' => 64*4096,
                 'exit' => 'threads_only',
                 'stringify');
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw( systemW readFasta);


sub convertMultAli2NT;
sub mergeMSAs;
sub synPosOnly;

my $doPhym= 0;

my $pal2nal = getProgPaths("pal2nal");#"perl /g/bork3/home/hildebra/bin/pal2nal.v14/pal2nal.pl";
#die $pal2nal;
my $clustaloBin = getProgPaths("clustalo");#= "/g/bork3/home/hildebra/bin/clustalo/clustalo-1.2.0-Ubuntu-x86_64";
my $fastq2phylip = getProgPaths("fastq2phylip_scr");
my $phymlBin = getProgPaths("phyml");
my $raxmlBin = getProgPaths("raxml");
my $msapBin = getProgPaths("msaprobs");
my $trimalBin = getProgPaths("trimal");
my $pigzBin  = getProgPaths("pigz");

#die "TODO $trimalBin\n";
#trimal -in /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/MSA/COG0185.faa -out /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/MSA/tst.fna -backtrans /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/inMSA0.fna -keepheader -keepseqs -noallgaps -automated1 -ignorestopcodon
#some runtim options...
#my $ncore = 20;#RAXML cores
 my $ntFrac =2; 
 my $clustalUse = 1; #do MSA with clustal (1) or msaprobs (0) 
 if ($clustalUse == 0){print "Warning:  MSAprobs with trimal gives warnings (ignore them)\n";}

my $ntCnt =0;
if (@ARGV < 2){die "Not enough input args!!!\n";}
my ($fnFna, $aaFna,$cogCats,$outD,$ncore,$Ete, $filt) = @ARGV;
if ($filt <1){$ntFrac=$filt; print "Using filter with $ntFrac fraction of nts\n";}
else {$ntCnt = $filt;}
my $tmpD = $outD;
system "mkdir -p $outD" unless (-d $outD);
my $cmd =""; my %usedGeneNms;

#------------------------------------------
#sorting by COG, MSA & syn position extraction

if (!$Ete){
	my $treeD = "$outD/TMCtree/";
	my $raxD = "$treeD";
	my $raxLogD = "$raxD/Logs/";
	system "rm -fr $treeD" if (-d $treeD);
	system "mkdir -p  $raxLogD" unless (-d $raxLogD);
	system "mkdir -p  $outD/MSA/" unless(-d "$outD/MSA/");
	my $multAli = "$outD/MSA/MSAli.fna";
	my $multAliSyn = $multAli.".syn.fna";

	my @MSAs; my @MSAsSyn;#full MSAs and MSAs with syn pos only
	my @MSrm;
	#my @xx = keys %FAA; die "$xx[0] $xx[1]\n$FAA{HM29_COG0185}\n";
	if (1 && $cogCats ne ""){
		my $hr = readFasta($aaFna); my %FAA = %{$hr};
		$hr = readFasta($fnFna); my %FNA = %{$hr};
		open I,"<$cogCats" or die "Can't open cogcats $cogCats\n";
		my $cnt = 0; 
		my %samples; my %genCats;
		while (<I>){
			chomp; my @spl = split /\t/;
			next if (@spl ==0);
			my @spl2 = split /_/,$spl[0] ;
			die "Double gene name in tree build pre-concat: $spl2[1]\n" if (exists($usedGeneNms{$spl2[1]}));
			$usedGeneNms{$spl2[1]} = 1;
			my $tmpInMSA = "$tmpD/inMSA$cnt.faa";
			my $tmpInMSAnt = "$tmpD/inMSA$cnt.fna";
			my $tmpOutMSA2 = "$outD/MSA/$spl2[1].$cnt.faa";
			my $tmpOutMSA = "$outD/MSA/$spl2[1].$cnt.fna";
			my $tmpOutMSAsyn = "$tmpD/MSA/$spl2[1].$cnt.syn.fna";
			open O,">$tmpInMSA" or die "Can;t open tmp faa file for MSA: $tmpInMSA\n";
			open O2,">$tmpInMSAnt" or die "Can;t open tmp fna file for MSA: $tmpInMSAnt\n";
			foreach my $seq (@spl){
				my @spl2 = split /_/,$seq; $genCats{$spl2[1]} = 1; $samples{$spl2[0]} = 1;
				die "can't find AA seq $seq\n" unless (exists ($FAA{$seq}));
				die "can't find fna seq $seq\n" unless (exists ($FNA{$seq}));
				$FAA{$seq} =~ s/\*//g if (!$clustalUse);
				print O ">$seq\n$FAA{$seq}\n";
				print O2 ">$seq\n$FNA{$seq}\n";
			}
			close O;close O2;
			if ($clustalUse){
				$cmd = $clustaloBin." -i $tmpInMSA -o $tmpOutMSA2 --outfmt=fasta --threads=$ncore --force\n";
			} else {
				$cmd = "$msapBin -num_threads $ncore $tmpInMSA > $tmpOutMSA2";
			}
			#die $cmd;
			system $cmd unless (-s $tmpOutMSA2);
			#die;
			convertMultAli2NT($tmpOutMSA2,$tmpInMSAnt,$tmpOutMSA);
			
			#die("XX\n") if ($spl2[1] eq "COG0081");
			synPosOnly($tmpOutMSA,$tmpOutMSA2,$tmpOutMSAsyn,0);
			system "rm -f $tmpInMSA $tmpInMSAnt";# $tmpOutMSA2";
			push (@MSAs,$tmpOutMSA);
			push (@MSAsSyn,$tmpOutMSAsyn);
			push (@MSrm,$tmpOutMSA2,$tmpOutMSA);
			$cnt ++;
			print "$cnt ";
		}
		close I;
		#die;
		#merge cogcats - can go to tree from here
		mergeMSAs(\@MSAs,\%samples,$multAli,0);
		mergeMSAs(\@MSAsSyn,\%samples,$multAliSyn,1);
		#die();
		system "rm -f ".join(" ",@MSrm)." ".join(" ",@MSAsSyn); 
	} elsif (1) {#no marker way, single gene
		my $tmpInMSA = $aaFna;
		my $tmpInMSAnt = $fnFna;
		my $tmpOutMSA2 = "$tmpD/outMSA.faa";
		my $tmpOutMSA = $multAli;#"$tmpD/outMSA.fna";
		my $tmpOutMSAsyn = $multAliSyn;#"$tmpD/outMSA.syn.fna";
		my $numFas = `grep -c '^>' $tmpInMSA`;
		chomp $numFas;
		if ($numFas <= 1){print "Not enough Sequences\n"; exit(0);}
		if ($clustalUse){
			$cmd = $clustaloBin." -i $tmpInMSA -o $tmpOutMSA2 --outfmt=fasta --threads=$ncore --force\n";
		} else {
			$cmd = "sed -i 's/\\*//g' $tmpInMSA\n";
			$cmd .= "$msapBin -num_threads $ncore $tmpInMSA > $tmpOutMSA2\n";
		}
		#die $cmd;
		system $cmd; print "finished MSA\n";
		convertMultAli2NT($tmpOutMSA2,$tmpInMSAnt,$tmpOutMSA);
		synPosOnly($tmpOutMSA,$tmpOutMSA2,$tmpOutMSAsyn,0);
		#system "rm $tmpInMSA $tmpInMSAnt $tmpOutMSA2";
		system "rm $tmpOutMSA2";
		push (@MSAs,$tmpOutMSA);
		#$multAli = $tmpOutMSA; $multAliSyn = $tmpOutMSAsyn;
	}

	
	#-------------------------------------------
	#Tree building part with RaxML
	#die $multAli."\n";
	#convert fasta again
	my $nwkFile = "$treeD/tree_phyml_all.nwk";
	my $nwkFile2 = "$treeD/tree_phyml_syn.nwk";
	
	my $raTmpF = "RXMtmp"; my $raxFile = "RXMall";
	my $raxFile2 = "RXMsyn";
	my @thrs;
	my $tcmd = "$fastq2phylip -c 50 $multAli > $multAli.ph\n";
	$tcmd .= "$fastq2phylip -c 50 $multAliSyn >$multAliSyn.ph\n";
	if (system $tcmd) {die "fasta2phylim failed:\n$tcmd\n";}

	my $raxDef = " --silent -m GTRGAMMA";
	#raxml - on all sites
	$tcmd="";
	if (1){
		$tcmd =  "$raxmlBin -T$ncore -f d -p 31416 -s $multAli.ph $raxDef -n $raTmpF -w $raxD > $raxLogD/ini.log\n";
		if (system $tcmd){#failed.. prob optimization problem, use other tree instead
			$raxDef = " --silent -m GTRGAMMAI";
			system "rm $raxD/*$raTmpF";
			systemW "$raxmlBin -T$ncore -f J -p 31416 -s $multAli.ph $raxDef -n $raxFile -w $raxD -t $raxD/RAxML_bestTree.$raTmpF > $raxLogD/optimized.log\n";
			#system "mv $raxD/RAxML_bestTree.$raTmpF $raxD/RAxML_fastTreeSH_Support.$raxFile";
		}
		systemW "$raxmlBin -T$ncore -f J -p 31416 -s $multAli.ph $raxDef -n $raxFile -w $raxD -t $raxD/RAxML_bestTree.$raTmpF > $raxLogD/optimized.log\n";
		systemW  "$raxmlBin -T$ncore -f x -p 31416 -s $multAli.ph $raxDef -n all -w $raxD -t $raxD/RAxML_bestTree.$raTmpF > $raxLogD/optimized.log\n";
		#push(@thrs, threads->create(sub{system $tcmd;}));
	}
	#reset to GTRGAMMA model
	$raxDef = " --silent -m GTRGAMMA";
	
	
	#raxml - on syn sites only
	$tcmd = "$raxmlBin -T$ncore -f d -p 31416 -s $multAliSyn.ph $raxDef -n $raTmpF.s -w $raxD  > $raxLogD/ini_syn.log\n";
	if (system $tcmd){#failed.. prob optimization problem, use other tree instead
		$raxDef = " --silent -m GTRGAMMAI";
		system "rm $raxD/*$raTmpF.s";
		$tcmd= "$raxmlBin -T$ncore -f d -p 31416 -s $multAliSyn.ph $raxDef -n $raTmpF.s -w $raxD  > $raxLogD/ini_syn.log\n";
		print $tcmd."\n\n\n";
		systemW $tcmd;
		#system "mv $raxD/RAxML_bestTree.$raTmpF.s $raxD/RAxML_fastTreeSH_Support.$raxFile2";
		#print "Fallback to original optimized tree\n";
	}
	systemW "$raxmlBin -T$ncore -f J -p 31416 -s $multAliSyn.ph $raxDef -n $raxFile2 -w $raxD -t $raxD/RAxML_bestTree.$raTmpF.s > $raxLogD/optimized_syn.log\n";
	systemW "$raxmlBin -T$ncore -f x -p 31416 -s $multAliSyn.ph $raxDef -n syn -w $raxD -t $raxD/RAxML_bestTree.$raTmpF.s > $raxLogD/optimized_syn.log\n";
#die $tcmd."\n";	
	#push(@thrs, threads->create(sub{system $tcmd;}));
	

	#phyml
	if ($doPhym){
		$tcmd = "$phymlBin --quiet -m GTR --no_memory_check -d nt -f m -v e -o tlr --nclasses 4 -b 2 -a e -i $multAli.ph > $nwkFile\n";
		push(@thrs, threads->create(sub{system $tcmd;}));
		$tcmd = "$phymlBin --quiet -m GTR --no_memory_check -d nt -f m -v e -o tlr --nclasses 4 -b 2 -a e -i $multAliSyn.ph > $nwkFile2\n";
		push(@thrs, threads->create(sub{system $tcmd;}));
	}
	for (my $t=0;$t<@thrs;$t++){
		my $state = $thrs[$t]->join();
		if ($state){die "Thread $t exited with state $state\nSomething went wrong with RaxML\n";}
	}
	system "rm  $multAli.ph $multAliSyn.ph $raxD/*_info* $raxD/*_log* $raxD/*_parsimony* $raxD/*_result*";
	system "mv $raxD/RAxML_fastTreeSH_Support.$raxFile2 $raxD/RXML_sym.nwk";
	system "mv $raxD/RAxML_fastTreeSH_Support.$raxFile $raxD/RXML_allsites.nwk";
	
	#die "$distTree_scr -d -a --dist-output $raxD/distance.syn.txt $raxD/RXML_sym.nwk\n";

	
	###################### ETE ######################3
} else {
	$cmd = "ete3 build -n $fnFna -a $aaFna -w clustalo_default-none-none-none  -m sptree_raxml_all --cpu $ncore -o $outD/tree --clearall --nt-switch 0.0 --noimg  --tools-dir /g/bork3/home/hildebra/bin/ete/ext_apps-latest"; #--no-seq-checks
	$cmd .= " --cogs $cogCats" unless ($cogCats eq "");
	print "Running tree analysis ..";
	system $cmd . "> $outD/tree/ETE.log";
	print " Done.\n$outD/tree\n";
}

print "All done \n\n";
exit(0);













##########################################################################################
##########################################################################################
sub mergeMSAs($ $ $ $){
	my ($MSAsAr,$samplesHr,$multAliF,$del) = @_;
	my @MSAs = @{$MSAsAr}; my %samples = %{$samplesHr};
	my %bigMSAFAA;foreach my $sm (keys %samples){$bigMSAFAA{$sm} ="";}
	foreach my $MSAf (@MSAs){
		#print $MSAf."\n"; 
		my $hit =0; my $miss =0;
		my $hr = readFasta($MSAf); my %MFAA = %{$hr};
		system "rm $MSAf" if ($del);
		my @Mkeys = keys %MFAA;
		#die "$Mkeys[0]\n";
		
		my @spl2 = split /_/,$Mkeys[0]; my $gcat = $spl2[1];
		my $len = length( $MFAA{$Mkeys[0]} );
		#die $len;
		foreach my $sm (keys %samples){
			my $curK = $sm."_".$gcat; #print $curK. " ";
			#print "$MFAA{$curK}\n";
			if (exists $MFAA{$curK}){
				$bigMSAFAA{$sm} .= $MFAA{$curK}; $hit++;
			} else {
				$bigMSAFAA{$sm} .= "-"x$len; $miss++;#print "nooooooo ";
			}
		}
		
		#die "$hit - $miss\n";
	}
	#filter part - count "-" in each seq
	
	my @ksMSAFAA = keys %bigMSAFAA;
	my $iniSeqNum = @ksMSAFAA; my $remSeqNum = 0;
	my %ntCnts; my $maxNtCnt=1;
	foreach my $kk (@ksMSAFAA){
		my $strCpy = $bigMSAFAA{$kk};
		my $num1 = $strCpy =~ tr/[\-N]//;
		#$num1++ while ($bigMSAFAA{$kk} !~ m/-/g);
		#now get actual num of bps
		$num1 = length($strCpy) - $num1;
		if ($num1 > $maxNtCnt){$maxNtCnt = $num1;}
		$ntCnts{$kk} = $num1;
	}
	foreach my $kk (@ksMSAFAA){
		my $num1 = $ntCnts{$kk};
		if ( ( ($num1 / $maxNtCnt ) < $ntFrac) || ($num1 < $ntCnt ) ){
			delete $bigMSAFAA{$kk}; $remSeqNum++; 
			print "$kk $num1  ".int($num1*1000 / ($maxNtCnt) )/1000 ." $ntFrac \n";
		}
		#print "$num1  $kk \n";#$bigMSAFAA{$kk}\n\n"; last;
	}
	open O,">$multAliF" or die "Can't open MSA outfile $multAliF\n";
	foreach my $kk (keys %bigMSAFAA){
		print O ">$kk\n$bigMSAFAA{$kk}\n";
	}
	close O;
	#die "$multAliF\n";
	print "Removed $remSeqNum of $iniSeqNum sequences\n";
}
sub convertMultAli2NT($ $ $){
	my ($inMSA,$NTs,$outMSA) = @_;
	my $tmpMSA=0;
	if ($inMSA eq $outMSA){$outMSA .= ".tmp"; $tmpMSA=1;}
	my $cmd = "";
	#"$pal2nal $inMSA $NTs -output fasta -nostderr -codontable 11 > $outMSA\n";
	
	#$cmd = "$trimalBin -in $inMSA -out $outMSA -backtrans $NTs -keepheader -keepseqs -noallgaps -automated1 -ignorestopcodon\n";
	$cmd = "$trimalBin -in $inMSA -out $outMSA -backtrans $NTs -keepheader -ignorestopcodon  -gt 0.1 -cons 60\n";
	#die "$inMSA,$NTs,$outMSA\n";
	#my $hr1= readFasta($inMSA);
	#my %MSA = %{$hr1};
	#$hr1= readFasta($NTs);
	#my %NTs = %{$hr1};
	if ($tmpMSA){$cmd .= "rm $inMSA;mv $outMSA $inMSA;\n";}
	#print $cmd;
	die "Can't execute $cmd\n" if (system $cmd);
}

sub synPosOnlyAA($ $){#only leaves "constant" AA positions in MSA file.. 
#stupid, don't know if pal2nal can handle this.. prob not
	my ($inMSA,$outMSA) = @_;
	print "Syn";
	my $hr = readFasta($inMSA); my %FNA = %{$hr};
	my @aSeq = keys %FNA;
	my $len = length ($FNA{$aSeq[0]});
	for (my $i=0; $i< $len; $i+=3){
		my $cod = substr $FNA{$aSeq[0]},$i,3;
		my $iniAA = "A";
		for (my $j=1;$j<@aSeq;$j++){
		}
	}
	print " only\n";

}

sub synPosOnly($ $ $ $){#not finished, I used the AA version instead
	my ($inMSA,$inAAMSA,$outMSA, $ffold) = @_;
	#print "Syn NT";
	my %convertor = (
    'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S',    # Serine
    'TTC' => 'F', 'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L', 'TTG' => 'L',    # Leucine
    'TAC' => 'Y',  'TAT' => 'Y',    # Tyrosine
    'TAA' => '*', 'TAG' => '*', 'TGA' => '*',    # Stop
    'TGC' => 'C', 'TGT' => 'C',    # Cysteine   
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L',    # Leucine
    'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',    # Proline
    'CAC' => 'H', 'CAT' => 'H',    # Histidine
    'CAA' => 'Q', 'CAG' => 'Q',    # Glutamine
    'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R',    # Arginine
    'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',    # Threonine
    'AAC' => 'N','AAT' => 'N',    # Asparagine
    'AAA' => 'K', 'AAG' => 'K',    # Lysine
    'AGC' => 'S', 'AGT' => 'S',    # Serine
    'AGA' => 'R','AGG' => 'R',    # Arginine
    'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',    # Valine
    'GCA' => 'A','GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',    # Alanine
    'GAC' => 'D', 'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E', 'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G','GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',    # Glycine
    );
	my %ffd;
	if ($ffold){ #calc 4fold deg codons in advance to real data
		foreach my $k (keys %convertor){
			my $subk = $k; my $iniAA = $convertor{$subk} ;
			my $cnt=0;
			foreach my $sNT ( ("A","T","G","C") ){
				
				substr ($subk,2,1) = $sNT;
				#print $subk ." " ;
				$cnt++ if ($convertor{$subk} eq $iniAA);
				
			}
#			if( $cnt ==4){ $ffd{$k} = 4;
#			} else {$ffd{$k} = 1;}
			if( $cnt ==4){ $ffd{$iniAA} = 4;
			} else {$ffd{$iniAA} = 1;}
		}
	}

	#assumes correct 3 frame for all sequences in inMSA
	my $hr = readFasta($inMSA); my %FNA = %{$hr};
	#my %FAA;
	#if (0  || !$ffold){
	#	$hr = readFasta($inAAMSA); %FAA = %{$hr};
	#}
	#print "$inMSA\n$inAAMSA\n$outMSA\n";
	my @aSeq = keys %FNA; my %outFNA;
	for (my $j=0;$j<@aSeq;$j++){$outFNA{$aSeq[$j]}="";}
	my $len = length ($FNA{$aSeq[0]});
	my $nsyn=0;my $syn=0;
	for (my $i=0; $i< $len; $i+=3){ #goes over every position
		my $j =0;
		my $iniAA = "-";
		while (1){ #check for first informative position
			my $iniCodon = substr $FNA{$aSeq[$j]},$i,3;
			if ($iniCodon =~ m/---/ || $iniCodon =~ m/N/){$j++; last if ($j >= @aSeq); next;}
			die "error: $iniCodon\n" if ($iniCodon =~ m/-/); #should not happen
			$iniAA = $convertor{$iniCodon};#substr $FAA{$aSeq[0]},$i,1; 
			last;
		}
		#die "$iniAA\n";
		my $isSame = 1;
		next unless (!$ffold || $ffd{$iniAA} == 4);
	#print $i." $iniAA ";
		for (;$j<@aSeq;$j++){
			my $newCodon = substr $FNA{$aSeq[$j]},$i,3;
			my $newAA = "-";
			if ($newCodon !~ m/-/ && $newCodon =~ m/[ACTG]{3}/i){
				die "Unkown AA $newCodon\n" unless (exists $convertor{$newCodon} );
				$newAA = $convertor{$newCodon} ; # substr $FAA{$aSeq[$j]},$i,1;
			} elsif ($newCodon =~ m/---/ || $newCodon =~ m/[NWYRSKMDVHB]/){
			} else {
				die "newCodon wrong $newCodon\n" ;
			}
			if ($iniAA ne $newAA && $newAA ne "-"){
				$isSame =0; last;
			}
		}
		if ($isSame){#add nts to file
			for (my $j=0;$j<@aSeq;$j++){
				if ($ffold){
					$outFNA{$aSeq[$j]} .= substr $FNA{$aSeq[$j]},($i)+2,1;
				} else {
					$outFNA{$aSeq[$j]} .= substr $FNA{$aSeq[$j]},$i,3;
				}
				#print substr $FNA{$aSeq[$j]},$i*3,3 . " ";
			}
			$syn++;
		} else {$nsyn++;}
	}
	#die $inMSA."\n";
	open O ,">$outMSA" or die "Can't open outMSA $outMSA\n";
	for (my $j=0;$j<@aSeq;$j++){
		print O ">$aSeq[$j]\n$outFNA{$aSeq[$j]}\n";
	}
	close O;
	$aSeq[0] =~ m/^.*_(.*)$/;
	#die "$outMSA\n";
	print "$1 ($syn / $nsyn) ".@aSeq." seqs \n";
	#print " only\n";
	#print "\n";
}






