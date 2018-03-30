package Mods::phyloTools;
use warnings;
#use Cwd 'abs_path';
use strict;
#use List::MoreUtils 'first_index'; 

#use Mods::GenoMetaAss qw(qsubSystem);

use Exporter qw(import);
our @EXPORT_OK = qw(runRaxML prep40MGgenomes prepNOGSETgenomes getE100 getGenoGenes getFMG renameFMGs 
			runFasttree runQItree fixHDs4Phylo getGenoName);
use Mods::GenoMetaAss qw(systemW readFasta renameFastHD gzipwrite gzipopen);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::FuncTools qw(assignFuncPerGene);

sub getGenoGenes;


sub fixHDs4Phylo ($){
	#routine to check that headers of fastas don't contain ":", ",", ")", "(", ";", "]", "[", "'"
	my ($inF) = @_;
	my $reqFix=0;
	if ($inF eq ""){return "";}
	my $hr = readFasta($inF,1); my %FAA = %{$hr};
	foreach my $hd (keys %FAA){
		if ($hd =~ m//){
			$reqFix=1;last;
		}
	}
	my $outF = $inF;
	if ($reqFix){
		$outF .= ".fix";
		print "Fixing headers in input file (to $outF)\n";
		my %newHDs;
		open O,">$outF";
		foreach my $hd (keys %FAA){
			my $hd2 = substr $hd,0,40; #cut to raxml length
			$hd2 =~ s/[:,\}\{;\]\[']/|/g;
			$newHDs{$hd2} ++;
			$hd2 .= $newHDs{$hd2};
			print O ">$hd2\n$FAA{$hd}\n";
		}
	}
	return $outF;
}

sub runQItree{
	my ($inMSA,$treeOut,$ncore,$outgr,$bootStrap,$useAA,$fast,$autoModel) = @_;
	my $iqTree  = getProgPaths("iqtree");
	$treeOut =~ s/\.nwk$//;
	my $treNM = "IQtree";
	my $cmd = "$iqTree -s $inMSA -nt $ncore -pre $treeOut ";
	$cmd .= "-o $outgr " unless ($outgr eq "" && $outgr !~ m/,/);
	unless ($fast == 0){$cmd .= "-fast "; print "IQtree - fast\n"; $treNM .= "_fast";}
	if ($autoModel){$treNM .= "_autoMOD";}
	if ($useAA){
		if ($autoModel){
			$cmd .= "-m TEST  "; 
		} else{
			$cmd .= "-m LG+F+G "; #needs to be HKY for nts
		}
	} else {
		if ($autoModel){
			$cmd .= "-m TEST ";#-mset HKY,HKY+F,HKY+F+I,HKY+F+I+G4,JC,F81,K2P,K3P,K81uf,GTR "; 
		} else {
			$cmd .= "-m HKY+F+G "; 
		}
	}
	if ($bootStrap >0){
		if ($bootStrap < 1000){
			print "standard non parametric bootstrap ($bootStrap). Use >1000 bootstraps to do ultrafast bootstrap\n";
			die  "normal bootstrap requires >= 100 iterations ($bootStrap given)\n" if ($bootStrap <100);
			$cmd .= "-b $bootStrap " ;
		} else {
			print " ultrafast bootstrap $bootStrap\n";
			$cmd .= "-bb $bootStrap " ;
		}
		#also consider -b >=100 for std bootstrap
	} else {
		$cmd .= "-alrt 1000 ";
	}
	$cmd .= "";
	#die $cmd;
	systemW $cmd;
	$treNM .= ".nwk";
	#"mv $treeOut/IQtree_fast_allsites.treefile $treeOut/$treNM";
}
sub runFasttree{
	my ($inMSA,$treeOut,$isAA,$ncore) = @_;
	my $fsttreeBin  = getProgPaths("fasttree");
	my $ntFlag = "";
	$ntFlag = "-nt -gtr" if ($isAA);
	my $cmd = "$fsttreeBin $ntFlag $inMSA > $treeOut\n";
	systemW $cmd;
}


sub getGenoName($){my $GenomeN = $_[0];$GenomeN =~ s/\.f[n]?a$//; $GenomeN =~ s/^.*\///;  $GenomeN =~ s/-/_/; return $GenomeN;}

#routine to format marker genes with naming etc to work with buildTree script
sub prep40MGgenomes{ 
	#\@refGenos,$rewrite,$ncore,$fnFna1,$aaFna1,$cogCats1)
	my @refGenos = @{$_[0]};my $finalD=$_[1]; my $tag= $_[2];
	my $rewrite = $_[3]; my $ncore=$_[4];
	my ($fnFna1,$aaFna1,$cogCats1) = ("","","");
	my $outGrp = "";
	$outGrp = $_[5] if (@_ > 5);
	if (@_ > 6){$fnFna1= $_[6];$aaFna1= $_[7];$cogCats1= $_[8];}
	print "Using outgroup $outGrp\n" if ($outGrp ne "");
	
	my $buildTreeScr = getProgPaths("buildTree_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/buildTree.pl";
	my $FMGd = getProgPaths("FMGdir");#"/g/bork5/hildebra/bin/fetchMG/";
	my $FMGrwkScr = getProgPaths("FMGrwk_scr");

	#system "$renameCtgScr $refGenos $GenomeN";
	system "mkdir -p $finalD" unless (-d $finalD);
	my $fnFna = $fnFna1; if ($fnFna eq ""){$fnFna = $finalD."$tag.fna";} else {$fnFna =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	my $aaFna = $aaFna1; if ($aaFna eq ""){$aaFna = $finalD."$tag.faa";} else {$aaFna =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	my $cogCats = $cogCats1; if ($cogCats eq ""){$cogCats = $finalD."$tag.cats.txt";} else {$cogCats =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	if (-e $fnFna1){system("cat $fnFna1 > $fnFna");} else {system "rm -f $fnFna";}#shortFNAhd($fnFna);}
	if (-e $aaFna1){system("cat $aaFna1 > $aaFna");} else {system "rm -f $aaFna";}
	my %COGgenes; my %allGenomes;
	#assumes ref genome, no genes called
	for (my $i=0;$i<@refGenos;$i++){
		my $refG = $refGenos[$i];
		next if ($refG eq "");
		my $GenomeN = getGenoName($refGenos[$i]);
		
		if (exists($allGenomes{$GenomeN})){
			print "Genome $GenomeN already exists; assumming double entry and skipping genome\n";next;
		} else {
			$allGenomes{$GenomeN} = 1;
		}
		#die "$GenomeN\n";
		my ($ntGenes,$proteins) = getGenoGenes($refG);
		#my $proteins = $refG; my $ntGenes = $refG; $proteins =~ s/\.[^\.]+$/\.genes\.faa/;
		#$ntGenes =~ s/\.[^\.]+$/\.genes\.fna/;
		#die "Can't find ref genome $refG\n" unless (-e $refG || (-e $ntGenes && -e $proteins));
		#my $prodigal_cmd .= "$prodigalBin -i $refG -a $proteins -d $ntGenes -f gff -p single > /dev/null\n";
		#die $prodigal_cmd."\n$rewrite || !-e $proteins || !-e $ntGenes\n";
		#system $prodigal_cmd if ($rewrite || !-e $proteins || !-e $ntGenes);
		#fetch mg
		my $cmd ;
		my $GenomeDir = $proteins; $GenomeDir =~ s/[^\/]+$//;
		my $outDFMG = $GenomeDir."FMGs$GenomeN/";
		#getFMG($outDFMG,$ntGenes,$proteins,$ncore,$rewrite+3);
		if (!-s "$outDFMG/COG0012.faa" && !-s "$outDFMG/COG0016.faa"){
			$cmd = "perl $FMGd/fetchMG.pl -m extraction -o $outDFMG -l $FMGd/lib -t $ncore -d $ntGenes $proteins"; # -x $FMGd/bin  -b $FMGd/lib/MG_BitScoreCutoffs.uncalibrated.txt
			system $cmd;
			system "rm  -rf $GenomeDir/*.cidx $outDFMG/temp $outDFMG/hmmResults";
		}
		my $hr = renameFMGs($outDFMG,$GenomeN,"allFMG",1);
		#system "$FMGrwkScr $outDFMG" ;
		my %COGg = %{$hr};
		foreach my $k (keys %COGg){
			push (@{$COGgenes{$k}},$COGg{$k});
		}
		#open I,"<$outDFMG/FMGids.txt";
		#while (<I>){
		#	chomp;my @spl = split /\s/;
		#	if (!exists($COGgenes{$spl[1]})){$COGgenes{$spl[1]}=$spl[0];}else{$COGgenes{$spl[1]}.= "\t".$spl[0];  }
		#}
		#close I;
		#now add this info to the existing ete input files
		system("cat $outDFMG/*.fna >> $fnFna");#shortFNAhd($fnFna);}
		system("cat $outDFMG/*.faa >> $aaFna") ;#shortFNAhd($aaFna);
	}
	#open genecat file 
	my $cogcatstr = "";my $catCnt=0;
	if (-e $cogCats1){
		open I,"<$cogCats1";
		while(<I>){
			chomp;my @spl = split /\t/;
			$spl[0] =~ m/_([^_]+)$/;
			die "Can't find category $1 in ref gene set\n" unless (exists($COGgenes{$1}));
			push (@spl,@{$COGgenes{$1}});
			$cogcatstr .= join ("\t",@spl) . "\n";
			$catCnt++;
		}
		close I;
	} else { #make anew
		foreach my $k (keys %COGgenes){
			$cogcatstr .= join ("\t",@{$COGgenes{$k}}) . "\n";
			$catCnt++;
		}
	}
	open O,">$cogCats";
	print O $cogcatstr;
	close O;
	print "Found $catCnt / ".scalar(keys %COGgenes)." gene cats\n";
	#and run tree building
	my $ogrpStr="";
	if ($outGrp ne ""){
		$ogrpStr = "-outgroup $outGrp";
	}
	my $cmd = "$buildTreeScr -fna $fnFna -aa $aaFna -cats $cogCats -outD $finalD -cores $ncore -NTfilt 0.15 -bootstrap 100  $ogrpStr > $finalD/tree_build.log";
	return $cmd;
}
#similar to prep40MGgenomes, but extended to more genes (as given by NOG annotation)
sub prepNOGSETgenomes{ 
	#\@refGenos,$rewrite,$ncore,$fnFna1,$aaFna1,$cogCats1)
	my @refGenos = @{$_[0]};my @NOGs = @{$_[1]};
	my $finalD=$_[2]; my $tag= $_[3];
	my $rewrite = $_[4]; my $ncore=$_[5];
	my ($fnFna1,$aaFna1,$cogCats1) = ("","","");
	my $outGrp = "";
	$outGrp = $_[6] if (@_ > 6);
	#these files will be added on top of the newly derrived genes
	if (@_ > 7){$fnFna1= $_[7];$aaFna1= $_[8];$cogCats1= $_[8];}
	print "Collating ".@NOGs." NOG genes for tree building\n";
	print "Using outgroup $outGrp\n" if ($outGrp ne "");
	
	my $buildTreeScr = getProgPaths("buildTree_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/buildTree.pl";
	my $FMGd = getProgPaths("FMGdir");#"/g/bork5/hildebra/bin/fetchMG/";
	my $FMGrwkScr = getProgPaths("FMGrwk_scr");

	#system "$renameCtgScr $refGenos $GenomeN";
	system "mkdir -p $finalD" unless (-d $finalD);
	my $tmpP = $finalD."tmp/";
	system "mkdir -p $tmpP" unless (-d $tmpP);
	
	my $fnFna = $fnFna1; if ($fnFna eq ""){$fnFna = $finalD."$tag.fna";} else {$fnFna =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	my $aaFna = $aaFna1; if ($aaFna eq ""){$aaFna = $finalD."$tag.faa";} else {$aaFna =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	my $cogCats = $cogCats1; if ($cogCats eq ""){$cogCats = $finalD."$tag.cats.txt";} else {$cogCats =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	if ( -e $fnFna1){system("cat $fnFna1 > $fnFna");} else {system "rm -f $fnFna";}#shortFNAhd($fnFna);}
	if ( -e $aaFna1){system("cat $aaFna1 > $aaFna");} else {system "rm -f $aaFna";}

	my %COGgenes; my %allGenomes;
	

	#assumes ref genome, no genes called
	for (my $i=0;$i<@refGenos;$i++){
		my $refG = $refGenos[$i];
		next if ($refG eq "");
		my $GenomeN = getGenoName($refGenos[$i]);
		
		if (exists($allGenomes{$GenomeN})){
			print "Genome $GenomeN already exists; assumming double entry and skipping genome\n";next;
		} else {
			$allGenomes{$GenomeN} = 1;
		}
		#die "$GenomeN\n";
		my ($ntGenes,$proteins) = getGenoGenes($refG);
		#die "$proteins\n";
		#functional annotation via NOG
		my $GenomeDir = $proteins; $GenomeDir =~ s/[^\/]+$//;
		my $outDNOG = $GenomeDir."NOGs$GenomeN/";
		my %optsDia = (eval=>1e-12,percID=>35,minPercSbjCov=>0.6,fastaSplits=>1,ncore=>$ncore,splitPath=>$tmpP,bacNOG=>1);
		my ($outF,$jdep) = assignFuncPerGene($proteins,$outDNOG,$tmpP,"NOG",\%optsDia);
		my $tarAnno = "${outF}geneAss.gz";
		#print "XX $tarAnno\n";
		
		my %COG2Gn;
		my ($I,$OK) = gzipopen($tarAnno,"Func anno",1); 
		while (<$I>){
			chomp;
			my @spl = split /\t/;
			$COG2Gn{$spl[1]} = $spl[0];
			#print;
		}
		close $I;
		#print "$tarAnno\n";
		#next; #DEBUG
		#now load AA and NT file, and rename genes to fit
		my $hr = readFasta($ntGenes,1);
		my %FAS = %{$hr};
		$hr = readFasta($proteins,1);
		my %FASAA = %{$hr};
		#print "Read NT AA\n";
		my %COGg ;
		open O,">>$fnFna" or die "Can't open $fnFna\n";
		open O2,">>$aaFna" or die "Can't open $aaFna\n";
		my $hitcnt=0;
		my $cnt=0;
		foreach my $tag (@NOGs){
			my $newName = "${GenomeN}_${tag}";
			#print $tag."\n" if ($cnt<5);; $cnt++;
			#if ($tag =~ m//)
			if (exists($COG2Gn{$tag})){
				#print "$COG2Gn{$tag}\n";
				die "Can't find $COG2Gn{$tag} in nt's\n" unless(exists $FAS{$COG2Gn{$tag}});
				die "Can't find $COG2Gn{$tag} in AA's\n" unless(exists $FASAA{$COG2Gn{$tag}});
				print O ">$newName\n$FAS{$COG2Gn{$tag}}\n";
				print O2 ">$newName\n$FASAA{$COG2Gn{$tag}}\n";
				$COGg{$tag} = $newName;
				$hitcnt++;
			}
		}
		close O;
		close O2;
		print "Found $hitcnt core genes in genome $GenomeN\n";
		
		foreach my $k (keys %COGg){
			push (@{$COGgenes{$k}},$COGg{$k});
		}
		#die "$fnFna\n";
		#last if ($i==10);
	}

	#open genecat file 
	my $cogcatstr = "";my $catCnt=0;
	if (-e $cogCats1){
		open I,"<$cogCats1" or die "Can't open $cogCats1\n";
		
		while(<I>){
			chomp;my @spl = split /\t/;
			$spl[0] =~ m/_([^_]+)$/;
			die "Can't find category $1 in ref gene set\n" unless (exists($COGgenes{$1}));
			push (@spl,@{$COGgenes{$1}});
			$cogcatstr .= join ("\t",@spl) . "\n";
			$catCnt++;
		}
		close I;
	} else { #make anew
		foreach my $k (keys %COGgenes){
			$cogcatstr .= join ("\t",@{$COGgenes{$k}}) . "\n";
			$catCnt++;
		}
	}
	open O,">$cogCats";
	print O $cogcatstr;
	close O;
	print "Found $catCnt / ".scalar(keys %COGgenes)." gene cats\n";
	#and run tree building
	#die "$cogCats\n";
	my $ogrpStr="";
	if ($outGrp ne ""){
		$ogrpStr = "-outgroup $outGrp";
	}
	my $cmd = "$buildTreeScr -fna $fnFna -aa $aaFna -cats $cogCats -outD $finalD -cores $ncore -NTfilt 0.15 -bootstrap 100  $ogrpStr > $finalD/tree_build.log";
	system "rm -r $tmpP" ;
	return $cmd;
}


sub runRaxML{
	my ($mAli,$bootStrap,$outgroup,$outTree,$ncore) = @_;
	my $cont = 0; my $useNT=1; my $useAA=0;
	$cont = $_[5] if (@_ > 5);
	$useNT = $_[6] if (@_ > 6);
	$useAA = !$useNT;#useless flag atm
	#die "$cont\n";
	my ($raTmpF,$raTmpF2,$raTmpF3) = ("RXMtmp","R2XM","RX3M"); #my $raxFile = "RXMall";
#	my $raxFile2 = "RXMsyn";
	my $raxmlBin = getProgPaths("raxml");
	my $trDist = getProgPaths("treeDistScr");
	my $raxD = $outTree;$raxD=~m/([^\/]+$)/; my $trNm = $1;$trNm=~ s/\..*$// if ($trNm =~ m/\./);$raxD =~ s/[^\/]+$/${trNm}tmp\//;
	my $raxLogD= $outTree; $raxLogD =~ s/[^\/]+$/${trNm}Logs\//;
	#die "$raxD";
	system "rm -rf $raxD;" unless ($cont);
	system "mkdir -p $raxD";
	system "mkdir -p  $raxLogD" unless (-d $raxLogD);
	my $outGrpRXML = "";
	if ($outgroup ne ""){$outGrpRXML = "-o $outgroup";}
	my $raxDef = " --silent -m GTRGAMMA -p 312413 ";
	if (!$useNT){
		$raxDef = " --silent -m PROTGAMMALG -p 312413 ";
	}
	
	#raxml - on all sites
	#die "$outTree\n";
	if (!$cont || !-e $outTree){
		my $tcmd="";
		$tcmd =  "$raxmlBin -T$ncore -f d -s $mAli $raxDef -n $raTmpF -w $raxD $outGrpRXML > $raxLogD/ini.log\n";
		#die $tcmd."\n";
		my $expTree1 = "$raxD/RAxML_bestTree.$raTmpF";
		$expTree1 = "$raxD/RAxML_result.$raTmpF" if ($useAA);
		if (!-e $expTree1){print "Calculating ML tree..\n" ; system "rm -f $raxD/*$raTmpF";system $tcmd;}
		if (!-e $expTree1 && !$useAA){#failed.. prob optimization problem, use other tree instead
			print "Calculating ML tree with GTRGAMMAI model..\n" ;
			$raxDef = " --silent -m GTRGAMMAI -p 3512413 ";	system "rm -f $raxD/*$raTmpF";
			systemW "$raxmlBin -T$ncore -f d -s $mAli $raxDef -n $raTmpF -w $raxD $outGrpRXML > $raxLogD/optimized.log\n";
		}elsif (!-e $expTree1){
			die "Failed $tcmd\n";
		}
		#decide which support vals to calc
		if ($bootStrap==0){
			
			print "Calculating SH support tree..\n";
			systemW "$raxmlBin -T$ncore -f J -s $mAli $raxDef -n $raTmpF3 -w $raxD -t $expTree1 $outGrpRXML > $raxLogD/optimized.log\n";
			system "mv $raxD/RAxML_fastTreeSH_Support.$raTmpF3 $outTree";
		} else {
			my $done=0; my $partBoots = "";
			if (-e "$raxD/RAxML_bootstrap.$raTmpF2" || -e "$raxD/RAxML_bootstrap.prev.x"){
				if (-e "$raxD/RAxML_bootstrap.$raTmpF2"){
					my $tmp = `wc -l $raxD/RAxML_bootstrap.$raTmpF2`; $tmp =~ m/^(\d+) /; $done += $1;
				}
				if (-e "$raxD/RAxML_bootstrap.prev.x"){
					my $tmp = `wc -l $raxD/RAxML_bootstrap.prev.x`; $tmp =~ m/^(\d+) /; $done += $1;
					$partBoots="$raxD/RAxML_bootstrap.prev.x";
				}
				if ($done < $bootStrap){system "cat $raxD/RAxML_bootstrap.$raTmpF2 >> $raxD/RAxML_bootstrap.prev.x;rm $raxD/RAxML_bootstrap.$raTmpF2" if (-e "$raxD/RAxML_bootstrap.$raTmpF2"); }
				
			}
			#die "$done < $bootStrap\n";
			if (!$cont || $done < $bootStrap){
				system "rm -f $raxD/*$raTmpF2"; 
				#if ($partBoots ne ""){system "cat $partBoots >> $raxD/RAxML_bootstrap.$raTmpF2";}
			}
			if ($done >= $bootStrap){
				system "mv $partBoots  $raxD/RAxML_bootstrap.$raTmpF2";
			}
			$bootStrap = $bootStrap - $done ;
			if (!-e "$raxD/RAxML_bootstrap.$raTmpF2"){
				my $bootCmd.= "-N $bootStrap -b ". int(rand(10000));
				print "Calculating bootstrap tree with $bootStrap bootstraps..\n";
				systemW  "$raxmlBin -T$ncore -f d  -s $mAli $raxDef -n $raTmpF2 -w $raxD $bootCmd $outGrpRXML > $raxLogD/ini.log\n";
				if ($partBoots ne ""){system "cat $partBoots >> $raxD/RAxML_bootstrap.$raTmpF2";}
			}
			#die;
			system "rm -f $raxD/*$raTmpF3";
			systemW  "$raxmlBin -T$ncore -f b  -s $mAli $raxDef -n $raTmpF3 -w $raxD -t $raxD/RAxML_bestTree.$raTmpF -z $raxD/RAxML_bootstrap.$raTmpF2 $outGrpRXML > $raxLogD/ini.log\n";
			system "mv $raxD/RAxML_bipartitions.$raTmpF3 $outTree";
		
		}
	} else {
		print "Tree already present\n";
	}
	#and create distance matrix
#	print "Calculating ML matrix from tree..\n";
#	system "rm -f $raxD/*all"; 
#	systemW  "$raxmlBin -T$ncore -f x -s $mAli $raxDef -n all -w $raxD -t $raxD/RAxML_bestTree.$raTmpF $outGrpRXML > $raxLogD/.log\n";
	#die "$outDist\n";
#	system "mv $raxD/RAxML_distances.all $outDist";
	
	print "Calculating tree distance matrix from tree..\n";
	my $outDist = $outTree; $outDist =~ s/\.[^\.]+$/\.dist/;
	#system "rm -f $raxD/*all"; 
	#this perl script has some errors!!
	#system  "$trDist $outTree > $outDist\n";
#	system "mv $raxD/RAxML_distances.all $outDist";

	
	#clean up
	system "rm -rf $raxD";
}

sub getGenoGenes{
	my $refG = $_[0];
	my $rewrite =0;
	$rewrite = $_[1] if (@_ >= 2);
	my $prodigalBin = getProgPaths("prodigal");
	my $proteins = $refG; my $ntGenes = $refG; $proteins =~ s/\.[^\.]+$/\.genes\.faa/;
	$ntGenes =~ s/\.[^\.]+$/\.genes\.fna/;
	die "Can't find ref genome $refG\n" unless (-e $refG || (-e $ntGenes && -e $proteins));
	my $prodigal_cmd .= "$prodigalBin -i $refG -a $proteins -d $ntGenes -f gff -p single > /dev/null\n";
	#die $prodigal_cmd."\n$rewrite || !-e $proteins || !-e $ntGenes\n";
	system $prodigal_cmd if ($rewrite || !-e $proteins || !-e $ntGenes);
	return ($ntGenes,$proteins);
}
sub getFMG{
	my ($oDess,$proteins,$genesNT) = @_;
	my $redo=0;my $ncore=1;my $rename="";
	$ncore = $_[3] if (@_ >= 4);
	$redo = $_[4] if (@_ >= 5);
	$rename = $_[5] if (@_ >= 6);
	
	if ($proteins eq ""){die "getFMG:: need protein file\n$proteins\n";}
	if ($oDess eq ""){
		$proteins=~m/^(.*)\/([^\/]+)$/;
		$oDess = $1;
		if ($2 =~ m/(.*)\b\.genes\b?\.faa$/){
			$oDess .= "/FMGs$1/";
		} else { die "getFMG:: protein name needs to end on .faa\n$2\n";}
	}
	#die $proteins."\n".$oDess."\n";
	my $GenomeDir = $proteins; $GenomeDir =~ s/[^\/]+$//;
	system "mkdir -p $oDess" unless (-d $oDess);
	#die $cmd;
	if ( $redo  || (!-s "$oDess/COG0012.faa" && !-s "$oDess/COG0016.faa")){
		my $FMGd = getProgPaths("FMGdir");#"/g/bork5/hildebra/bin/fetchMG/";
		my $FMGrwkScr = getProgPaths("FMGrwk_scr");
		my $cmd = "";
		$cmd = "perl $FMGd/fetchMG.pl -m extraction -o $oDess -l $FMGd/lib -t $ncore -d $genesNT $proteins "; # -x $FMGd/bin  -b $FMGd/lib/MG_BitScoreCutoffs.uncalibrated.txt
		#print $cmd;
		systemW $cmd;
		system "rm  -rf $GenomeDir/*.cidx $oDess/temp $oDess/hmmResults";
		system "$FMGrwkScr $oDess" if ($redo);
	}
	#die "$redo\n";
	return $oDess;
}
sub renameFMGs{
	my ($oDess,$rename,$outN,$renameShrt) = @_;
	#create an "allFMG" file with the correct fasta header names
	my %catGe;
	my $allFine=0;
	opendir(DIR, $oDess) or die "Can't find: $oDess\n";	
	my @fnas = sort ( grep { /COG\d+\.fna/  && -e "$oDess/$_"} readdir(DIR) );	close(DIR);
	my @fnaT=("","");my $aaT="";
	foreach my $nf (@fnas){
		$nf =~ m/(.*)\.fna/;
		my $cat = $1;
		for (my $j=0;$j<2;$j++){
			$nf =~ s/\.fna/\.faa/ if ($j==1);
			my $hdFnd=0;
			open I ,"<$oDess$nf";
			while (my $li=<I>){
				if ($li =~ m/^>/){ 
					if ($hdFnd){print "Found > 1 seq in $nf\n";last;}
					$hdFnd=1;
					my $newH = $rename."_". $cat ;
					$fnaT[$j] .=  ">".$newH  . "\n";
					$catGe{$cat} = $newH;
					if ($li =~ m/>$newH/){$allFine=1;}
				} else {
					$fnaT[$j] .= $li;
				}
			}
			close I;
			if (!$allFine && $renameShrt==1){
				open O ,">$oDess$nf"; print O $fnaT[$j]; close O;
				$fnaT[$j]=""; 
			}
		}
	}
	if ($renameShrt==0){
		print "renaming FMGs to name $rename\n";
		open O,">$oDess/$outN.fna" or die "Can't open renameFMGs outfile $oDess/$outN.fna\n";
		print O $fnaT[0]; close O;
		open O,">$oDess/$outN.faa" or die "Can't open renameFMGs outfile $oDess/$outN.faa\n";
		print O $fnaT[0]; close O;
	}
	return (\%catGe);
}
sub getE100($ $ $ $){
	my ($oDess,$proteins,$genesNT,$ncore) = @_;
	my $hmmbin = getProgPaths("hmmsearch");#"/g/bork5/hildebra/bin/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmsearch";
	my $essDB = getProgPaths("essentialHMM");#"/g/bork5/hildebra/bin/multi-metagenome-master/R.data.generation/essential.hmm";
	my $essEukDB = getProgPaths("essentialEUK");#"/g/bork3/home/hildebra/DB/HMMs/eukCore/eukCore.hmm"; #TODO
	systemW("mkdir -p $oDess");
	my $cmd = "";
	$cmd .= "$hmmbin --domtblout $oDess/assembly.hmm.orfs.txt -o $oDess/assembly.hmmsout.txt --cut_tc --cpu $ncore --notextw $essDB $proteins\n";
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
}
