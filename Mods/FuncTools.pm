package Mods::FuncTools;
use warnings;
use Cwd 'abs_path';
#use File::chdir;

use strict;
#use List::MoreUtils 'first_index'; 
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::TamocFunc qw( getSpecificDBpaths);
use Mods::GenoMetaAss qw( clenSplitFastas splitFastas  qsubSystem emptyQsubOpt systemW );

use Exporter qw(import);
our @EXPORT_OK = qw(assignFuncPerGene calc_modules);



#unifying script to diamond AA gene file against required DB & assign functions via the dia filter scripts of MATAFILER
#for now mostly diamond focused: NOG,KGM,KGE,KGB,ABRc,PTV,CZy,ACL,TCDB,PAB
#now-diamond: mp3 (pathogenic genes)
sub assignFuncPerGene{
	my ($query,$outD,$tmpD,$curDB) = @_;

	my $secCogBin = getProgPaths("secCogBin_scr");#MATAFILER script to filter raw mappings by diamond/blast etc
	my $diaBin = getProgPaths("diamond");#diamond binary
	my $prsMP3_scr = getProgPaths("mp3Prs");

	my $otpsHR = {};
	$otpsHR = $_[4] if (@_ > 4);
	my $doQsub=0; #if 0, then local execution
	my $QSBoptHR = {};
	my $qsubDir = $outD."qsubLOG/";;
	if (@_ > 5 && $_[5] != 0){
		$doQsub = 1;
		$QSBoptHR = $_[5] ;
		$qsubDir = $QSBoptHR->{qsubDir} if (exists($QSBoptHR->{qsubDir}) && $QSBoptHR->{qsubDir} ne "" );
	}
	my $localTmp = 0;#local tmp, requires different file copying
	#$localTmp = 1 if ($doQsub);
	my $fastaSplits = 1; #defaults
	my $ncore = 10;
	my %opts = %{$otpsHR};
	$ncore = $opts{ncore} if (exists($opts{ncore}));
	$fastaSplits = $opts{fastaSplits} if (exists($opts{fastaSplits}));
	$opts{splitPath} = $tmpD if (!exists($opts{splitPath}));
	$opts{redo} = 0 if (!exists($opts{redo}));
	$opts{keepSplits} = 0 if (!exists($opts{redo}));
	my $redo = $opts{redo};
	my $globalDiamondDependence = "";
	system "mkdir -p $outD" unless (-d $outD);
	if (!$opts{keepSplits}){
		clenSplitFastas($query,$opts{splitPath}."/");
	}
	
	
	
	#build DB
	my ($DBpath ,$refDB ,$shrtDB) = getSpecificDBpaths($curDB,1);
	if ($DBpath ne "" && !-e "$DBpath$refDB.db.dmnd"){
		my $DBcmd .= "$diaBin makedb --in $DBpath$refDB -d $DBpath$refDB.db -p $ncore\n";
		if ($doQsub){
			my ($jN, $tmpCmd) = qsubSystem($qsubDir."DiamondDBprep.sh",$DBcmd,$ncore,"2G","diaDB","","",1,[],$QSBoptHR);
			$globalDiamondDependence = $jN;
		} else {
			systemW $DBcmd;
		}
	}
	
	my @subFls ;
#	print "Splitting catalog into $fastaSplits\n" if ($fastaSplits>1);
	my @jdeps; my @allFiles;
	my $allAss = "$outD/DIAass_$shrtDB.srt.gz";
	my $tarAnno = "${allAss}geneAss.gz";
	system "rm -f $allAss" if ($redo);
	system "rm -f $tarAnno" if ($redo);
	my $calcDia = 1;$calcDia = 0 if (-e $allAss);
	my $interpDia = 1;$interpDia = 0 if (-e $tarAnno);
	#my $N = 20;
	my $jdep=""; my $qCmd = "";
	
	#die "!$calcDia && !$interpDia $curDB\n$allAss\n$tarAnno\n";
	return $allAss,$jdep if (!$calcDia && !$interpDia);

	
	#die "$calcDia\n$allAss\n";
	#alignment options (and defaults)
	$opts{eval} = 1e-7 if (!exists($opts{eval}));
	$opts{percID} = 25 if (!exists($opts{percID}));
	$opts{minAlignLen} = 60 if (!exists($opts{minAlignLen}));
	$opts{minBitScore} = 60 if (!exists($opts{minBitScore}));
	$opts{minPercSbjCov} = 0.3 if (!exists($opts{minPercSbjCov}));
	$opts{bacNOG} = 0 if (!exists($opts{bacNOG}));

	print "$query assigned to $curDB ($fastaSplits splits, $ncore cores)\n" if ($calcDia || $interpDia);

	
	if ($calcDia){
		print "Diamond pars: eval=$opts{eval}, percID=$opts{percID}, minAlLength=$opts{minAlignLen}, minBitScore=$opts{minBitScore}, minPercSbjCov=$opts{minPercSbjCov}}\n" if ($calcDia);
		my $ar = splitFastas($query,$fastaSplits,$opts{splitPath}."/");
		@subFls = @{$ar};
		for (my $i =0 ; $i< @subFls;$i++){
			my $cmd = "mkdir -p $tmpD\n";
			my $outF = "";
			if ($curDB eq "mp3"){
				my $mp3Dir = getProgPaths("mp3");
				#$CWD = $mp3Dir;
				chdir $mp3Dir;	
				$cmd .= "$mp3Dir/./mp3 $subFls[$i] 1 50 0.2\n";#1=complete genes,2=metag prots
				if ($localTmp){
					$subFls[$i] =~ m/\/([^\/]+$)/; my $fnm=$1;
					$outF = "$outD/$fnm.Hybrid.result.gz"; #this is sometimes shared, sometimes local tmp
					#$cmd .= "cp $subFls[$i].Hybrid.result $outF\n";
					$cmd .= "$prsMP3_scr $subFls[$i].Hybrid.result $outF\n";
					#$cmd .= "rm $tmpD\n";
				} else {
					$outF = "$subFls[$i].Hybrid.result.gz";
				}
			} else {
				#my $outF = "$GCd/DiaAssignment.sub.$i";
				if ($localTmp){
					$outF = "$outD/DiaAs.sub.$i.$shrtDB.gz";
				} else {
					$outF = "$tmpD/DiaAs.sub.$i.$shrtDB.gz";
				}
				$cmd .= "$diaBin blastp -f tab --compress 1 --quiet -t $tmpD -d $DBpath$refDB.db -q $subFls[$i] -k 5 -e 1e-5 --sensitive -o $outF -p $ncore\n";
				#$cmd = "$diaBin blastp -f tab --compress 1 --sensitive --quiet -d $eggDB.db -q $subFls[$i] -k 3 -e 0.001 -o $outF -p $ncore\n";
				#$cmd .= "$diaBin view -a $outF.tmp -o $outF -f tab\nrm $outF.tmp* $subFls[$i] \n";
				$cmd .= "rm -r $tmpD\n" if ($localTmp);
			}
			my $tmpCmd;
			system "rm -f $outF" if ($redo);
			if ($calcDia && !-e $outF){
						#die "$cmd\n";
				if ($doQsub){
					my ($jobName,$mptCmd) = qsubSystem($qsubDir."D$i$shrtDB.sh",$cmd,$ncore,"3G","D$shrtDB$i",$globalDiamondDependence,"",1,[],$QSBoptHR); #$jdep.";".
					push(@jdeps,$jobName);
				} else {
					systemW $cmd;
				}
				#print $qsubDir."Diamond.sh\n";
			}
			push(@allFiles,$outF);
			#if ($i==5){die;}
		}
	}
	#last job that converges all
	my $cmd = "";
	my $catcmd = "cat";
	if (!$calcDia){
	
	} elsif (@allFiles == 1){
		if ($allFiles[0] =~ m/\.gz$/){
			$cmd .= "mv $allFiles[0] $allAss\n";
		} else {
			$cmd .= "gzip -c $allFiles[0] > $allAss\nrm $allFiles[0]\n";
		}
	} else {
		if ($allAss !~ m/\.gz$/ && $allFiles[0] =~ m/\.gz$/){
			$allAss.=".gz" ;
		} elsif($allAss =~ m/\.gz$/ && $allFiles[0] !~ m/\.gz$/){ #zip output
			$catcmd = "gzip -c";
		}
		#my $cmd= "cat ".join(" ",@allFiles). " > $allAss\n";   #
		$cmd .= "$catcmd ".join(" ",@allFiles). " > $allAss\n";
		$cmd .= "rm -f ".join(" ",@allFiles) . "\n";
	}
	#die "$cmd";
	if ($curDB eq "mp3"){
		$cmd .= "mv $allAss ${allAss}geneAss.gz mp3\ntouch $allAss\n";
	} else {
		$cmd .= "$secCogBin -i $allAss -DB $shrtDB -singleSpecies 1  -bacNOG $opts{bacNOG} -KOfromNOG 0 -eggNOGmap 1 -calcGeneLengthNorm 0 -lenientCardAssignments 2 -mode 2 -CPU $ncore -percID 25 -LF $DBpath/$refDB.length -DButil $DBpath -tmp $tmpD -eggNOGmap 0 -minPercSbjCov 0.3 -minBitScore 60 -minAlignLen 60 -eval 1e-7\n";
	}
	$cmd .= "rm -f ".join(" ",@subFls)."\n" unless($opts{keepSplits});
	if ($interpDia && $doQsub){
		($jdep,$qCmd) = qsubSystem($qsubDir."colDIA$shrtDB.sh",$cmd,1,"40G","ColDIA",join(";",@jdeps),"",1,[],$QSBoptHR);
	} elsif ($interpDia) {
		systemW $cmd;
	}
	return $allAss,$jdep;
}



sub calc_modules{
	my ($inMat,$outD,$ModCompl,$EnzCompl) = @_;
	my $rareBin = getProgPaths("rare");
	my $modD = "/g/bork3/home/hildebra/DB/FUNCT/myModules/Feb16/";
	my @modDBfs = ("module_new.list","module_c.list","modg.list","module_s.list");
	my @modDBodir = ("modules","metaCyc","BSB","SEED");
	my @modDescr = ("mod.descr","modc.descr","modg.descr","mods.descr");
	my @modHiera = ("mod_hiera.txt","modc_hiera.txt","modg_hiera.txt","mods_hiera.txt");
	my @modShort = ("modKEEG","modMC","modGMM","modSEED");

	for (my $k=0;$k<@modDBfs;$k++){
		next if ($k==1); #skyp metacyc
		my $keggDB = $modD.$modDBfs[$k];
		my $outD2 = "$outD/$modDBodir[$k]";
		system "mkdir -p $outD2" unless (-d $outD2);
		my $outMat = "/".$outD2."KEGG";
		#die "$inMat\n";
		my $cmdMod = "$rareBin module -i $inMat -o $outMat -refMods $keggDB -description $modD/$modDescr[$k] -hiera  $modD/$modHiera[$k] -redundancy 5 -writeExtraModEstimates -moduleCompl $ModCompl -enzymeCompl $EnzCompl -collapseDblModules\n";
		die "Failed module calc:\n$cmdMod\n" if (system "$cmdMod");

	}
}
