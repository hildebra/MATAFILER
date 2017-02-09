#!/usr/bin/env perl
#main TAMOC routine
#ex
#./TAMOC.pl map2tar /g/scb/bork/hildebra/SNP/test/refCtg.fasta,/g/scb/bork/hildebra/SNP/test/refCtg.fasta test1,test2
#./TAMOC.pl map2tar /g/bork3/home/hildebra/results/TEC2/v5/TEC2.MM4.BEE.GF.rn.fa TEC2
#./TAMOC.pl map2tar /g/bork3/home/hildebra/results/prelimGenomes/TEC3/MM3.TEC3.scaffs.fna,/g/bork3/home/hildebra/results/prelimGenomes/TEC3/TEC3ref.fasta,/g/bork3/home/hildebra/results/prelimGenomes/TEC4/MM3.TEC4.scaffs.fna,/g/bork3/home/hildebra/results/prelimGenomes/TEC4/TEC4ref.fasta,/g/bork3/home/hildebra/results/prelimGenomes/TEC5/MM4.TEC5.scaffs.fna,/g/bork3/home/hildebra/results/prelimGenomes/TEC5/TEC5ref.fasta,/g/bork3/home/hildebra/results/prelimGenomes/TEC6/MM29.TEC6.scaffs.fna,/g/bork3/home/hildebra/results/prelimGenomes/TEC6/TEC6ref.fasta TEC3,TEC3r,TEC4,TEC4r,TEC5,TEC5r,TEC6,TEC6r
#/g/bork5/hildebra/results/TEC2/v5/T6/TEC6.ctgs.rn.fna
#./TAMOC.pl map2DB /g/bork3/home/hildebra/DB/freeze11/freeze11.genes.representatives.fa frz11; ./TAMOC.pl map2DB /g/bork3/home/hildebra/DB/GeneCats/Tara/Tara.fna Tara; ./TAMOC.pl map2DB /g/bork3/home/hildebra/DB/GeneCats/IGC/1000RefGeneCat.fna IGC
#./TAMOC.pl map2tar /g/bork3/home/hildebra/results/TEC2/v5/TEC2.MM4.BEE.GF.rn.fa,/g/bork3/home/hildebra/results/TEC2/v5/T6/TEC6.ctgs.rn.fna T2d,T6d

#The Metagenomic Assembly, Genomic Recovery and Assembly Independent Mapping Tool


use warnings;
use strict;
use File::Basename;
use Cwd 'abs_path';
use POSIX;
use Mods::GenoMetaAss qw(readMap qsubSystem emptyQsubOpt readFastHD prefix_find);
use Mods::IO_Tamoc_progs qw(getProgPaths jgi_depth_cmd inputFmtSpades createGapFillopt  buildMapperIdx);

#use Mods::TamocFunc qw(runDiamond);

sub seedUnzip2tmp; sub clean_tmp; sub cleanInput; #unzipping reads; removing these at later stages ; remove tmp dirs
sub readG2M; sub check_matesL;
sub mapReadsToRef;  sub bamDepth;
sub ReadsFromMapping;
sub spadesAssembly; #main Assembler
sub createPsAssLongReads; #pseudo assembler

sub scaffoldCtgs;  #scaffold assemblies / external contigs
sub GapFillCtgs;   #post scaffolding using reads from samples
sub filterSizeFasta;
sub sdmClean; sub check_sdm_loc; #qual filter reads
sub mergeReads; #merge reads via flash
sub SEEECER; 
sub run_prodigal_augustus; sub run_prodigal; #gene prediction
#sub randStr; 
sub contigStats;sub smplStats;
sub adaptSDMopt;#sub readMap; sub qsubSystem;
sub checkDrives;sub bam2cram;
sub mocat_reorder; sub postSubmQsub;
sub detectRibo; 
sub runDiamond; sub DiaPostProcess; sub IsDiaRunFinished;
sub nopareil; sub calcCoverage;
sub prepKraken;sub krakHSap; sub krakenTaxEst;
sub d2metaDist;
sub metphlanMapping;sub mergeMP2Table;

#bjobs | awk '$3=="CDDB" {print $1}' |xargs bkill
#bjobs | grep '_ABR' | cut -f1 -d' ' | xargs -t -i bkill {}
#bhosts | cut -f1 -d' ' | grep -v HOST_NAME | xargs -t -i ssh {} 'killall -u hildebra'
#hosts=`bhosts | grep ok | cut -d" " -f 1 | grep compute | tr "\\n" ","`; pdsh -w $hosts "rm -rf /tmp/hildebra"

#program configuration
my $Spades_Cores=48; my $Spades_Memory = 200; #in GB
my $Spades_HDspace = 100; #required space in GB
my $Spades_Kmers = "-k 27,33,55,71";
my $bwt_Cores = 12; my $map_DoConsensus = 1; my $doRmDup = 1; #mapping cores; ??? ; remove Dups (can be costly if many ref seqs present)
my $diaEVal = "0.0000001"; my $dia_Cores = 16; my $krakenCores = 9;
my $bwtMapMem = "3G";
#----------------- defaults ----------------- 
my $rawFileSrchStr1 = '.*1\.f[^\.]*q\.gz$';
my $rawFileSrchStr2 = '.*2\.f[^\.]*q\.gz$';
my $rawFileSrchStrSingl = "";
my $rawFileSrchStrXtra1= '.*1_sequence\.f[^\.]*q\.gz$';
my $rawFileSrchStrXtra2= '.*2_sequence\.f[^\.]*q\.gz$';
my %sdm_opt; #empty object that can be used to modify default sdm parameters

my $mateInsertLength =20000; #controls expected mate insert size , import for bowtie2 mappings
my %jmp=();
my $logDir = "";
my $sharedTmpDirP = "/scratch/bork/hildebra/";
$sharedTmpDirP = "/g/scb/bork/hildebra/tmp/" if (`hostname` !~ m/submaster/);
my $nodeTmpDirBase = "/tmp/hildebra/";


#die $sharedTmpDirP;

#programs --------------------------

my $smtBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";
my $metPhl2Bin = getProgPaths("metPhl2");#"/g/bork3/home/hildebra/bin/metaphlan2/metaphlan2.py";
my $metPhl2Merge = getProgPaths("metPhl2Merge");#"/g/bork3/home/hildebra/bin/metaphlan2/utils/merge_metaphlan_tables.py";
#my $spadesBin = "/g/bork5/hildebra/bin/SPAdes-3.7.0-dev-Linux/bin/spades.py";
my $spadesBin = getProgPaths("spades");#"/g/bork3/home/hildebra/bin/SPAdes-3.7.1-Linux/bin/spades.py";
my $usBin = getProgPaths("usearch");#"/g/bork5/hildebra/bin/usearch/usearch8.0.1421M_i86linux32_fh";
my $bwt2Bin = getProgPaths("bwt2");#"/g/bork5/hildebra/bin/bowtie2-2.2.9/bowtie2";
my $bwaBin = getProgPaths("bwa");#"/g/bork3/home/hildebra/bin/bwa-0.7.12/bwa";
my $d2metaBin = getProgPaths("d2meta");#"/g/bork3/home/hildebra/bin/d2Meta/d2Meta/d2Meta.out";
my $prodigalBin = getProgPaths("prodigal");#"/g/bork5/hildebra/bin/Prodigal-2.6.1/prodigal";
my $augustusBin = getProgPaths("augustus");#"/g/bork3/home/hildebra/bin/augustus-3.2.1/bin/augustus";
my $sambambaBin = getProgPaths("sambamba");#"/g/bork3/home/hildebra/bin/sambamba/sambamba_v0.5.9";
my $npBin = getProgPaths("nonpareil");#"/g/bork5/hildebra/bin/nonpareil/nonpareil";
my $krkBin = getProgPaths("kraken");#"/g/scb/bork/hildebra/DB/kraken/./kraken";
my $diaBin = getProgPaths("diamond");#"/g/bork5/hildebra/bin/diamond/./diamond";
my $nxtrimBin = getProgPaths("nxtrim");#"/g/bork3/home/hildebra/bin/NxTrim/./nxtrim";
my $besstBin = getProgPaths("BESST");#"/g/bork3/home/hildebra/bin/BESST/./runBESST";
#my $novosrtBin = getProgPaths("novosrt");#"/g/bork5/hildebra/bin/novocraft/novosort";
my $GFbin = getProgPaths("gapfiller");#"perl /g/bork5/hildebra/bin/GapFiller/GapFiller_n.pl";
my $flashBin = getProgPaths("flash");
my $pigzBin  = getProgPaths("pigz");
my $sdmBin = getProgPaths("sdm");#"/g/bork3/home/hildebra/dev/C++/sdm/./sdm";
my $rareBin = getProgPaths("rare");#"/g/bork3/home/hildebra/dev/C++/rare/rare";
my $readCov_Bin =getProgPaths("readCov");

my $baseDir = ""; my $baseOut = "";
my $mapF = ""; my $baseID = "empty_ass_metag";
#my $krkFiltBin = "/g/scb/bork/hildebra/DB/kraken/./kraken-filter";


#local TAMOC scripts --------------------------
my $cLSUSSUscript = getProgPaths("cLSUSSU_scr");#"perl /g/bork3/home/hildebra/dev/Perl/16Stools/catchLSUSSU.pl";
my $lotusLCA_cLSU = getProgPaths("lotusLCA_cLSU_scr");#"perl /g/bork3/home/hildebra/dev/Perl/16Stools/lotus_LCA_blast2.pl";
my $krakCnts1 = getProgPaths("krakCnts_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/krak_count_tax.pl";
my $genelengthScript = getProgPaths("genelength_scr");#= "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/geneLengthFasta.pl";
my $secCogBin = getProgPaths("secCogBin_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/parseBlastFunct.pl";
my $KrisABR = getProgPaths("KrisABR_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/ABRblastFilter.pl";
my $sepCtsScript = getProgPaths("sepCts_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/separateContigs.pl";
my $assStatScr = getProgPaths("assStat_scr");#"perl /g/bork3/home/hildebra/dev/Perl/assemblies/assemblathon_stats.pl";
my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/renameCtgs.pl";
my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl /g/bork3/home/hildebra/dev/Perl/assemblies/sizeFilterFas.pl";
my $sizSplitScr = getProgPaths("sizSplit_scr");#"perl /g/bork3/home/hildebra/dev/Perl/assemblies/splitFNAbyLength.pl";
my $splitKgdContig = getProgPaths("contigKgdSplit_scr");
my $decoyDBscr = getProgPaths("decoyDB_scr"); #
my $bamHdFilt_scr = getProgPaths("bamHdFilt_scr");
my $mrgDiScr = getProgPaths("mrgDia_scr");


my $mergeTblScript = "/g/bork3/home/hildebra/bin/metaphlan2/utils/merge_metaphlan_tables.py";
my $mergeMiTagScript = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/miTagTaxTable.pl";



#databases --------------------------
my $globalKraTaxkDB = "";
my %globalRiboDependence =(DBcp => "");
my %globalDiamondDependence = (CZy=>"",MOH=>"",NOG=>"",ABR=>"",ABRc=>"",KGB=>"",KGE=>"",ACL=>"",KGM=>"" );




#more specific control: unfiniRew=rewrite unfinished sample dir; $redoCS = redo ContigStats completely; 
#removeInputAgain=remove unzipped files from scratch, after sdm; remove_reads_tmpDir = leave cleaned reads on scratch after everything finishes
my $unfiniRew=0; my $redoCS=0; my $removeInputAgain=1; my $remove_reads_tmpDir=0;
my $readsRpairs=1; #are reads given in pairs?
my $splitFastaInput = 0; #assembly as input..
my $importMocat = 0; my  $mocatFiltPath = "reads.screened.screened.adapter.on.hg19.solexaqa/"; 
my $alwaysDoStats = 1;
my $humanFilter = 1; #use kraken to filter out human reads
my $DoMetaPhlan = 0; my $metaPhl2FailCnts=0; #non-asembly based tax + functional assignments
my $DoRibofind = 0; my $doRiboAssembl = 0; my $RedoRiboFind = 0; my $riboFindFailCnts=0; #ITS/SSU/LSU detection
my $RedoRiboAssign = 0;
my $DoBinning = 1;
my $doReadMerge = 0;
my $DoAssembly = 1;  my $SpadesAlwaysHDDnode = 1;my $spadesBayHam = 0; my $useSDM = 2;my $spadesMisMatCor = 0; my $redoAssembly =0 ;
my $doBam2Cram= 1;
my $DoNP=0; #non-pareil
my $doDateFileCheck = 0; #very specific option for Moh's reads that were of different dates..
my $DoDiamond = 0; my $rewriteDiamond =0; my $redoDiamondParse = 0; #redoes matching of reads; redoes interpretation
my $MapperProg = 1;#1=bowtie2, 2=bwa
my $DO_EUK_GENE_PRED = 0;
my $calcD2s = 0; 
my $DoKraken = 0; my $RedoKraken = 0; my $KrakTaxFailCnts=0;
my $pseudoAssembly = 0; #in case no assembly is possible (soil single reads), just filter for reads X long 
my $DoFreeGlbTmp = 0; my $defaultReadLength = 100;
my $maxReqDiaDB = 6; #max number of databases supported by TAMOC
my $reqDiaDB = "ABR";#,NOG,MOH,ABR,ABRc,ACL,KGM";#,ACL,KGM,ABRc,CZy";#"NOG,CZy"; #"NOG,MOH,CZy,ABR,ABRc,ACL,KGM"   #old KGE,KGB




#control broad flow
my $from=0; my $to=500;

my $doSubmit=1; 
my $rewrite=0;my $rewrite2ndMap = 0;


if ( 1 ){ #real data 
	if ( 0 ){ #MM samples
		$importMocat=0;	
		die "Sure you want to rewrite?\n" if ($rewrite);
		#$jmp{$_} = 1 for qw(11 );
	#$baseID = "GNMass_mocat";
		if (1){ # compound assembly
			$Spades_Cores=70;$Spades_Memory = 170;$Spades_HDspace=100;
			$baseID = "GNMass3";
			#MH samples: my $from=137; my $to=150;
			#$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MM_at_v3_T2subset.txt";
			$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MM_at_v5_T2subset.txt";
			#$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/T2_MetaHit.txt"; #no read access...
			#$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/T2_HMP.txt";
			#$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MM_at_v3.txt";
			#testing assembly options
			my $SpadesTest=0;
			if ($SpadesTest==1){
				$baseID.="bHmm__meta";$spadesBayHam = 1; $useSDM = 0;
			} elsif($SpadesTest==2){
				$baseID.="_meta";$spadesBayHam = 0; $useSDM = 1;
			} elsif($SpadesTest==3){
				$baseID.="_bHmm";$spadesBayHam = 1; $useSDM = 0;
			} elsif($SpadesTest==4){
				$baseID.="_meta_bHmSdm";$spadesBayHam = 1; $useSDM = 1;
			}
		} elsif (0){ #run for Ana
			$Spades_Cores=16;$Spades_Memory = 80;
			$baseID = "Ana";$importMocat = 1;
			$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MM_ana.txt";
		} else {
			$Spades_Cores=18;$Spades_Memory = 140;# single assembly
			$baseID = "GNMass2_singl"; #GNMass2_singl GNMass3_singl
			$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MM_at_singl.txt";
			#$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MM_at_v3_T2subset.txt";
			
		}
		#
		$baseOut = "/g/scb/bork/hildebra/SNP/$baseID/";
		$baseDir = "/g/bork5/kultima/MM/";
#		$Spades_Kmers = "-k 27,33,55,71,99";
		$Spades_Kmers = "-k 25,43,67,87,127";
	} elsif (0) { #ColPark 
		$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/PD_cologne.txt";
		$importMocat=1; $useSDM =0; $DoBinning=0; $globalKraTaxkDB = "virusDB";
		$DoDiamond = 0; $redoDiamondParse=0; $rewriteDiamond=0;
		$DoRibofind =1; $doRiboAssembl = 0; $RedoRiboFind = 0;  $RedoRiboAssign=0;		$DoAssembly=0;  $humanFilter =0;
		$DoMetaPhlan = 0; $DoKraken = 1; $globalKraTaxkDB = "virusDB";#"soilOrgs";#"minikraken_2015/";
		$RedoKraken=0;
	} elsif (0) { #Louis dog samples
		$baseID = "Dogs";
		$baseOut = "/g/scb/bork/hildebra/Tamoc/$baseID/";
		$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/Dogs.txt";
		$baseDir = "/g/bork5/mocat/DOGS/";
	} elsif (0) { #HMP samples
		#$jmp{$_} = 1 for qw(201 188 209 248);
		$baseID = "HMPtrial";$Spades_Memory = 100;
		$baseOut = "/g/scb/bork/hildebra/Tamoc/$baseID/";
		$Spades_Kmers = "-k 21,33,67,77,89";
		$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/HMP_ana.txt";
		$baseDir = "/g/bork5/mocat/HMP_stool//";
	} elsif (0) { #FruitBodies Moh
		$baseID = "FruitBodies";
		$baseOut = "/g/scb/bork/hildebra/Tamoc/$baseID/";
		$mapF = "/g/scb/bork/hildebra/data2/FruitBodies/Fruits.txt";
		$baseDir = "/g/scb/bork/hildebra/data2/FruitBodies/";		
		$DoNP=0;$DoRibofind = 1;$DoDiamond = 1; 
		$rawFileSrchStr1 = '.*R1_00\d\.f[^\.]*q\.gz$';
		$rawFileSrchStr2 = '.*R2_00\d\.f[^\.]*q\.gz$';
	} elsif (0) { #DeepSeq 40
		$baseID = "DeSe40";
		$baseOut = "/g/scb/bork/hildebra/Tamoc/$baseID/";$mapF = "/g/scb/bork/hildebra/data2/Soil_finland/DeepSeq40/map.txt";
		$baseDir = "/g/scb/bork/hildebra/data2/Soil_finland/DeepSeq40/";$DoNP=0;$DoRibofind = 0;$DoDiamond = 0; $DoAssembly=1;
		$rawFileSrchStr1 = '.*1\.f[^\.]*q\.gz$';$rawFileSrchStr2 = '.*2\.f[^\.]*q\.gz$';
	} elsif (0){ #TARA
		#$from=0;$to = 300;#
		#$to = 60; #restrict myself to first 60 samples for now..
		$mapF = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/maps/Tara_bac.map";
		#$mocatFiltPath = "reads.screened.screened.adapter.on.hg19.solexaqa/";
		$importMocat = 1; $DoBinning=0;
		#$jmp{$_} = 1 for qw(102 118);
		my @dvals = ("1e-8","1e-9","1e-10","1e-11","1e-12","1e-13","1e-14","1e-15");$diaEVal=join(",",@dvals);
		#$rawFileSrchStr1 = '.*_1\.f[^\.]*q\.gz$'; $rawFileSrchStr2 = '.*_2\.f[^\.]*q\.gz$';
		$DoDiamond = 1; $redoDiamondParse=0; $rewriteDiamond=0;
		$DoRibofind = 0; $doRiboAssembl = 0; $RedoRiboFind = 0;  $RedoRiboAssign=0;
		$DoAssembly=0;  $humanFilter =0;
		$DoMetaPhlan = 0; $DoKraken = 0; $globalKraTaxkDB = "soilOrgs";#"minikraken_2015/";
		$RedoKraken=0;
	} elsif ( 0 ) { #Anna single cell seq
		$rawFileSrchStr1 = '.*1_sequence\.txt\.gz$';$rawFileSrchStr2 = '.*2_sequence\.txt\.gz$';
		$DoNP=0;$humanFilter=0; $calcD2s=0; $DoDiamond = 0;  $redoDiamondParse = 0;  $rewriteDiamond=0;
		$DoRibofind = 0; $doRiboAssembl = 0; $RedoRiboFind = 0; $RedoRiboAssign=0;
		$DoAssembly=1; $DoBinning=1; $Spades_Kmers = "-k 27,33,55,71,101,127";
		$mapF = "/g/bork3/home/hildebra/data/AnnaPry/singleCell.map";
		#$mapF = "/g/bork3/home/hildebra/data/AnnaPry/singleCell_comb.map";
	} elsif ( 0 ) { #Yongkyu fungi samples
		$rawFileSrchStr1 = '.*1_sequence\.f[^\.]*q\.gz$';
		$rawFileSrchStr2 = '.*2_sequence\.f[^\.]*q\.gz$';
		$DoNP=0;$humanFilter=0; $calcD2s=0; $DoDiamond = 1;  $redoDiamondParse = 0;  $rewriteDiamond=0;
		$DoAssembly=1;$DoRibofind = 1; $doRiboAssembl = 0; $RedoRiboFind = 0; $RedoRiboAssign=0; $DoBinning=0;
		$DoMetaPhlan = 0; $DoKraken = 0; $globalKraTaxkDB = "minikraken_2015";#= "minikraken_2015/";
		my @dvals = ("1e-8","1e-9","1e-10","1e-11","1e-12","1e-13","1e-14","1e-15");$diaEVal=join(",",@dvals);
		$mapF = "/g/bork3/home/hildebra/data/Kieran/Yongkyu/kefir1.txt"; 
	} elsif (0) { #IHMS samples
		#$from=305;
		$rawFileSrchStr1 = '.*1\.f[^\.]*q\.gz$';$rawFileSrchStr2 = '.*2\.f[^\.]*q\.gz$'; $importMocat=1;
		$DoNP=0;$humanFilter=0;$DoAssembly=1; $DoBinning=1; $Spades_Cores=48;$Spades_Memory = 480; $dia_Cores = 14;$DoBinning=0; $redoAssembly=0;
		$Spades_Kmers = "-k 21,33,67,111,127";
		$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/IHMS.txt"; 
	} elsif (0) { #MetaHit samples
		$DoNP=0;$humanFilter=0; $calcD2s=0; $DoDiamond = 1;  $redoDiamondParse = 0;  $rewriteDiamond=0;
		$DoRibofind = 0; $doRiboAssembl = 0; $RedoRiboFind = 0; $RedoRiboAssign=0;$DoAssembly=0; $DoBinning=0;
		$DoMetaPhlan = 0; $DoKraken = 0; $globalKraTaxkDB = "virusDB";#= "minikraken_2015/";
		my @dvals = ("1e-8","1e-9","1e-10","1e-11","1e-12","1e-13","1e-14","1e-15");$diaEVal=join(",",@dvals);
		$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MetaHit.txt"; 
	} elsif ( 0 ) { #NIH & HMP skin samples
		$from = 0;
		$DoNP=0;$humanFilter=1; $calcD2s=0; $DoDiamond = 1;  $redoDiamondParse = 0;  $rewriteDiamond=0;
		$DoRibofind = 0; $doRiboAssembl = 0; $RedoRiboFind = 0; $RedoRiboAssign=0;$DoAssembly=0; $DoBinning=0;
		$DoMetaPhlan = 0; $DoKraken = 0; $globalKraTaxkDB = "virusDB";#= "minikraken_2015/";
		my @dvals = ("1e-8","1e-9","1e-10","1e-11","1e-12","1e-13","1e-14","1e-15");$diaEVal=join(",",@dvals);
		$mapF = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/maps/HMPNIH_skin.map";
	} else { #soil samples
		#$jmp{$_} = 1 for qw(201 188 209 248);
		#$baseID = "FinSoil";
		#die "TODO: restart NOG dia from 191\n"; 
		#$from=11; 
		#$to=12;
		
		$DO_EUK_GENE_PRED=0;
		$doRmDup=0;
		$removeInputAgain=0; $doReadMerge=1; $remove_reads_tmpDir=0;
		$bwtMapMem = "5G"; #for mapping to refgenomes
		$DoNP=0;$humanFilter=0; $calcD2s=0; 
		$DoDiamond = 1;  $redoDiamondParse = 0;  $rewriteDiamond=0;
		$DoRibofind = 0; $doRiboAssembl = 0; $RedoRiboFind = 0; $RedoRiboAssign=0;
		$defaultReadLength = 250;
		$DoMetaPhlan = 0; $DoKraken = 0; $globalKraTaxkDB = "virusDB";#= "minikraken_2015/";
		$DoAssembly=1; $Spades_Cores=48;$Spades_Memory = 800; $dia_Cores = 14;$DoBinning=0; $redoAssembly=0;
		$Spades_Kmers = "-k 21,33,67,111,127";
		$rawFileSrchStr1 = '.*R1_00\d\.f[^\.]*q\.gz$|.*_1\.f[^\.]*q\.gz$';
		$rawFileSrchStr2 = '.*R2_00\d\.f[^\.]*q\.gz$|.*_2\.f[^\.]*q\.gz$';
		$DoBinning=0;
		my @dvals = ("1e-8","1e-9","1e-10","1e-11","1e-12","1e-13","1e-14","1e-15");$diaEVal=join(",",@dvals);
		$baseOut = "/g/scb/bork/hildebra/Tamoc/$baseID/";
		$mapF = "/g/scb/bork/hildebra/data2/Soil_finland/soil_map.txt"; 
		#$mapF = "/g/scb/bork/hildebra/data2/refData/Soil/PNAS_refs/other_soil.map";$DoAssembly=0;$pseudoAssembly =1;$rawFileSrchStrSingl='.*\.fq\.gz$'; $readsRpairs=0;$sdm_opt{maxSeqLength}=10000000;$sdm_opt{minSeqLength}=90; $splitFastaInput = 1000;
		#$mapF = "/g/scb/bork/hildebra/data2/refData/Soil/howe2014/iowa.map";$Spades_Kmers = "-k 21,33,67,77";$doReadMerge=0; $from=0;$to=10;$DoDiamond=1;$defaultReadLength=100;

		#$mapF = "/g/scb/bork/hildebra/data2/Soil_finland/soil_map_d.txt"; $doDateFileCheck = 1; $reqDiaDB = "CZy"; $DoRibofind = 0;$DoMetaPhlan = 0; $DoKraken = 0;
	}
} else {
	$from=0; $to=50;$rewrite =0;
	$doBam2Cram =0; $DoAssembly=1;
	$Spades_Cores=12;$Spades_Memory = 24;
	$Spades_HDspace=100;
	$bwt_Cores = 5;$humanFilter=0;
	$doReadMerge=0;
	if (1){
		$baseID = "SimuA";  $DoRibofind=0; $RedoRiboFind = 0;$RedoRiboAssign=0; $DO_EUK_GENE_PRED=1;
		$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/simus2.txt";#_singl
		$baseOut = "/g/scb/bork/hildebra/SNP/$baseID/"; 
		$RedoKraken=0; $DoKraken = 0; $globalKraTaxkDB = "soilOrgs";#$globalKraTaxkDB = "minikraken_2015/";
		$baseDir = "/g/bork5/hildebra/data/Simulations/";#"/g/bork5/kultima/MM/"
	} else { #fungi & 2 bacs, 8 times repeated
		$DO_EUK_GENE_PRED = 1;
		
		#$reqDiaDB = "KGE,KGB";
		$DoBinning=0;$DoAssembly = 1; $DoMetaPhlan=0;
		$DoDiamond = 0; $redoDiamondParse = 0; $rewriteDiamond=0;
		$DoRibofind = 0; $doRiboAssembl = 0;$RedoRiboFind = 0;$RedoRiboAssign=0;
		$DoKraken = 0;
		my @dvals = ("1e-7","1e-8","1e-9","1e-10","1e-11","1e-12","1e-13","1e-14","1e-15","1e-16","1e-17","1e-18","1e-19","1e-20");
		#my @dvals = ("0.0000001","0.000001","0.00001","0.0001","0.001","0.00000001","0.000000001");
		$diaEVal=join(",",@dvals);
		$globalKraTaxkDB = "minikraken_2015/";
		$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/simu3.txt";#200, 250 bp
		#$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/simu4.txt";#200, 100 bp
		#$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/simus_Fungi2.txt";#_singl
		$baseOut = "/g/scb/bork/hildebra/SNP/$baseID/";
		$baseDir = "/g/bork5/hildebra/data/Simulations/";#"/g/bork5/kultima/MM/"
	}
}

#die $mapF;

#die $baseOut;
#queing capability
my $JNUM=0;
#my $LocationCheckStrg=""; #command that is put in front of every qsub, to check if drives are connected, sub checkDrives
my $QSBoptHR = emptyQsubOpt($doSubmit,"");
my %QSBopt = %{$QSBoptHR};
my @Spades_Hosts = (); my @General_Hosts = ();
#figure out if only certain node subset has enough HDD space
if (0 && $Spades_HDspace > 100){
	my $locHosts = `bhosts | grep ok | cut -d" " -f 1 | grep compute | tr "\\n" ","`;
	my $tmpStr = `pdsh -w $locHosts -u 4 "df -l" `;
	my $srchTerm = '/$';
	if (`hostname` =~ m/submaster$/){
		$srchTerm = '/tmp';
	}
	#print "psdh done\n";
	#die $tmpStr;
	#54325072 = 52G
	foreach my $l (split(/\n/,$tmpStr)){
		next if ($l !~ m/compute\S+:/);
		next unless ($l =~ m/$srchTerm/);
		#print $l."\n";
		my @hosts = split(/:/,$l);
		my @spl = split(/\s+/,$hosts[1]);
		#die $spl[1]." XX ".$spl[2]." XX ".$spl[3]." XX ".$spl[4]." \n ";
		#print $spl[4] / 1024 / 1024 . "\n";
		if (($spl[4] / 1024 / 1024) > $Spades_HDspace){
			push(@Spades_Hosts,$hosts[0]);
		}
		if (($spl[4] / 1024 / 1024) > 40){
			push(@General_Hosts,$hosts[0]);
		}
	}
	print "Found ".scalar @Spades_Hosts." host machines with > $Spades_HDspace G space\n";
	print "Found ".scalar @General_Hosts." host machines with > 40 G space\n";
	if (scalar @Spades_Hosts ==0 ){die "Not enough hosts for spades temp space found\n";}
	sleep (2);
}



#do I need to redo d2s?
if ($calcD2s) {$calcD2s = !-e "$baseOut/d2StarComp/d2meta.stone";}

#die();
#----------- map all reads to a specific reference ---------
my $prevDeps= "";my $mapModeActive=0; my $mapModeCovDo=1;
my @bwt2outD =(); my @DBbtRefX = (); my @DBbtRefGFF=(); #my $bwt2Name=""; 

#here the map and some base parameters (base ID, in path, out path) can be (re)set
my ($hr,$hr2) = readMap($mapF);
my %AsGrps = %{$hr2};
my %map = %{$hr};
#$baseDir = $map{inDir} if (exists($map{inDir} ));
$baseOut = $map{outDir} if (exists($map{outDir} ));
$baseID = $map{baseID} if (exists($map{baseID} ));
my $runTmpDirGlobal = "$sharedTmpDirP/$baseID/";
my $runTmpDBDirGlobal = "$runTmpDirGlobal/DB/";
system "mkdir -p $runTmpDBDirGlobal" unless (-d $runTmpDBDirGlobal);
my $globaldDiaDBdir = $runTmpDBDirGlobal."DiamDB/";
system "mkdir -p $globaldDiaDBdir" unless (-d $globaldDiaDBdir);


#here are scaffolding external contigs parameters
my $scaffTarExternal = "";my $scaffTarExternalName = ""; my @scaffTarExternalOLib1; my @scaffTarExternalOLib2;
my $scaffTarExtLibTar = ""; my $bwt2ndMapDep = ""; my @bwt2ndMapNmds;
my %make2ndMapDecoy;
if (@ARGV > 2 && ($ARGV[0] eq "map2tar" || $ARGV[0] eq "map2DB") ){
#in this case primary focus is on mapping and not on assemblies
	if ($ARGV[0] eq "map2DB"){$mapModeCovDo=0;}
	my $refDBall = $ARGV[1];
	my $bwt2NameAll = $ARGV[2];
	my @refDB = split(/,/,$refDBall);
	my @bwt2Name = split(/,/,$bwt2NameAll);
	$MapperProg=1;  $bwtMapMem = "6G";$bwt_Cores=12;
	
	#decoy mapping setup (only required in map2tar
	$make2ndMapDecoy{Lib} = "";
	#die "decoy mapping not ready for multi fastas\n" if (@refDB > 1);
	
	print "\n=======================\nmap to $refDBall\n=======================\n\n";
	$mapModeActive =1;
	for (my $i=0;$i<@refDB; $i++){
		#$refDB[$i] =~ m/(.*\/)[^\/]+/;
		#my $refDir = $1;
		my $bwt2outDl = "$baseOut/GlbMap/$bwt2Name[$i]/";
		push @bwt2ndMapNmds , $bwt2Name[$i];
		push(@bwt2outD,$bwt2outDl);
		#die $bwt2Name."\n";
		#my $mapF = "$baseOut/LOGandSUB/inmap.txt";
		if ($rewrite || $rewrite2ndMap){
			print "Deleting previous mapping results..\n";
			$refDB[$i] =~ m/.*\/([^\/]+)$/;
			system("rm -r -f $bwt2outDl $runTmpDBDirGlobal/$1*");#;mkdir -p $bwt2outDl
		}
		$refDB[$i] =~ m/.*\/([^\/]+)$/;
		#print "\n$refDB[$i]\n";
		#die "$bwt2outDl/$1\n";
		system "mkdir -p $bwt2outDl/LOGandSUB" unless (-d "$bwt2outDl/LOGandSUB");
		system "cp $refDB[$i] $bwt2outDl" if ($mapModeCovDo && !-e "$bwt2outDl/$1");
		my $bwtDBcore = 20; my $largeDB = 1;
		my ($cmd,$DBbtRef) = buildMapperIdx($refDB[$i],$bwtDBcore,$largeDB,$MapperProg) ;
		$DBbtRef =~ s/\.bw2$//;
		$DBbtRef =~ m/.*\/([^\/]+)$/;
		$DBbtRef = "$runTmpDBDirGlobal/$1";#set up to scratch dir to map onto
		$cmd.= "\ncp $refDB[$i]* $runTmpDBDirGlobal\n" unless (-e "$DBbtRef.pak" || ($MapperProg==1 && -e "$DBbtRef"));# if (!$mapModeCovDo && !-e "$bwt2outDl/$1");
		#print $cmd."\n";
		#die $DBbtRef."\n$runTmpDBDirGlobal/\n";
		unless (-e "$DBbtRef.bw2.0.sa" || -e  "$DBbtRef.bw2.rev.1.bt2l"){
			#system $cmd 
			($bwt2ndMapDep,$cmd) = qsubSystem($bwt2outDl."/LOGandSUB/builBwtIdx$i.sh",$cmd,$bwtDBcore,int(400/$bwtDBcore)."G","BWI".$i,"","",1,[],\%QSBopt) ;
		}
		
		#$DBbtRefX = $DBbtRef;
		push(@DBbtRefX,$DBbtRef);
		if($mapModeCovDo){
			my $gDir = $bwt2outDl."";
			my $nativeGFF = $refDB[$i];$nativeGFF =~ s/\.[^\.]+/\.gff/;
			my $gffF = "genes.$bwt2Name[$i].gff";
			if (-e $nativeGFF){
				system "cp $nativeGFF $gDir/$gffF";
			} else {
				system "mkdir -p $gDir";
				$logDir = $gDir;
				my $dEGP = $DO_EUK_GENE_PRED; $DO_EUK_GENE_PRED = 0;
				my $tmpDep1 = run_prodigal_augustus($refDB[$i],$gDir,"",$gDir,"iGP$i","");
				$DO_EUK_GENE_PRED = $dEGP;
				my ($tmpDep,$tmpCmd) = qsubSystem( $bwt2outDl."/LOGandSUB/cpGenes.sh",  "cp $gDir/genes.gff $nativeGFF; mv $gDir/genes.gff $gDir/$gffF",
				1,"1G","genecop".$i,$tmpDep1,"",1,[],\%QSBopt);
				$prevDeps .= ";".$tmpDep;
			}
			push(@DBbtRefGFF,$gDir."/$gffF");#"genePred/genes.gff"
			# die "C";
		}
		$logDir="";
		
		my $aref = readFastHD($refDB[$i]);
		push(@{$make2ndMapDecoy{regions}}, join(" ",@{$aref}) );
		#die "@{$aref}\n".prefix_find($aref)."\n";
		push(@{$make2ndMapDecoy{region_lcs}},  prefix_find($aref)  );
	}
	#die @bwt2outD."  @bwt2outD\n";
	#die $prevDeps;
} elsif (@ARGV > 2 && $ARGV[0] eq "scaffold"){
	$scaffTarExternal = $ARGV[1];
	if (!-f $scaffTarExternal){
		die "Could not find scaffold file:\n$scaffTarExternal\n";
	}
	$scaffTarExternalName = $ARGV[2];
	if (@ARGV>3){
		$scaffTarExtLibTar = $ARGV[3];
	}
}

#die "@DBbtRefGFF\n";
#die "@{$make2ndMapDecoy{regions}}\n";
#die "$bwt2ndMapDep\n";
#my $baseOut = "/g/bork1/hildebra/SNP/GNMass/";
my $dir2rd=""; 
my $sequencer = "hiSeq";#plattform the algos have to deal with
my $DBpath=""; my $assDir=""; 
my $continue_JNUM = 0; #debug, set to 0 for full run
my $prevAssembly = ""; my $shortAssembly = "";#files with full length and short length assemblies
my $AssJob = "";#job name for assembly, restriction for mapping to this assembly
my $curDir = "";my $curOutDir = ""; 
my $mmpuOutTab = "";

my $dir_MP2 = $baseOut."pseudoGC/Phylo/MP2/"; #metaphlan 2 dir
my $dir_RibFind = $baseOut."pseudoGC/Phylo/RiboFind/"; #ribofinder dir
my $dir_KrakFind = $baseOut."pseudoGC/Phylo/KrakenTax/$globalKraTaxkDB/"; #kraken dir


system("mkdir -p $baseOut");
my $globalLogDir = $baseOut."LOGandSUB/";
system("mkdir -p $globalLogDir/sdm");
open $QSBopt{LOG},">",$globalLogDir."qsub.log";# unless ($doSubmit == 0);
print $globalLogDir."qsub.log\n";
my $collectFinished = $baseOut."runFinished.log\n";
my $globalNPD = $baseOut."NonPareil/";
system "mkdir -p $globalNPD";
system "cp $mapF $globalLogDir/inmap.txt";

my $present = 0; my $totalChecked=0;
my @samples = @{$map{smpl_order}}; my @allSmplNames;
my @allFilter1; my @allFilter2;
if ($to > @samples){
	print "Reset range of samples to ". @samples."\n";
	$to = @samples;
	#die();
}
my $statStr = ""; my $statStr5 = "";
my %sampleSDMs; 


my $baseSDMopt = getProgPaths("baseSDMopt_rel"); #"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/data/sdm_opt_inifilter_relaxed.txt";
if ($useSDM ==2 ){$baseSDMopt = getProgPaths("baseSDMopt");}#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/data/sdm_opt_inifilter.txt";}
my $baseSDMoptMiSeq = getProgPaths("baseSDMoptMiSeq_rel");#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/data/sdm_opt_miSeq.txt";	
if ($useSDM ==2 ){$baseSDMoptMiSeq = getProgPaths("baseSDMoptMiSeq");}#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/data/sdm_opt_miSeq_relaxed.txt";	}

my @unzipjobs; my $curUnzipDep = "";
my $sdmjNamesAll = "";
my $waitTime = 0;

my $emptCmd = "sleep 333";
#qsubSystem($globalLogDir."emptyrun.sh","sleep 333",1,"1G",0,"empty","","",1);


#set up kraken human filter
my $krakDeps = ""; my $krakenDBDirGlobal = $runTmpDirGlobal;
if ($DoKraken && $globalKraTaxkDB eq ""){die "Kraken tax specified, but no DB specified\n";}
if ($humanFilter || ($DoKraken) || $DO_EUK_GENE_PRED){
	$krakDeps = prepKraken();
}






#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#die $to."\n";
for ($JNUM=$from; $JNUM<$to;$JNUM++){
	#die $samples[$JNUM]."\n";
	next if (exists($jmp{$JNUM}));
	$dir2rd = $map{$samples[$JNUM]}{dir};
	
	#print "$samples[$JNUM]  $map{$samples[$JNUM]}{dir}  $map{$samples[$JNUM]}{rddir}  $map{$samples[$JNUM]}{wrdir}\n";	next;
	my $SmplName = $map{$samples[$JNUM]}{SmplID};
	push (@allSmplNames,$SmplName);
	if ($dir2rd eq "" ){#very specific read dir..
		if ($map{$samples[$JNUM]}{SupportReads} ne ""){
			$curDir = "";	$curOutDir = "$baseOut$SmplName/";	
		} else {
			die "Can;t find valid path for $SmplName\n";
		}
		$dir2rd = $SmplName;
	} else {
		$curDir = $map{$samples[$JNUM]}{rddir};	$curOutDir = $map{$samples[$JNUM]}{wrdir};	
	}

	#die "$curOutDir\n$baseOut\n";
	
	my $samplReadLength = $defaultReadLength; #some default value
	if (exists $map{$samples[$JNUM]}{readLength}){
		$samplReadLength = $map{$samples[$JNUM]}{readLength};
	}
	my $curSDMopt = $baseSDMopt;
	#my $iqualOff = 33; #62 for 1st illu

	if (exists($map{$samples[$JNUM]}{SeqTech})){
		if ($map{$samples[$JNUM]}{SeqTech} eq "GAII_solexa" || $map{$samples[$JNUM]}{SeqTech} eq "GAII"){
			#$iqualOff = 59; #really that old??
		} elsif ($map{$samples[$JNUM]}{SeqTech} eq "miSeq"){ 
			$baseSDMopt = $baseSDMoptMiSeq; 
		}
	} 
	if (!exists($sampleSDMs{$samplReadLength})){
		$sampleSDMs{$samplReadLength} = adaptSDMopt($baseSDMopt,$globalLogDir,$samplReadLength);
	}
	my $cAssGrp = $JNUM;
	#print $map{$samples[$JNUM]}{AssGroup}."\n";
	if ($map{$samples[$JNUM]}{AssGroup} ne "-1"){ $cAssGrp = $map{$samples[$JNUM]}{AssGroup};}
	my $cMapGrp = $map{$samples[$JNUM]}{MapGroup};
	
	$curSDMopt = $sampleSDMs{$samplReadLength};
	$totalChecked++;
	
	#die $SmplName."\n";
	print "\n======= $dir2rd - $JNUM - $SmplName =======\n";
	#my $dir2rd = "alien2-11-0-0";
	my $curTmpDir = "$runTmpDirGlobal$dir2rd/";
	my $nodeSpTmpD = "$nodeTmpDirBase/$SmplName";

	my $GlbTmpPath = "$curTmpDir/tmp/";
	my @checkLocs = ($GlbTmpPath,$curTmpDir,$curDir);
	if ($curDir eq ""){
		@checkLocs = ($GlbTmpPath,$curTmpDir);
	} 
	#$GlbTmpPath = "$curOutDir/tmp/";
	

	my $mapOut = "$GlbTmpPath/mapping/";
	$DBpath="$curOutDir/readDB/";
	my $finalCommAssDir = "$curOutDir/assemblies/metag/";
	if (-d "$finalCommAssDir/mismatch_corrector" || -d "$finalCommAssDir/tmp"){
		print "Removing temporary Spades Dirs..\n";
		system "rm -rf $finalCommAssDir/mismatch_corrector/ $finalCommAssDir/tmp";
	}
	
	$AsGrps{$cAssGrp}{AssemblSmplDirs} .= $curOutDir."\n";
	my $AssemblyGo=0; my $MappingGo=0; #controls if assemblies / mappings are done in respective groups
	if ($cAssGrp eq "-1"){
		die "Should not be here";
		$assDir="$GlbTmpPath/assemblies/"; #"$curOutDir/assemblies"
		$AsGrps{$cAssGrp}{AssemblJobName} = "_ASS$JNUM";

		$AssemblyGo = 1;
	} else {
		$assDir="$runTmpDirGlobal/AssmblGrp_$cAssGrp/";
		die $baseID."\n" if ($cAssGrp eq "");
		$finalCommAssDir = "$baseOut/AssmblGrp_$cAssGrp/metag/" if ( !exists($AsGrps{$cAssGrp}{CntAimAss}) || $AsGrps{$cAssGrp}{CntAimAss}>1);
		#assign job name (dependency) only ONCE
		if ( !exists($AsGrps{$cAssGrp}{CntAss}) || $AsGrps{$cAssGrp}{CntAss} == 0){
			$AsGrps{$cAssGrp}{AssemblJobName} = "_XXASpl$cAssGrp"."XX_" ;
			$AsGrps{$cAssGrp}{CSfinJobName} = "_XXCSpl$cAssGrp"."XX_" ;
		}
		$AsGrps{$cAssGrp}{CntAss} ++;
		print $AsGrps{$cAssGrp}{CntAss} .":".$AsGrps{$cAssGrp}{CntAimAss};
		if ($AsGrps{$cAssGrp}{CntAss}  >= $AsGrps{$cAssGrp}{CntAimAss} ){
			#print "running assembluy";
			$AssemblyGo = 1;
			if ($AsGrps{$cAssGrp}{CntAimAss}==1){
				$assDir="$GlbTmpPath/assemblies/";
			}
		}
		
		#mapping groups?
		$AsGrps{$cMapGrp}{CntMap} ++;
		print "  ".$AsGrps{$cMapGrp}{CntMap} .":".$AsGrps{$cMapGrp}{CntAimMap}."\n";
		if (!exists($AsGrps{$cMapGrp}{CntMap})){ die "Can;t find CntMap for $cMapGrp";}
		if ($AsGrps{$cMapGrp}{CntMap}  >= $AsGrps{$cMapGrp}{CntAimMap} ){
			#print "running assembluy";
			$MappingGo = 1;
		}
	}
	if ($DoAssembly ==0 ){$AssemblyGo=0;}
	$logDir = "$curOutDir/LOGandSUB/";
	my $statsDone=0;
	#maybe just do stats on whatever is there??
	if ($rewrite){
		print "Deleting previous results..\n";
		system("rm -f -r $curOutDir $GlbTmpPath $collectFinished ");
		system ("rm -r -f $assDir $finalCommAssDir");
	} 
	if ($alwaysDoStats){
		my ($statsHD,$curStats,$statsHD5,$curStats5) = smplStats($curOutDir,$assDir);
		$statsDone=1;
		if ($statStr eq ""){
			$statStr.="SMPLID\tDIR\t".$statsHD."\n".$samples[$JNUM]."\t$dir2rd\t".$curStats."\n";
			$statStr5.="SMPLID\tDIR\t".$statsHD5."\n".$samples[$JNUM]."\t$dir2rd\t".$curStats5."\n";
		} else {$statStr.=$samples[$JNUM]."\t$dir2rd\t".$curStats."\n"; $statStr5.=$samples[$JNUM]."\t$dir2rd\t".$curStats5."\n";}
	}
	if ($redoAssembly){
		system "rm -fr $finalCommAssDir";
	}

	my $boolGenePredOK=0;
	if ($DO_EUK_GENE_PRED){
		$boolGenePredOK = 1 if (-s "$finalCommAssDir/genePred/proteins.bac.shrtHD.faa" || ($pseudoAssembly && -e "$finalCommAssDir/genePred/proteins.bac.shrtHD.faa"));
	} else {
		$boolGenePredOK = 1 if (-s "$finalCommAssDir/genePred/proteins.shrtHD.faa" || ($pseudoAssembly && -e "$finalCommAssDir/genePred/proteins.shrtHD.faa") );
	}
	#die "$boolGenePredOK\n$finalCommAssDir/genePred/proteins.bac.shrtHD.faa\n";

	my $boolAssemblyOK=0;
	#die $finalCommAssDir;
	#quick fix
	if (!-e "$finalCommAssDir/genePred/genes.per.ctg" && -e "$finalCommAssDir/genePred/genes.gff"){my $tmpGene = "$finalCommAssDir/genePred"; 
		system "cut -f1 $tmpGene/genes.gff | sort | uniq -c | grep -v '#' |  awk -v OFS='\\t' {'print \$2, \$1'} > $tmpGene/genes.per.ctg"} 
	$boolAssemblyOK=1 if ($boolGenePredOK && -e "$finalCommAssDir/scaffolds.fasta.filt" && 
			-e "$curOutDir/mapping/$SmplName-smd.bam.coverage.gz" && 
			(-s "$curOutDir/mapping/$SmplName-smd.bam" || -s "$curOutDir/mapping/$SmplName-smd.cram"));
	#die "$boolAssemblyOK $AssemblyGo ass $finalCommAssDir\n";
	my $boolScndMappingOK = 0; my $iix =0;
	foreach my $bwt2outDTT (@bwt2outD){
	#die "$map{$samples[$JNUM]}{mapFinSmpl} $map{$samples[$JNUM]}{assFinSmpl} \n";
	#"$bwt2outD[$i]/$SmplName-0-smd.bam.coverage.gz"
		#my $expectedMapCovGZ = "$bwt2outDTT/$map{$samples[$JNUM]}{mapFinSmpl}"."-0-smd.bam.coverage.gz";
		my $expectedMapCovGZ = "$bwt2outDTT/$bwt2ndMapNmds[$iix]"."_".$SmplName."-0-smd.bam.coverage.gz";
		$iix++;
		print $expectedMapCovGZ."\n";
		if ( -e $expectedMapCovGZ ){
			$boolScndMappingOK=1;
		}else{
			$boolScndMappingOK=0; last;
		}
		#print $expectedMapCovGZ."\n";
		#die "$expectedMapCovGZ";
	}
	
	$boolScndMappingOK = 1 if (@bwt2outD == 0 );#|| !$MappingGo);	
	#die $boolScndMappingOK."\n$boolAssemblyOK\n";
	
	my $KrakenOD = $curOutDir."Tax/kraken/$globalKraTaxkDB/";
	system "rm -r $KrakenOD" if ($RedoKraken && -d $KrakenOD);
	my $calcKraken =0;
	$calcKraken = 1 if ($DoKraken && (!-d $KrakenOD || !-e "$KrakenOD/krakDone.sto"));
	if (!$calcKraken && $DoKraken){
		opendir D, $KrakenOD; my @krkF = grep {/krak\./} readdir(D); closedir D;
		foreach my $kf (@krkF){
			$kf =~ m/krak\.(.*)\.cnt\.tax/; my $thr = $1;# die $thr."  $kf\n";
			system "mkdir -p $dir_KrakFind/$thr" unless (-d "$dir_KrakFind/$thr"); #system "mkdir -p $dir_RibFind/SSU/" unless (-d "$dir_RibFind/SSU/"); system "mkdir -p $dir_RibFind/LSU/" unless (-d "$dir_RibFind/LSU/");
			system "cp $KrakenOD/$kf $dir_KrakFind/$thr/$SmplName.$thr.krak.txt";
		}
	} else {$KrakTaxFailCnts++;}
	
	if ($RedoRiboFind){system "rm -rf $curOutDir/ribos";}
	if ($RedoRiboAssign){system "rm $curOutDir/ribos//ltsLCA/*.sto";}
	#system "rm -f $curOutDir/ribos//ltsLCA/LSU_ass.sto $curOutDir/ribos//ltsLCA/Assigned.sto"; fix for new LSU assignments
	my $calcRibofind = 0; my $calcRiboAssign = 0;
	$calcRibofind = 1 if ($DoRibofind && (!-e "$curOutDir/ribos//ITS_pull.sto"|| !-e "$curOutDir/ribos//SSU_pull.sto"|| !-e "$curOutDir/ribos//LSU_pull.sto" || ($doRiboAssembl && !-e "$curOutDir/ribos/Ass/allAss.sto" )));
	$calcRiboAssign = 1 if ($DoRibofind && ( !-e "$curOutDir/ribos//ltsLCA/ITS_ass.sto"
			|| !-e "$curOutDir/ribos//ltsLCA/Assigned.sto" || !-e "$curOutDir/ribos//ltsLCA/LSU_ass.sto" || !-e "$curOutDir/ribos//ltsLCA/SSU_ass.sto") );
	#die "$calcRibofind $calcRiboAssign\n";
	if ($calcRibofind){
		#system "rm -r $curOutDir/ribos";
		$riboFindFailCnts ++ ;
	} elsif ($DoRibofind &&!$calcRiboAssign) { #copy files to central dir
		my @RFtags = ("ITS","SSU","LSU");
		foreach my $RFtag (@RFtags){
			system "mkdir -p $dir_RibFind/$RFtag/" unless (-d "$dir_RibFind/$RFtag/"); #system "mkdir -p $dir_RibFind/SSU/" unless (-d "$dir_RibFind/SSU/"); system "mkdir -p $dir_RibFind/LSU/" unless (-d "$dir_RibFind/LSU/");
			my $fromCp = "$curOutDir/ribos/ltsLCA/${RFtag}riboRun_bl.hiera.txt"; my $toCpy = "$dir_RibFind/$RFtag/$SmplName.$RFtag.hiera.txt";
			if (!-e $toCpy  || ( (-e $toCpy  || -e "$toCpy.gz" ) && -e $fromCp && -s $fromCp != -s $toCpy)){
				if (-e "$fromCp.gz"){
					system "zcat $fromCp.gz > $toCpy" ;
				} else {
					system "cp $fromCp $toCpy" ;
				}
			}
			system "gzip $fromCp" unless (-e "$fromCp.gz");
			system "rm -f $curOutDir/ribos/ltsLCA/inter${RFtag}riboRun_bl.fna";
		}
	}
	my ($calcDiamond,$calcDiaParse) = IsDiaRunFinished($curOutDir);
	#die "XX $calcDiamond $calcDiaParse\n";
	my $calcMetaPhlan=0;$calcMetaPhlan=1 if ($DoMetaPhlan &&  !-e $dir_MP2."$SmplName.MP2.sto");
	if (!-e $dir_MP2."$SmplName.MP2.sto"){$metaPhl2FailCnts++;}
	if ($redoCS){
		system("rm -r -f $finalCommAssDir/ContigStats/ $curOutDir/assemblies/metag/ContigStats/ $curOutDir/Binning/");
	}
	if ( $boolAssemblyOK ){#causes a lot of overhead but mainly to avoid unpacking reads again..
		print "present: $curOutDir \n"; $present++;
		#base is present, but is the additions? 
		my $CRAMsto = "$curOutDir/mapping/$SmplName-smd.cram.sto";
		#my $boolMappingOK = -e "$curOutDir/mapping//$SmplName-smd.bam.coverage.gz";
		#die $boolMappingOK." gg\n";
		#my ($map2Ctgs,$delaySubmCmd) = mapReadsToRef("$curOutDir/mapping/","",\@cfp1,\@cfp2,$GlbTmpPath."/toMGctgs/",$nodeSpTmpD,$bwt_Cores,
		#	"$finalCommAssDir/scaffolds.fasta.filt",$AsGrps{$cAssGrp}{AssemblJobName},$mapOut."unaligned/",1,"$curOutDir/mapping/",1,$SmplName,\@libsCFP);#$localAssembly);

		if ($doBam2Cram && !-e $CRAMsto){#.cram.sto to check that everything went well
			$QSBopt{LocationCheckStrg}="";
			bam2cram("$curOutDir/mapping/$SmplName-smd.bam", "$finalCommAssDir/scaffolds.fasta.filt",1,1,$doBam2Cram,$CRAMsto,4);
		}
		#if ($AssemblyGo && !-s "$finalCommAssDir/genePred/proteins.shrtHD.faa"){
		#	run_prodigal($metaGassembly,$geneDir,"");
		#}
		#print $AssemblyGo."\n";
		if ($AssemblyGo && (!-s "$finalCommAssDir/ContigStats/scaff.pergene.GC" || !-s "$finalCommAssDir/ContigStats/scaff.4kmer.gz" 
				|| !-s "$curOutDir/assemblies/metag/ContigStats/Coverage.pergene" || !-s "$finalCommAssDir/ContigStats/FMG/proteins.shrtHD.all.marker_genes_scores.table" 
				|| !-s "$curOutDir/assemblies/metag/ContigStats/Coverage.count_pergene" || !-s "$finalCommAssDir/ContigStats/GeneStats.txt") ){
			print "Running Contig Stats on assembly\n";
			#print "$finalCommAssDir/ContigStats/scaff.pergene.GC\n";
			contigStats($curOutDir ,"",$curTmpDir,$finalCommAssDir,"agkems",1,1,$samplReadLength);
		} elsif (!$AssemblyGo && (!-s "$curOutDir/assemblies/metag/ContigStats/Coverage.pergene" || !-s "$curOutDir/Binning/MaxBin/MB.summary" 
									|| !-s "$curOutDir/Binning/MetaBat/MeBa.sto" )){
			print "Running Gene Coverage & Binning\n";
			#die -s "$curOutDir/assemblies/metag/ContigStats/Coverage.pergene" .-s "$curOutDir/Binning/MaxBin/MB.summary" .-s "$curOutDir/Binning/MetaBat/MeBa.sto"."\n"; 
		#my ($jn,$delaySubmCmd2) = contigStats($curOutDir ,$AsGrps{$cAssGrp}{CSfinJobName},$curTmpDir,$finalCommAssDir,"a",$AssemblyGo);
			contigStats($curOutDir ,"",$curTmpDir,$finalCommAssDir,"ams",1,1,$samplReadLength);
		} elsif (!$statsDone) {#ready to collect some stats
			my ($statsHD,$curStats,$statsHD5,$curStats5) = smplStats($curOutDir,$assDir);
			$statsDone=1;
			if ($statStr eq ""){
				$statStr.="SMPLID\tDIR\t".$statsHD."\n".$samples[$JNUM]."\t$dir2rd\t".$curStats."\n";
				$statStr5.="SMPLID\tDIR\t".$statsHD5."\n".$samples[$JNUM]."\t$dir2rd\t".$curStats5."\n";
			} else {$statStr.=$samples[$JNUM]."\t$dir2rd\t".$curStats."\n"; $statStr5.=$samples[$JNUM]."\t$dir2rd\t".$curStats5."\n";}
		}
		#die "$boolScndMappingOK && !$calcD2s && !$calcRibofind && !$calcDiamond && !$calcMetaPhlan && !$calcKraken\n";
		if ($boolScndMappingOK && !$calcD2s && !$calcRibofind && !$calcRiboAssign && !$calcDiamond && !$calcDiaParse && !$calcMetaPhlan && !$calcKraken && $scaffTarExternal eq ""){
			#free some scratch
			system "rm -r $GlbTmpPath/rawRds" if ($DoFreeGlbTmp);
			print "next";next;
		}
	} elsif ($unfiniRew==1){
		print "Deleting previous results..\n";
		system("rm -f -r $curOutDir $GlbTmpPath $collectFinished");
	} elsif (0&&$doSubmit == 0){
		#system("mkdir -p $GlbTmpPath\nmkdir -p $DBpath\nmkdir -p $assDir\n mkdir -p $logDir");
		next;
	}
#die "X";
	
	#next;
	system("mkdir -p $GlbTmpPath\nmkdir -p $DBpath\nmkdir -p $assDir\n mkdir -p $logDir");
	#435590.58253.NC_009614
	$QSBopt{LocationCheckStrg} = checkDrives(\@checkLocs);
#print $finalCommAssDir."\n";

	#--------------------------------------------------------------
	# get the raw fasta files in Jens metagenomic dirs in ordinary structure
	#unzip and change ifastap & cfp1/cfp2
	if (@unzipjobs > 12){#only run 8 jobs in parallel, lest the cluster breaks down..
		#my @last_n = @unzipjobs[-8..-1]; 
		$curUnzipDep = $unzipjobs[-12];#join(";",@last_n);
		$waitTime=0;
	}
	#set up dirs for this sample
	my $metagAssDir = $assDir."metag/";
	my $geneDir = $metagAssDir."genePred/";
	system("mkdir -p $metagAssDir $geneDir");
	my $metaGassembly=$metagAssDir."scaffolds.fasta.filt"; 
	my $finAssLoc = "$finalCommAssDir/scaffolds.fasta.filt";
	my $finalCommScaffDir = "$finalCommAssDir/scaffolds/";
	my $metaGscaffDir = "$metagAssDir/scaffolds/";
	my $finScaffStone = "$finalCommScaffDir/scaffDone.sto";
	my $pseudoAssFile = "$metagAssDir/longReads.fasta.filt";
	my $pseudoAssFileFinal = "$finalCommAssDir/longReads.fasta.filt";
	my $NPD = $curOutDir."nonpareil/";
	#my $finScaffLoc = "";
	my $assemblyDep = "";

#	#-----------------------  FLAGS  ------------------------  
#	#and some more flags for subprocesses
	my $nonPareilFlag = !-s "$NPD/$SmplName.npo" && $DoNP ;
	my $scaffoldFlag = 0; $scaffoldFlag = 1 if (( !-e $finScaffStone) && $map{$samples[$JNUM]}{"SupportReads"} =~ /mate/i );
	my $assemblyFlag = 0; $assemblyFlag = 1 if (!-s $finAssLoc && $DoAssembly);
	my $mapAssFlag = 0; $mapAssFlag = 1 if (!-e "$curOutDir/mapping/$SmplName-smd.bam.coverage.gz" && $assemblyFlag );
	my $pseudAssFlag = 0; $pseudAssFlag = 1 if ($pseudoAssembly && $map{$samples[$JNUM]}{ExcludeAssem} eq "0" && (!-e $pseudoAssFileFinal.".sto" || !$boolGenePredOK));

	my $seqCleanFlag = 0; $seqCleanFlag =1 if (!-e "$GlbTmpPath/seqClean/filterDone.stone" && 
				($scaffTarExternal ne "" || $assemblyFlag  || $pseudAssFlag || $scaffoldFlag || !$boolScndMappingOK || $nonPareilFlag || $calcDiamond || $calcD2s || $calcKraken || $calcRibofind || $calcMetaPhlan) );
#die "$assemblyFlag\n$seqCleanFlag\n";
	my $calcUnzip=0;
	$calcUnzip=1 if ($seqCleanFlag  || $mapAssFlag || !$boolScndMappingOK); #in these cases I need raw reads anyways..
	#die "$calcUnzip = $seqCleanFlag  || $mapAssFlag || !$boolScndMappingOK\n";
	#die "!$boolScndMappingOK || (!-e $GlbTmpPath/seqClean/filterDone.stone && ( $nonPareilFlag || $calcDiamond || $calcD2s || $calcKraken || $calcRibofind || $calcMetaPhlan)) );\n";

#	#-----------------------  END FLAGS  ------------------------  


	#die "$pseudAssFlag\n$boolGenePredOK\n$pseudoAssFileFinal\n";
	if ($scaffTarExternal ne "" &&  $map{$samples[$JNUM]}{"SupportReads"} !~ /mate/i && $scaffTarExtLibTar ne $samples[$JNUM] ){print"scNxt\n";next;}
	#has this sample extra reads (e.g. long reads?)
	#die "$assemblyFlag || !$boolScndMappingOK || $nonPareilFlag || $calcDiamond || $calcD2s || $calcKraken\n";
	my ($jdep,$cfp1ar,$cfp2ar,$cfpsar,$WT,$rawFiles, $mmpuNum, $libInfoRef, $jdepUZ) = 
		seedUnzip2tmp($curDir,$map{$samples[$JNUM]}{SupportReads},$curUnzipDep,$nodeSpTmpD,
		$GlbTmpPath,$waitTime,$importMocat,$AsGrps{$cMapGrp}{CntMap},$calcUnzip);
	if($scaffTarExtLibTar eq $samples[$JNUM]){
		@scaffTarExternalOLib1 = @{$cfp1ar}; @scaffTarExternalOLib2 = @{$cfp2ar};
		next unless ($map{$samples[$JNUM]}{"SupportReads"} =~ /mate/i );
		#print "@scaffTarExternalOLib1\n";
		#die;
	}
		$mmpuOutTab .= $dir2rd."\t".$mmpuNum."\n";
	$waitTime = $WT;
	push (@unzipjobs,$jdepUZ) unless ($jdepUZ eq "");
	$AsGrps{$cMapGrp}{SeqUnZDeps} .= $jdep.";";
	$AsGrps{$cAssGrp}{UnzpDeps} .= $jdep.";";
	$AsGrps{$cAssGrp}{readDeps} = $AsGrps{$cAssGrp}{UnzpDeps};
	my @cfp1 = @{$cfp1ar}; my @cfp2 = @{$cfp2ar}; my @libsCFP = ();#stores the read files used plus which library they come from
	if (@cfp1!=@cfp2){print "Fastap path not of equal length:\n@@cfp1\n@@cfp2\n"; die();}
	#die "@cfp1\n";
	push(@{$AsGrps{$cMapGrp}{RawSeq1}},@cfp1);
	push(@{$AsGrps{$cMapGrp}{RawSeq2}},@cfp2);
	#add info for library used
	my @libInfo = @{$libInfoRef};#local info from unzip step
	push(@libsCFP,@libInfo);
	#idea was to do this check before redoing the unzip step, if filtered already present
	#my $filterReadsPresent = check_sdm_loc($cfp1ar,$cfp2ar,\@libInfo,$GlbTmpPath."seqClean/");

	push(@{$AsGrps{$cMapGrp}{Libs}},@libsCFP);
	open O,">$curOutDir/input_raw.txt"; print O $rawFiles; close O;
	
	#empty links for assembler and nonpareil
	my($arp1,$arp2,$singAr,$matAr,$sdmjN) = ([],[],[],[],"");
	#empty links and objects for merging of reads
	my($mergRdsHr,$mergJbN) = ({},"");
	if (1){ #needs to run to get link to files
		if ($importMocat==1 ){#function actually does nothing
			($arp1,$arp2,$singAr,$sdmjN) = mocat_reorder($cfp1ar,$cfp2ar,$cfpsar,$jdep);
			#my @ta1 = @{$arp1}; my @ta2 = @{$arp2}; die "@{$cfp1ar}XX\n@ta1\n";
			#rewrite ifastp in case sdm takes over from here..
			$cfp1ar = $arp1; $cfp2ar = $arp2;$cfpsar = $singAr;
		} 
		#still punsh the whole thing through sdm..
		if ($useSDM!=0) {
			#$ifastp
			($arp1,$arp2,$singAr,$matAr,$sdmjN) = sdmClean($cfp1ar,$cfp2ar,$cfpsar,\@libInfo,$nodeSpTmpD."sdm/",$GlbTmpPath."mateCln/",
					,$GlbTmpPath."seqClean/",$curSDMopt,$jdep.";$sdmjN") ;
			
		}
		#die "@{$singAr}";
		($mergRdsHr,$mergJbN) = mergeReads($arp1,$arp2,$sdmjN,$GlbTmpPath."merge_clean/",$doReadMerge);
		#my %mrr = %{$mergRdsHr}; die "$mrr{pair1}\n$mrr{pair2}\n$mrr{mrg}\n";
		#raw files only required for mapping reads to assemblies, so delete o/w
		#die "!$DoAssembly && $importMocat";
		if (!$DoAssembly && $importMocat==0 && $removeInputAgain){ $sdmjN = cleanInput($cfp1ar,$cfp2ar,$sdmjN,$GlbTmpPath);}
		$AsGrps{$cAssGrp}{readDeps} .= ";$mergJbN";
	}
	#keeps track of all sdm jobs
	$sdmjNamesAll .= ";".$sdmjN;
	push(@allFilter1,@{$arp1});
	push(@allFilter2,@{$arp2});
	
	
	#%%%%%%%%%%%%%%%%   functions only dependent on FILTERED reads   #%%%%%%%%%%%%%%%%
	
	#take long reads and filter for very long reads (454 etc might be long enough)
	#then do gene predictions *instead* of gene predictions on assemblies
	#die $pseudAssFlag."\n";
	if ($pseudAssFlag ){
		my ($psAssDep, $psFile, $metagDir) = createPsAssLongReads($arp1,$arp2,$singAr, $mergJbN.";".$sdmjN, $pseudoAssFile, $finalCommAssDir, $SmplName);#pseudoAssFileFinal
		if ($psAssDep ne ""){
			$AsGrps{$cAssGrp}{pseudoAssmblDep} = $psAssDep;
		}
		#die "$metagDir\n";
		my $prodRun = run_prodigal_augustus($psFile,$metagDir."/genePred/",$psAssDep,$finalCommAssDir,"",$GlbTmpPath);
		if ($prodRun ne ""){
			$AsGrps{$cAssGrp}{pseudoAssmblDep}  = $prodRun;
		#hijack structure originally used by spadesAssmebly routines..
		}
		push(@{$AsGrps{$cAssGrp}{PsAssCopies}}, $assDir."/metag/*",$finalCommAssDir);
		$AsGrps{$cAssGrp}{readDeps} .= ";$prodRun";
	}
	if ($calcDiamond || $calcDiaParse){
		my ($djname,$djCln) = runDiamond($arp1,$arp2,$singAr,$mergRdsHr,$curOutDir."diamond/",$globaldDiaDBdir,$nodeSpTmpD."/diaRefDB/",
					$mergJbN.";".$sdmjN,$reqDiaDB); #GlbTmpPath
		$AsGrps{$cAssGrp}{DiamDeps} = $djname.";";
		$AsGrps{global}{DiamDeps} .= ";$djname";
		$AsGrps{global}{DiamCln} = $djCln unless($djCln eq "");
		$AsGrps{$cAssGrp}{readDeps} .= ";$djname";
	}
	
	#non pareil (estimate community size etc)
	if ($nonPareilFlag){
		nopareil($arp1,$NPD, $globalNPD, $SmplName,$sdmjN);
		
		#debug
		next;
	}
	
	#kraken (estimate taxa abundance
	if ($calcKraken){
		my $krJdep = krakenTaxEst($arp1,$arp2,$singAr,$KrakenOD, $nodeSpTmpD."krak/",$SmplName,$sdmjN);
		$AsGrps{$cAssGrp}{readDeps} .= ";$krJdep";
	}
	if ($calcRibofind || $calcRiboAssign){
		my $ITSrun = detectRibo($arp1,$arp2,$singAr,$nodeSpTmpD."ITS/",$curOutDir."ribos/",$sdmjN,$SmplName,$runTmpDirGlobal); #GlbTmpPath  \@cfp1,\@cfp2
		#$AsGrps{$cAssGrp}{ITSDeps} .= $ITSrun.";";
		$AsGrps{$cAssGrp}{readDeps} .= ";$ITSrun";
	}
	
	#metaphlan2 - taxa abudnance estimates
	if ($calcMetaPhlan){
		my $MP2jname = metphlanMapping($arp1,$arp2,$nodeSpTmpD."MP2/",$dir_MP2,$SmplName,$bwt_Cores,$sdmjN); #\@cfp1,\@cfp2
		$AsGrps{$cAssGrp}{readDeps} .= ";$MP2jname";
	}
	
	my $SmplNameX = $SmplName;
	if ($AsGrps{$cAssGrp}{CntAss} > 1){$SmplNameX .= "M".$AsGrps{$cAssGrp}{CntAss};}
	#my @tmp = @{$arp1};die ("ASSflag: ".$assemblyFlag."\n@tmp\n");
	if ($assemblyFlag){
		die "Can't do assembly and pseudoassembly on the same sample!\n" if ($pseudAssFlag || $pseudoAssembly);
		#add filtered seq to assembly
		push(@{$AsGrps{$cAssGrp}{FilterSeq1}},@{$arp1});
		push(@{$AsGrps{$cAssGrp}{FilterSeq2}},@{$arp2});
		push(@{$AsGrps{$cAssGrp}{FilterSeqS}},@{$singAr});
		#add jdep to current assembl group
		$AsGrps{$cAssGrp}{SeqClnDeps} .= $sdmjN.";";
		if ($AssemblyGo && $assemblyFlag){
			print "Assembly step";
			#print $AsGrps{$cAssGrp}{AssemblJobName}."\n";		die "@{$AsGrps{$cAssGrp}{FilterSeq1}}"."\n";
			my $hostFilter = 0;
			$hostFilter = 1 if ($AsGrps{$cAssGrp}{CntAimAss} > 3);#reset required HDD space
			my $tmpN = spadesAssembly( \%AsGrps,$cAssGrp,"$nodeSpTmpD/ass",$metagAssDir,$spadesBayHam ,
				$shortAssembly, $SmplNameX,$hostFilter,$scaffoldFlag) ;
			#if mates available, do them here
			postSubmQsub("$logDir/MultiMapper.sh",$AsGrps{$cAssGrp}{PostAssemblCmd},$AsGrps{$cAssGrp}{AssemblJobName},$tmpN);
			$AsGrps{$cAssGrp}{PostAssemblCmd} = "";$AsGrps{$cAssGrp}{AssemblJobName} = $tmpN.";$jdep"; #always add in dep on read extraction
			#die "$AsGrps{$cAssGrp}{AssemblJobName} deps\n";
			push(@{$AsGrps{$cAssGrp}{AssCopies}}, $assDir."/metag/*",$finalCommAssDir);
		}
		#die ("asas $AssemblyGo $assemblyFlag\n");
	} elsif ($AssemblyGo && $AsGrps{$cAssGrp}{PostAssemblCmd} ne ""){
			print "assembly exists, but postassembly jobs unfinished\n";
			#die "POSTCMD: ".$AsGrps{$cAssGrp}{PostAssemblCmd}."\n";
			postSubmQsub("$logDir/MultiMapper.sh",$AsGrps{$cAssGrp}{PostAssemblCmd},$AsGrps{$cAssGrp}{AssemblJobName},"");
			$AsGrps{$cAssGrp}{PostAssemblCmd} = "";$AsGrps{$cAssGrp}{AssemblJobName} = $jdep; #always add in dep on read extraction
			push(@{$AsGrps{$cAssGrp}{AssCopies}}, $assDir."/metag/*",$finalCommAssDir);
			$metaGassembly = $finAssLoc;
#		die $AsGrps{$cAssGrp}{PostAssemblCmd}."\n";
	} else {
		print "No Assembly routines required\n";
		$metaGassembly = $finAssLoc;
		$AssJob = $jdep;
		$AsGrps{$cAssGrp}{AssemblJobName} = $jdep;
	}
	
	
	# 2nd assembly step: scaffolding; maybe move later further down?
	if ($scaffoldFlag){
	#$finalCommScaffDir "$finalCommScaffDir/scaffDone.sto" $finAssLoc $metaGassembly
		my $curAssLoc = $metaGassembly;
		$curAssLoc = $finAssLoc if (-e $finAssLoc);
		#die "@{$AsGrps{$cMapGrp}{RawSeq1}},@{$AsGrps{$cMapGrp}{RawSeq2}},@{$AsGrps{$cMapGrp}{Libs}} scaff\n";
		my ($newScaff,$sdep) = scaffoldCtgs($AsGrps{$cMapGrp}{RawSeq1},$AsGrps{$cMapGrp}{RawSeq2},$AsGrps{$cMapGrp}{Libs},
				[],[],
				$metagAssDir,$nodeSpTmpD,$metaGscaffDir,$AsGrps{$cAssGrp}{AssemblJobName},$bwt_Cores, $SmplNameX,1,"");
		unless ($sdep eq ""){
			$AsGrps{$cAssGrp}{AssemblJobName} = $sdep;
		}
	}
	
	#external contigs to be scaffolded (e.g. TEC2 extracts)
	if ($scaffTarExternal ne ""){
	#die "inscaff\n";
	#die "@scaffTarExternalOLib1\n";
		my $metaGscaffDirExt = "$curOutDir/scaffolds/$scaffTarExternalName/";
		my ($newScaff,$sdep) = scaffoldCtgs($AsGrps{$cMapGrp}{RawSeq1},$AsGrps{$cMapGrp}{RawSeq2},$AsGrps{$cMapGrp}{Libs},
				#\@scaffTarExternalOLib1,\@scaffTarExternalOLib2,
				[],[],
				$scaffTarExternal,$nodeSpTmpD."/SCFEX$scaffTarExternalName/",$metaGscaffDirExt,$sdmjN,$bwt_Cores, $SmplNameX,0,
				$scaffTarExternalName);
	#my($ar1,$ar2,$scaffolds,$GFdir_a) = @_; #.= "_GFI1";
		if (@scaffTarExternalOLib1 > 0 ){
			#die "in gapfill\n$newScaff\n";
			GapFillCtgs(\@scaffTarExternalOLib1,\@scaffTarExternalOLib2,$newScaff,$metaGscaffDirExt."GapFill/",$sdep,$scaffTarExternalName);
		}

		last;
	}
#	die "noscaff";
	
	
	
	#%%%%%%%%%%%%%%%%   functions only dependent on reads   #%%%%%%%%%%%%%%%%
	
	
	#------------   secondary mapping -----------------
	#print "@bwt2outD && $MappingGo && !$boolScndMappingOK\n";
	if (@bwt2outD>0 && $MappingGo && !$boolScndMappingOK){#map reads to specific tar
		#collate strings..
		my @mapOutXS;my @bamBaseNameS;
		for (my $i=0;$i<@bwt2outD; $i++){
			$bwt2outD[$i] =~ m/\/([^\/]+)\/*$/;
			my $smplXDB = $1;
			push @mapOutXS, $GlbTmpPath."xtraMaps/$smplXDB/";
			#my $p1ar = $AsGrps{$cAssGrp}{RawSeq1};
			#bwt2Name
			#$AsGrps{$cMapGrp}{CntMap}
			push @bamBaseNameS, $bwt2ndMapNmds[$i]."_".$SmplName."-0";
		}
		$make2ndMapDecoy{Lib} = $curOutDir;
		my $cramthebam=0;
		#map to all refs at once
		my %dirset = 	(nodeTmp=>$nodeSpTmpD,
						outDir => join(",",@bwt2outD),
						glbTmp => $GlbTmpPath."/xtraMaps/",
						glbMapDir => join(",",@mapOutXS));
		my ($map2CtgsX,$delaySubmCmdX,$mapOptHr) = mapReadsToRef(\%dirset,join(",",@bamBaseNameS), 1,
			$AsGrps{$cMapGrp}{RawSeq1},$AsGrps{$cMapGrp}{RawSeq2},$bwt_Cores,
			join(",",@DBbtRefX),$AsGrps{$cMapGrp}{SeqUnZDeps}.";".$bwt2ndMapDep, "", 1, $cramthebam,
			$SmplName, $AsGrps{$cMapGrp}{Libs});#$localAssembly);
	# ($dirsHr,$outName, $is2ndMap, $par1,$par2,$Ncore,#hm,m,x,x,x,x,x
	#	$REF,$jDepe, $unaligned,$immediateSubm,#m,x,x,x,m
	#	$doCram,$smpl,$libAR) = @_;#x,x,x
		#------------  and calc coverage for each separate
		for (my $i=0;$i<@bwt2outD; $i++){
			 my $mapOutX = $mapOutXS[$i]; 
			$dirset{outDir} = $bwt2outD[$i];$dirset{glbMapDir} =$mapOutXS[$i];
			my ($map2CtgsY,$delaySubmCmdY)  = bamDepth(\%dirset,					#$mapOutXS[$i],$GlbTmpPath."/xtraMaps/",$nodeSpTmpD,
						$bamBaseNameS[$i],$DBbtRefX[$i],$map2CtgsX,$cramthebam,$mapOptHr);
			
			#call abundance on newly made genes
			#prevDeps being gff file
			if ($map2CtgsY eq ""){
				my $empty = calcCoverage("$bwt2outD[$i]/"."$bamBaseNameS[$i]-smd.bam.coverage.gz",$DBbtRefGFF[$i] ,
					$SmplName,$map2CtgsY.";$prevDeps") if ($mapModeCovDo);
			} else {
				$map2CtgsY = calcCoverage("$mapOutX/"."$bamBaseNameS[$i]-smd.bam.coverage.gz",$DBbtRefGFF[$i] ,
					$SmplName,$map2CtgsY.";$prevDeps") if ($mapModeCovDo);
				push(@{$AsGrps{$cMapGrp}{MapCopiesNoDel}},$mapOutX."/*",$bwt2outD[$i]);
			}
			$AsGrps{$cMapGrp}{MapDeps} .= $map2CtgsY.";";



			#my ($jn,$delaySubmCmd2) = contigStats($curOutDir ,$AsGrps{$cAssGrp}{CSfinJobName},$curTmpDir,$finalCommAssDir,"a",$AssemblyGo,0);

			#only used for now for the mapping to spec ref
			my $cln1 = clean_tmp([],[], $AsGrps{$cMapGrp}{MapCopiesNoDel},$AsGrps{$cMapGrp}{MapDeps},"",
				"_mcl$i"."_$JNUM");
			#print "XX $AsGrps{$cMapGrp}{MapDeps}\n";
			$AsGrps{$cMapGrp}{MapCopiesNoDel} = [];$AsGrps{$cMapGrp}{MapDeps}="";
			$AsGrps{$cAssGrp}{scndMapping} .= $cln1.";";#tmp dir shouldn't be deleted before this is done
			#@cleans = {}; #just deactivate clean up for sec mapping..
		}
	} elsif ($mapModeActive) {
		#still needs delays in cleaning command
		$AsGrps{$cMapGrp}{ClSeqsRm} .= ";".$GlbTmpPath;
		next;
	}
	#print "end of ref map \n";next;
	
	#%%%%%%%%%%%%%%%%   functions dependent on assembly -> submit post-assembly   #%%%%%%%%%%%%%%%%
	if ($DoAssembly || $mapAssFlag){
		my $cramthebam = 1;
		my %dirset = 	(nodeTmp=>$nodeSpTmpD,
						outDir => "$curOutDir/mapping/",
						glbTmp => $GlbTmpPath."/toMGctgs/",
						glbMapDir => $mapOut);
						
		my ($map2Ctgs,$delaySubmCmd,$mapOptHr) = mapReadsToRef(\%dirset,$SmplName,0, 
			\@cfp1,\@cfp2,$bwt_Cores,
			$metaGassembly,$AsGrps{$cAssGrp}{AssemblJobName},$mapOut."unaligned/",$AssemblyGo,$cramthebam,
			$SmplName,\@libsCFP);#$localAssembly);
		my ($map2Ctgs_2,$delaySubmCmd_2)  = bamDepth(\%dirset,					#$mapOutXS[$i],$GlbTmpPath."/xtraMaps/",$nodeSpTmpD,
			$SmplName,$metaGassembly,$map2Ctgs,$cramthebam,$mapOptHr);
			$delaySubmCmd .= "\n".$delaySubmCmd_2;

		if ($map2Ctgs_2 ne ""){
			$AsGrps{$cAssGrp}{PostAssemblCmd} .= $delaySubmCmd;
			push(@{$AsGrps{$cAssGrp}{MapCopies}},$mapOut."/*",$curOutDir."/mapping");
			$AsGrps{$cAssGrp}{MapDeps} .= $map2Ctgs_2.";";
		}
	}

	#cleaning of filtered reads & copying of assembly / mapping files
	$AsGrps{$cAssGrp}{ClSeqsRm} .= $curOutDir."seqClean/".";";
	
	
	if ($pseudAssFlag ){
		my @cleans = ();my @copiesNoDels=@{$AsGrps{$cAssGrp}{PsAssCopies}};
		my @copies = ();
		 
		my $cln1="";
		$cln1 = clean_tmp(\@cleans,\@copies, \@copiesNoDels,
			#all dependencies before deleting tmp dirs
			$AsGrps{$cAssGrp}{UnzpDeps}.";".$AsGrps{$cAssGrp}{pseudoAssmblDep},
			$curOutDir,"");#$AsGrps{$cAssGrp}{CSfinJobName}); #.";".$contRun
	}
	
	if ($AssemblyGo && $DoAssembly) {
	
		#call genes, depends on assembly
		#print "DE $AsGrps{$cAssGrp}{AssemblJobName}\n";
		my $geneDirX = $geneDir;
		if (@{$AsGrps{$cAssGrp}{AssCopies}} == 0){$geneDirX = $finalCommAssDir."/genePred/";} #in case assembly already copied over, this needs to be set to copy correctly to finalD
		#die "@{$AsGrps{$cAssGrp}{AssCopies}}\n";
		
		my $prodRun = run_prodigal_augustus($metaGassembly,$geneDirX,$AsGrps{$cAssGrp}{AssemblJobName},$finalCommAssDir,"",$GlbTmpPath);
		#clean / copy results to server
		my @cleans =();
		if ($MappingGo && @bwt2outD==0){ #sync clean of tmp with scnd mapping..
			push(@cleans , $GlbTmpPath);# $GlbTmpPath,$metagAssDir."seqClean/filter*",
		}
		#print $AsGrps{$cAssGrp}{ClSeqsRm}." FF ".$AsGrps{$cMapGrp}{ClSeqsRm}."\n";
		push(@cleans,split(/;/,$AsGrps{$cAssGrp}{ClSeqsRm} ));
		if ($MappingGo){#can only delete at mapping step
			push(@cleans,split(/;/,$AsGrps{$cMapGrp}{ClSeqsRm})) ; $AsGrps{$cMapGrp}{ClSeqsRm} = "";
		}
		my @copies = ();
		#die "@cleans";
		
		push(@copies, @{$AsGrps{$cAssGrp}{MapCopies}},@{$AsGrps{$cAssGrp}{AssCopies}});
		my @copiesNoDels = ();
		#dependent on: (1) all mappings to the assembly finished (2) riboDetect
		#die ("clean dep: ".$AsGrps{$cAssGrp}{MapDeps}.";".$AsGrps{$cAssGrp}{ITSDeps}."\n");
		
		
		#tmp deactivate
		my $cln1="";
		if (1){
			if (!$remove_reads_tmpDir){@cleans =();}
			$cln1 = clean_tmp(\@cleans,\@copies, \@copiesNoDels,
				#all dependencies before deleting tmp dirs
				$AsGrps{$cAssGrp}{MapDeps} . ";".
				$AsGrps{$cAssGrp}{scndMapping}.";".$AsGrps{$cAssGrp}{readDeps}.";".$prodRun,
				$curOutDir,"");#$AsGrps{$cAssGrp}{CSfinJobName}); #.";".$contRun
		}
		$QSBopt{LocationCheckStrg}="";
		
		#clean up assembly groups
		$AsGrps{$cAssGrp}{ClSeqsRm} = ""; @{$AsGrps{$cAssGrp}{MapCopies}} = ();
		$AsGrps{$cAssGrp}{MapDeps} = "";  $AsGrps{$cAssGrp}{readDeps} = ""; $AsGrps{$cAssGrp}{DiamDeps} = "";
		# calc statsitics concercing readqual, mappings, genes & contigs
		my $subprts = "agkse";
		if ($DoBinning){$subprts .= "m";}
		my ($contRun,$tmp33) = contigStats($curOutDir ,$cln1.";".$prodRun,$curTmpDir,$finalCommAssDir,$subprts,1,0,$samplReadLength);

		#run contig stats
		postSubmQsub("$logDir/MultiContigStats.sh",$AsGrps{$cAssGrp}{PostClnCmd},$AsGrps{$cAssGrp}{CSfinJobName},$contRun);
		$AsGrps{$cAssGrp}{PostClnCmd} = "";$AsGrps{$cAssGrp}{CSfinJobName} = $contRun;
	
	} else {
		#calculate solely abundance / gene, has to be run after clean & assembly contigstat step
		my ($jn,$delaySubmCmd2) = contigStats($curOutDir ,$AsGrps{$cAssGrp}{CSfinJobName},$curTmpDir,$finalCommAssDir,"a",$AssemblyGo,0,$samplReadLength);
		$AsGrps{$cAssGrp}{PostClnCmd} .= $delaySubmCmd2;
	}
	#debug
	#if ($AssemblyGo){die();}
	#die();

}

#global clean up cmds (like DB removals from scratch)
DiaPostProcess("",$baseOut);

d2metaDist(\@allSmplNames,\@allFilter1,\@allFilter2,$sdmjNamesAll,$baseOut."/d2StarComp/");

close $QSBopt{LOG};




print "\n$sharedTmpDirP\n".$baseOut."\nFINISHED TAMOC\n";
if ($statStr ne ""){
	open O,">$baseOut/metagStats.txt";
	print O $statStr;
	close O;
	print "Stats in $baseOut/metagStats.txt\n";
}
if ($statStr5 ne ""){
	open O,">$baseOut/metagStats_500.txt";
	print O $statStr5;
	close O;
	#print "Stats in $baseOut/metagStats_500.txt\n";
}


	#my $dir_MP2 = $baseOut."pseudoGC/Phylo/MP2/"; #metaphlan 2 dir
	#my $dir_RibFind = $baseOut."pseudoGC/Phylo/RiboFind/"; #metaphlan 2 dir
	#my $dir_KrakFind = $baseOut."pseudoGC/Phylo/KrakenTax/"; #metaphlan 2 dir
print "Found ".$present." of ".$totalChecked." samples already assembled\n";
if ($DoRibofind ){ if( $riboFindFailCnts){print "$riboFindFailCnts samples with incomplete RiboFind\n";} 
	else {
		print "All samples have assigned RiboFinds\n";
		my $mrgCmd = "";
		my @lvls = ("Domain","Phylum","Class","Order","Family","Genus","Species","Hit2DB");
		$mrgCmd .= "$mergeMiTagScript ".join(",",@lvls)." $dir_RibFind/ITS.miTag $dir_RibFind/ITS/*.hiera.txt \n";
		$mrgCmd .= "$mergeMiTagScript ".join(",",@lvls)." $dir_RibFind/SSU.miTag $dir_RibFind/SSU/*.hiera.txt \n";
		$mrgCmd .= "$mergeMiTagScript ".join(",",@lvls)." $dir_RibFind/LSU.miTag $dir_RibFind/LSU/*.hiera.txt \n";
		die $mrgCmd."\n";
		system $mrgCmd."\n";
		#system "$mrgCmd";
}}
if ($DoKraken ){if( $KrakTaxFailCnts){print "$KrakTaxFailCnts samples with incomplete KrakenTax\n";} 
	else {print "All samples have assigned Kraken Taxonomy\n";
		opendir D, $dir_KrakFind; my @krkF = grep {/0.*/ && -d $dir_KrakFind.$_} readdir(D); closedir D;
		my $mrgCmd = "";
		foreach my $kf (@krkF){
			#$kf =~ m/krak\.(.*)\.cnt\.tax/; my $thr = $1;# die $thr."  $kf\n";
			$mrgCmd .= "$mergeTblScript $dir_KrakFind/$kf/*.krak.txt > $dir_KrakFind/Krak.$kf.mat\n";
		}
		#die $mrgCmd."\n$dir_KrakFind\n";
		system "$mrgCmd";
}}
if (0&&$DoMetaPhlan){if ($metaPhl2FailCnts){print "$metaPhl2FailCnts samples with incomplete Metaphlan assignments\n";}
	else { print "All samples have metaphlan assignments.\n";
		mergeMP2Table($dir_MP2);

}}
open O,">$baseOut/MMPU.txt";print O $mmpuOutTab;close O;

exit(0);
die();







































#--------------------------------------------------------------
#get ref genomes
#print("/g/bork5/hildebra/bin/cdbfasta/cdbfasta $PaulRefGenomes");
#if (!-f $thisRefSeq){system("/g/bork5/hildebra/bin/cdbfasta/cdbyank -a $GID $PaulRefGenomes.cidx > $thisRefSeq");}

#--------------------------------------------------------------
#and the assembly
#$cmd = "/g/bork5/hildebra/dev/Perl/assemblies/./multAss.pl $DBpath2 $assDir2 $assDir2/tmp/ $Ref 4 $CutSeq";

#####################################################
#
#####################################################
sub DiaPostProcess(){
	my ($cmd,$baseOut) = @_;
	#also do higher level summaries of diamond to DB mappings
	return if (1 || ! $DoDiamond ); #temporary deactivated
	$cmd .= "\n\n$mrgDiScr $baseOut\n";
	$cmd .= $AsGrps{global}{DiamCln};
	qsubSystem($globalLogDir."DiaDBcleanup.sh",$cmd ,1,"1G","CDDB",
	$AsGrps{global}{DiamDeps} ,"",1,[],\%QSBopt); #unless ($AsGrps{global}{DiamCln} eq "");

}
sub mergeMP2Table($){
	my ($dir_MP2) = @_;
	my @slvl = ("k__","p__","c__","o__","f__","g__"); my @llvl = ("kingdom","phylum","class","order","family","genus");
	my $mrgCmd = "$mergeTblScript $dir_MP2/*MP2.txt > $dir_MP2/tmp.All.MP2.mat\n";
	for (my $i=0;$i<@slvl;$i++){
		$mrgCmd .= "cat  $dir_MP2/tmp.All.MP2.mat | grep \"".$slvl[$i]."[^\\|]*\\s\" > $dir_MP2/All.$llvl[$i].MP2.mat\n";
	}
	$mrgCmd .= "rm -f $dir_MP2/All.MP2.mat\n\n";
	$mrgCmd .= "$mergeTblScript $dir_MP2/*MP2.noV.noB.txt > $dir_MP2/All.MP2.noV.noB.mat\n";
	for (my $i=0;$i<@slvl;$i++){
		$mrgCmd .= "cat $dir_MP2/All.MP2.noV.noB.mat | grep \"".$slvl[$i]."[^\\|]*\\s\" > $dir_MP2/All.$llvl[$i].MP2.noV.noB.mat\n";
	}
	$mrgCmd .= "rm -f $dir_MP2/All.MP2.noV.noB.mat\n\n";
	$mrgCmd .= "$mergeTblScript $dir_MP2/*MP2.noV.txt > $dir_MP2/All.MP2.noV.mat\n";
	for (my $i=0;$i<@slvl;$i++){
		$mrgCmd .= "cat $dir_MP2/All.MP2.noV.mat | grep \"".$slvl[$i]."[^\\|]*\\s\" > $dir_MP2/All.$llvl[$i].MP2.noV.mat\n";
	}
	$mrgCmd .= "rm -f $dir_MP2/All.MP2.noV.mat\n\n";
	#viruses
	$mrgCmd .= "$mergeTblScript $dir_MP2/*MP2.VirusOnly.txt > $dir_MP2/All.MP2.VirusOnly.mat\n";
	for (my $i=0;$i<@slvl;$i++){
		$mrgCmd .= "cat $dir_MP2/All.MP2.VirusOnly.mat | grep \"".$slvl[$i]."[^\\|]*\\s\" > $dir_MP2/All.$llvl[$i].MP2.VirusOnly.mat\n";
	}
	$mrgCmd .= "rm -f $dir_MP2/All.MP2.noV.mat\n\n";
	
	
	system "$mrgCmd";
	#die $mrgCmd."\n";
}

#calculates the hmm based freq estimates and divergence from these -> used for dist matrix
sub d2metaDist{
	my ($arSmpls,$arPaths1,$arPaths2,$deps,$outPath) = @_;
	
	my @paths = @{$arPaths1}; my @Smpls = @{$arSmpls};
	if(!$calcD2s){return;}
	if (@paths < 1){print "Not enough samples for d2s!\n";return;}
	print "Calculating kmer distances for ".@paths." samples\n";
	system "mkdir -p $outPath/LOGandSUB";
	my $d2MK = 6;
	my $smplFile = "$outPath/mapd2s.txt";
	open O,">$smplFile";
	for (my $i=0;$i<@paths;$i++){
		print O "$paths[$i] $Smpls[$i]\n";
	}
	close O;
	my $cmd = "";
	$cmd .= "cd $outPath\n";
	$cmd .= "$d2metaBin $d2MK $smplFile Q\n"; #or Q for fastQ
	$cmd .= "touch $outPath/d2meta.stone\n";
	my $jobN = ""; my $tmpCmd;
	unless (-e "$outPath/d2meta.stone"){
		$jobN = "_d2met";
		($jobN, $tmpCmd) = qsubSystem($outPath."/LOGandSUB/d2Met.sh",$cmd,1,"80G",$jobN,"$deps","",1,\@General_Hosts,\%QSBopt) ;
	}
	return $jobN;
}


sub postSubmQsub(){#("$logDir/MultiMapper.sh",$AsGrps{$cAssGrp}{PostAssemblCmd},$AsGrps{$cAssGrp}{AssemblJobName},$SmplName);
	my ($outf,$cmdM,$placeh,$jobN) = @_;
	return if ($cmdM eq "");
	print "Replacing Job $placeh with $jobN\n";
	$cmdM =~ s/$placeh/$jobN/g;
	$cmdM =~ s/ -w "done\(\)"//g;
	#die $AsGrps{$cAssGrp}{PostAssemblCmd}."\n";
	open O,">$outf"; print O $cmdM; close O;
	print "$outf\n";
	system "bash $outf";
}

sub detectRibo(){
	my ($ar1,$ar2,$sar,$tmpP,$outP,$jobd,$SMPN,$glbTmpDDB) = @_;
	my $numCore = 12;
	my $numCore2 = 12;
	
	my @re1 = @{$ar1}; my @re2 = @{$ar2}; my @singl = @{$sar};
	#print "ri"; 
	if (@re1 > 1 || @re2 > 1){
		#print"\nWARNING::\nOnly the first read file will be searched for ribosomes\n";
		#first create new tmp file
	}
	
	#copy DB to server
	my $DBrna = "$glbTmpDDB/rnaDB/";
	my $DBrna2 = "$glbTmpDDB/LCADB/";
	if ($globalRiboDependence{DBcp} eq "" ){
		my $DBcmd = "";
		$globalRiboDependence{DBcp}="alreadyCopied";
		if ( !-d $DBrna || !-e "$DBrna/ITS_comb.idx.kmer_0.dat" ||  !-e "$DBrna/silva-euk-18s-id95.fasta"  || !-e "$DBrna/silva-euk-28s-id98.idx.kmer_0.dat"){
			$DBcmd .= "mkdir -p $DBrna\n";
			my $srtMRNA_path = getProgPaths("srtMRNA_path");
			$DBcmd .= "cp $srtMRNA_path/rRNA_databases/silva* $DBrna\n";
			#$DBcmd .= "cp /g/bork3/home/hildebra/DB/MarkerG/ITS_fungi/sh_general_release_30.12.2014.* $DBrna\n";
			my $ITS_DB_path = getProgPaths("ITS_DB_path");
			$DBcmd .= "cp $ITS_DB_path/ITS_comb* $DBrna\n";
			#has to be noted that this doesn't need to happen again
			#print "ribo DB already present\n";
			
			#$DBcmd = "";
		} 
		#and get flash DBs as well over to that dir
		my @LCAdbs = ("/g/bork3/home/hildebra/DB/MarkerG/SILVA/123/SLV_123_LSU_sorted_097_centroids.fasta*","/g/bork3/home/hildebra/dev/lotus//DB//SLV_123_LSU.tax",
			"/g/bork3/home/hildebra/DB/MarkerG/SILVA/123//SLV_123_SSU_sorted_0.97_centroids*", "/g/bork3/home/hildebra/dev/lotus//DB//SLV_123_SSU.tax",
			"/g/bork3/home/hildebra/DB/MarkerG/ITS_combi/ITS_comb.fa*","/g/bork3/home/hildebra/DB/MarkerG/ITS_combi/ITS_comb.tax",
			#"/g/bork3/home/hildebra/dev/lotus//DB//UNITE/sh_refs_qiime_ver7_99_02.03.2015.fasta", "/g/bork3/home/hildebra/dev/lotus//DB//UNITE/sh_taxonomy_qiime_ver7_99_02.03.2015.txt",
			"/g/bork3/home/hildebra/DB/MarkerG/PR2/gb203_pr2_all_10_28_99p.fasta*", "/g/bork3/home/hildebra/DB/MarkerG/PR2/PR2_taxonomy.txt");
		if (!-d $DBrna2 || -s "/g/bork3/home/hildebra/DB/MarkerG/SILVA/123/SLV_123_LSU_sorted_097_centroids.fasta" != -s $DBrna2."/SLV_123_LSU_sorted_097_centroids.fasta"
		|| !-e "$DBrna2/gb203_pr2_all_10_28_99p.fasta" || !-e "$DBrna2/ITS_comb.fa.dna5.fm.sa.val" ){
			$DBcmd .= "mkdir -p $DBrna2\n";
			$DBcmd.= "cp ". join(" ",@LCAdbs) . " $DBrna2";
		}
		my $jN = "_RRDB$JNUM"; my $tmpCmd;
		#die $DBcmd."\n";
		#die $DBcmd."$DBrna/ITS_comb.idx.kmer_0.dat";
		if ($DBcmd ne ""){
			($jN, $tmpCmd) = qsubSystem($logDir."RiboDBprep.sh",$DBcmd,1,"1G",$jN,"","",1,\@General_Hosts,\%QSBopt);
			$globalRiboDependence{DBcp} = $jN;
		}
	}
#	die;
	
	#first detect LSU/SSU/ITS in metag
	my $tmpDY = "$tmpP/SMRNA/";
	my $cmd = "\n\n";
	if ($readsRpairs){
		$cmd .= "$cLSUSSUscript '".join(";",@re1)."' '". join(";",@re2)."' $tmpDY $outP $numCore $SMPN $doRiboAssembl $DBrna\n\n"; #tmpP < scratch too slow
	} else {
		$cmd .= "$cLSUSSUscript '".join(";",@singl)."' '-1' $tmpDY $outP $numCore $SMPN $doRiboAssembl $DBrna\n\n"; #tmpP < scratch too slow
	}
	#this part assigns tax
	my $tmpDX = "$tmpP/LCA/";
	my $numCoreL = $numCore2+3 ;
	my $cmd2 = "$lotusLCA_cLSU $outP $SMPN $numCoreL $DBrna2 $tmpDX $readsRpairs\n\n";
	#die $cmd2."\n";
	my $jobName="";
	my $Scmd = "";
	#die $doRiboAssembl."\n".$cmd."\n";
	#die "$outP/reads_LSU.log\n";
	my $allLCAstones = 0; 
	$allLCAstones = 1 if ( -e "$outP//ltsLCA/LSU_ass.sto" && -e "$outP//ltsLCA/SSU_ass.sto" && -e "$outP//ltsLCA/ITS_ass.sto");
	if (-d $outP && -e "$outP/ITS_pull.sto" && -e "$outP/SSU_pull.sto" && -e "$outP/LSU_pull.sto" && 
		-s "$outP//ltsLCA/LSUriboRun_bl.hiera.txt" && $allLCAstones &&
		($doRiboAssembl && -e $outP."/Ass/allAss.sto" ) ){
		#really everything done
		$jobName = $jobd;
	} else {
		$jobd.=";".$globalRiboDependence{DBcp} unless ($globalRiboDependence{DBcp} eq "alreadyCopied");
		my $tmpCmd; my $mem = "3G";
		if ( !-e "$outP/ITS_pull.sto"|| !-e "$outP/SSU_pull.sto"|| !-e "$outP/LSU_pull.sto" || ($doRiboAssembl && !-e $outP."/Ass/allAss.sto" )){
			$jobName = "_RF$JNUM"; 
			#$mem = "5G" if ($doRiboAssembl);
			($jobName, $tmpCmd) = qsubSystem($logDir."RiboFinder.sh",$cmd,$numCore,$mem,$jobName,$jobd,"",1,[],\%QSBopt);
		} else {
			$jobName = $jobd;
		}
		if (!-e "$outP//ltsLCA/LSUriboRun_bl.hiera.txt" || !-e "$outP//ltsLCA/SSUriboRun_bl.hiera.txt" || !-e "$outP//ltsLCA/ITSriboRun_bl.hiera.txt" || !$allLCAstones ){
			$jobd=$jobName; $mem="3G";
			$jobd .= ";".$globalRiboDependence{DBcp} unless ($globalRiboDependence{DBcp} eq "alreadyCopied");
			($jobName, $tmpCmd) = qsubSystem($logDir."RiboLCA.sh",$cmd2,$numCore2,$mem,"_RA$JNUM",$jobName,"",1,[],\%QSBopt);
		}
	}
	return $jobName;
}
sub IsDiaRunFinished($){
	my ($curOutDir) = @_;
	my @alldbs = split /,/,$reqDiaDB;
	if (!$DoDiamond){return (0,0);}
	if ($rewriteDiamond && @alldbs == $maxReqDiaDB){ system "rm -r $curOutDir/diamond/" if (-d "$curOutDir/diamond/");}
	my $cD = 0; my $pD = 0;
	if ($rewriteDiamond){$redoDiamondParse = 1;}
	if ($redoDiamondParse && @alldbs == $maxReqDiaDB){
		system "rm -r $curOutDir/diamond/CNT*";
	}
	foreach my $term (@alldbs){
		#print $term."   $cD, $pD\n";
		if ($redoDiamondParse ){#&& ( -e "$curOutDir/diamond/dia.$term.blast.gz.stone" || -e "$curOutDir/diamond/dia.$term.blast.srt.gz.stone") ){ 
			system "rm -f $curOutDir/diamond/dia.$term.blast.*.stone" ;
			system "$secCogBin $curOutDir/diamond/XX $term $diaEVal 4";
			system "rm -fr $curOutDir/diamond/ABR/" if ($term eq "ABR");
		}
		system "rm -f $curOutDir/diamond/dia.$term.blast*" if ($rewriteDiamond );
		#die "$rewriteDiamond\n";
#print "$curOutDir/diamond/dia.$term.blast.gz\n";
		if ((!-e "$curOutDir/diamond/dia.$term.blast.gz" && !-e "$curOutDir/diamond/dia.$term.blast.srt.gz")){$cD = 1; }
		$pD = 1 if (!-e "$curOutDir/diamond/dia.$term.blast.gz.stone");#  <- last version always requires .srt.gz
	}
	#die "$cD, $pD\n";
	return ($cD,$pD);
}

sub getSpecificDBpaths($){
	my ($curDB) = @_;
	my $DBpath = "";	my $refDB = ""; my $shrtDB = "";
	#if ($curDB eq "NOG"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/eggNOG10/";	$refDB = "eggnog4.proteins.all.fa"; $shrtDB = $curDB;}
	#elsif ($curDB eq "MOH"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/MohFuncts/"; $refDB = "Extra_functions.fna";$shrtDB = $curDB;}
	#elsif ($curDB eq "CZy"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/CAZy/"; $refDB = "Cazys_2015.fasta";$shrtDB = $curDB;}
	#elsif ($curDB eq "ABR"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/ABR_FORS/"; $refDB = "ardb_and_reforghits.fa";$shrtDB = $curDB;}
	#elsif ($curDB eq "ABRc"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/ABR_Card/"; $refDB = "f11_and_card.faa";$shrtDB = $curDB; }
	#elsif ($curDB eq "KGE"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/KEGG/"; $refDB = "genus_eukaryotes.pep";$shrtDB = $curDB; }
	#elsif ($curDB eq "KGB"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/KEGG/"; $refDB = "species_prokaryotes.pep";$shrtDB = $curDB; }
	#elsif ($curDB eq "ACL"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/Aclame/"; $refDB = "aclame_proteins_all_0.4.fasta";$shrtDB = $curDB; }
	#elsif ($curDB eq "KGM"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/KEGG/"; $refDB = "euk_pro.pep";$shrtDB = $curDB; }
	if ($curDB eq "NOG"){$DBpath = getProgPaths("eggNOG40_path_DB");	$refDB = "eggnog4.proteins.all.fa"; $shrtDB = $curDB;}
	elsif ($curDB eq "MOH"){$DBpath = getProgPaths("Moh_path_DB"); $refDB = "Extra_functions.fna";$shrtDB = $curDB;}
	elsif ($curDB eq "CZy"){$DBpath = getProgPaths("CAZy_path_DB"); $refDB = "Cazys_2015.fasta";$shrtDB = $curDB;}
	elsif ($curDB eq "ABR"){$DBpath = getProgPaths("ABRfors_path_DB"); $refDB = "ardb_and_reforghits.fa";$shrtDB = $curDB;}
	elsif ($curDB eq "ABRc"){$DBpath = getProgPaths("ABRcard_path_DB"); $refDB = "f11_and_card.faa";$shrtDB = $curDB; }
	elsif ($curDB eq "KGE"){$DBpath = getProgPaths("KEGG_path_DB"); $refDB = "genus_eukaryotes.pep";$shrtDB = $curDB; }
	elsif ($curDB eq "KGB"){$DBpath = getProgPaths("KEGG_path_DB"); $refDB = "species_prokaryotes.pep";$shrtDB = $curDB; }
	elsif ($curDB eq "ACL"){$DBpath = getProgPaths("ACL_path_DB");$refDB = "aclame_proteins_all_0.4.fasta";$shrtDB = $curDB; }
	elsif ($curDB eq "KGM"){$DBpath = getProgPaths("KEGG_path_DB"); $refDB = "euk_pro.pep";$shrtDB = $curDB; }
	else {die"Unknown DB for Diamond: $curDB\n";}
	return ($DBpath ,$refDB ,$shrtDB );
}

sub prepDiamondDB($ $ $){#takes care of copying the respective DB over to scratch
	my ($curDB,$CLrefDBD,$ncore) = @_;
	my ($DBpath ,$refDB ,$shrtDB) = getSpecificDBpaths($curDB);
	system "mkdir -p $CLrefDBD" unless (-d $CLrefDBD);
	my $clnCmd = "";
	if ($globalDiamondDependence{$curDB} eq "" ){
		my $DBcmd = ""; 
		# check for pr.db.dmnd
		my $ncoreDB=1;
		unless (-e "$DBpath$refDB.db.dmnd"){
			$ncoreDB = $ncore;
			$DBcmd .= "$diaBin makedb --in $DBpath$refDB -d $DBpath$refDB.db -p $ncoreDB\n";
		}
		unless (-e "$DBpath$refDB.length"){
				$DBcmd .= "$genelengthScript $DBpath$refDB $DBpath$refDB.length\n";
		}
		#$clnCmd .= "rm -rf $CLrefDBD;" if (length($CLrefDBD)>6);
		#idea here is to copy to central hdd (like /scratch)
		#system "rm $CLrefDBD/$refDB.db.dmnd"  if ($rewriteDiamond && -d $CLrefDBD);
		#print " !-e $CLrefDBD/$refDB.db.dmnd && !-e $CLrefDBD/$refDB.length\n ";
		if ( -e "$CLrefDBD/$refDB.db.dmnd" && -e "$CLrefDBD/$refDB.length" 
			#&& ($curDB eq "NOG" && !-e "$CLrefDBD/NOG.members.tsv") &&
			#(($curDB ne "KGB" && $curDB ne "KGM" && $curDB ne "KGE")|| -s "$CLrefDBD/genes_ko.list")
			){
			#has to be noted that this doesn't need to happen again
			$globalDiamondDependence{$curDB}="$shrtDB-1";
		} else {
			$DBcmd .= "mkdir -p $CLrefDBD\n";
			$DBcmd .= "cp $DBpath$refDB.db.dmnd $DBpath$refDB.length $CLrefDBD\n";
			if ($curDB eq "NOG" && !-s "$CLrefDBD/NOG.members.tsv"){
				$DBcmd .= "cp $DBpath/NOG.members.tsv $DBpath/NOG.annotations.tsv $DBpath/all_species_data.txt $CLrefDBD\n";
			}
			if ($curDB eq "CZy" && !-s "$CLrefDBD/MohCzy.tax"){
				$DBcmd .= "cp $DBpath/MohCzy.tax $CLrefDBD\n";
			}
			if (($curDB eq "KGB" || $curDB eq "KGM" || $curDB eq "KGE") && !-s "$CLrefDBD/genes_ko.list"){
				$DBcmd .= "cp $DBpath/genes_ko.list $DBpath/kegg.tax.list $CLrefDBD\n";
			}
		}
		my $jN = "_DIDB$shrtDB$JNUM"; my $tmpCmd;
#		die "$DBcmd";
		if ($DBcmd ne ""){
			($jN, $tmpCmd) = qsubSystem($logDir."DiamondDBprep$shrtDB.sh",$DBcmd,$ncoreDB,"1G",$jN,"","",1,\@General_Hosts,\%QSBopt);
			$globalDiamondDependence{$curDB} = $jN;
		}
		#die $DBcmd."\n";
	}
	return ($refDB,$shrtDB,$clnCmd);
}

sub runDiamond(){
	my ($ar1,$ar2,$sa1,$mrgHshHR,$outD,$CLrefDBD,$tmpP,$jdep,$curDB_o) = @_;
	my @reads1 = @{$ar1}; my @reads2 = @{$ar2}; my @singlRds = @{$sa1};
	my %mrgHsh = %{$mrgHshHR};
	#die();
	
	#my $DBpath = "/g/bork3/home/hildebra/DB/eggNOGv9/";	my $refDB = "eggnogv4.proteins.core_periphery.fa";
	my $clnCmd="";my $jobN2="";
	my $mrgMode = 0;
	$mrgMode =1 if (exists ($mrgHsh{mrg}) && $mrgHsh{mrg} ne "");
	#die;
	system "mkdir -p $outD $tmpP" unless (-d $outD && -d $tmpP);

	foreach my $curDB (split /,/,$curDB_o){
	#print "$curDB";
		my $ncore = $dia_Cores;
		#my ($DBpath ,$refDB ,$shrtDB) = getSpecificDBpaths($curDB);
		my ($refDB,$shrtDB,$clnCmd) = prepDiamondDB($curDB,$CLrefDBD,$ncore);
		my $doInterpret = 1;
		#$doInterpret = 0 if ($shrtDB eq "ABR");


		#run actual diamond
		my @collect = (); my @collectSingl=();
		my $cmd ="mkdir -p $tmpP\n";
		my $useLibArrays = 3;
		if ($mrgMode) {$useLibArrays=4;}
		for (my $kk=0;$kk<$useLibArrays;$kk++){
			my @rds = @singlRds;
			if ($mrgMode){
				if ($kk==1){@rds = $mrgHsh{pair2};}
				if ($kk==2){@rds = $mrgHsh{pair1};}
			} else {
				if ($kk==1){@rds = @reads2;}
				if ($kk==2){@rds = @reads1;}
			}
			if ($kk==3){@rds = $mrgHsh{mrg};}
			my $rdsUsed = 0;
			for (my $ii=0;$ii<@rds;$ii++){
				my $query = $rds[$ii];	
				my $outF = "$tmpP/DiaAssignment.sub.$shrtDB.$kk.$ii";
				#my $tmpcnt = `grep -c '^>' $query`; chomp $tmpcnt
				$cmd .= "$diaBin blastx -f tab --compress 1 --quiet -t $tmpP -d $CLrefDBD$refDB.db -q $query -k 5 -e 1e-4 --sensitive -o $outF -p $ncore\n";
				#$cmd .= "$diaBin view -a $outF.tmp -o $outF -f tab\nrm $outF.tmp.daa\n";
				if ($kk==0 || $kk == 3){#single or ext fragments, doesn't need to be sorted
					push(@collectSingl,$outF.".gz");
				} else {
					push(@collect,$outF.".gz");
				}
			}
		}
		
		#die "$cmd\n";
		
		my $out = $outD."dia.$shrtDB.blast";
		my $outgz = "$out.srt.gz";
		#die "$outgz\n";
		#unzip, sort, zip
		$cmd .= "zcat ".join( " ",@collect) ." | sort -T $tmpP | $pigzBin --stdout -p $ncore > $outgz \n";
		$cmd .= "rm -f ". join( " ",@collect) . "\n";
		if (@collectSingl >= 1){
			#append on gzip, can be done with gzip
			$cmd .= "cat ".join( " ",@collectSingl) ." >> $outgz\nrm -f " . join( " ",@collectSingl) ."\n"; #$out.srt
		}
		$cmd.= "rm -r $tmpP\n";
		#die $cmd."\n";
		my $cmd2 = "$secCogBin $outgz $shrtDB $diaEVal 0 $CLrefDBD/$refDB.length $CLrefDBD $tmpP\n";
		#die $cmd2;
		if ($curDB eq "ABR"){
			$cmd2 = "$KrisABR $outgz $outD/ABR/ABR.genes.txt $outD/ABR/ABR.cats.txt\n";
		}
		#if ($doInterpret == 2){
		#	$cmd2 = "$secCogBin $out.gz $shrtDB $diaEVal\n";
		#}
		#$cmd .= "gzip $query"; #WHY???
		
		#check if the secondary routines for parsing blast out still need to be run
		if ($doInterpret){
			if (!(-e  "$out.gz.stone" || -e  "$out.srt.gz.stone" || -e  "$out.stone")){$doInterpret=1;
			} else {$doInterpret=0;
			}
		}
		my $jobName = $jdep;
		
		my $memu = "7G";my $tmpCmd;
		if (!-d $outD || !(-e "$out" || -e "$out.gz"|| -e "$out.srt.gz") ){ #diamond alignments
			$jobName = "_D$shrtDB$JNUM"; 
			my $globDep = $globalDiamondDependence{$curDB};
			$globDep = "" if ($globalDiamondDependence{$curDB} eq "$shrtDB-1");
			
			($jobName, $tmpCmd) = qsubSystem($logDir."Diamo$shrtDB.sh",$cmd,$ncore,$memu,$jobName,$jdep.";".$globDep,"",1,\@General_Hosts,\%QSBopt);
		} else {
			$jobName = "";
		}
		if ($doInterpret){ #parsing of dia output
			my $jobName2 =  "_DP$shrtDB$JNUM";
			$ncore=1; $memu = "30G";
			($jobName, $tmpCmd) = qsubSystem($logDir."Diamo_parse$shrtDB.sh",$cmd2,$ncore,$memu,$jobName2,$jobName,"",1,\@General_Hosts,\%QSBopt) ;
			#die $logDir."ContigStats.sh\n";
			#die "bo";
		}
		if ($jobN2 eq ""){ $jobN2 = $jobName; } else {$jobN2 .= ";".$jobName;}
	}
	return ($jobN2,$clnCmd);
}




sub nopareil(){
	my ($ar1,$outD,$Gdir,$name,$jobd) = @_;
	my $numCore = 4;
	my @re1 = @{$ar1};

	my $sumOut = "$name.npo";
	my $cmd = "mkdir -p $outD\n";
	$cmd .= "$npBin -s $re1[0] -f fastq -t 20 -m 40000 -b $outD$name -n 10240 -i 0.1 -m 0.2\n";#-t $numCore  -o $sumOut
	$cmd .= "cp $outD/$sumOut $Gdir";
	#R part
	#source('/g/bork5/hildebra/bin/nonpareil/utils/Nonpareil.R');
	#Nonpareil.curve('$outD/$sumOut');
	my $jobName = "_NP$JNUM"; my $tmpCmd;
	if (!-d $outD || !-e "$outD/$sumOut"){
		($jobName,$tmpCmd) = qsubSystem($logDir."NonPar.sh",$cmd,1,"42G",$jobName,$jobd,"",1,\@General_Hosts,\%QSBopt);
	}
	return $jobName;
}
sub contigStats(){
	my ($path,$jobd,$cwd,$assD,$subprts,$immSubm,$force, $readLength) = @_;
	
	my $jobName = ""; $QSBopt{LocationCheckStrg}=""; my $tmpCmd="";
	#return;
	my $cmd = "$sepCtsScript $path $assD $subprts $readLength";
	if ($force || (!-d $path."/assemblies/metag/ContigStats/"  || !-s $path."/assemblies/metag/ContigStats/Coverage.count_pergene" || !-s "$path/assemblies/metag/ContigStats/scaff.4kmer" 
			|| !-s "$path/assemblies/metag/ContigStats/ess100genes/ess100.id.txt" && !-s "$path/assemblies/metag/ContigStats/GeneStats.txt")){
		$jobName = "_CS$JNUM";
		($jobName,$tmpCmd) = qsubSystem($logDir."ContigStats.sh",$cmd,1,"20G",$jobName,$jobd,$cwd,$immSubm,[],\%QSBopt);
		#print $logDir."ContigStats.sh\n";
	} else {
		$jobName = $jobd;
	}
	return ($jobName,$tmpCmd);
}
sub calcCoverage(){
	my ($cov,$gff,$cstNme,$jobd) = @_;
	my $jobName = ""; $QSBopt{LocationCheckStrg}=""; my $tmpCmd="";
	my $cmd = "$readCov_Bin $cov $gff";
	if (!-s $cov.".pergene" ){
		$jobName = "_COV$JNUM";
		$jobName = "$cstNme"."_$JNUM" if ($cstNme ne "");
		($jobName,$tmpCmd) = qsubSystem($logDir."COV$cstNme.sh",$cmd,1,"20G",$jobName,$jobd,0,1,[],\%QSBopt);
		#print $logDir."ContigStats.sh\n";
	} else {
		$jobName = $jobd;
	}
	return ($jobName);
}

sub checkDrives{
	my ($aref) = @_;
	my @locs = @{$aref};
	my $retStr = "\n####### BEGIN file location check ######\n";
	foreach my $llo (@locs){$retStr.="ls -l $llo > /dev/null \n";}
	$retStr.=" \nsleep 3\n"; my $cnt=0;
	foreach my $llo (@locs){ $cnt++;
		$retStr.="if [ ! -d \"$llo\" ]; then echo \'Location $cnt does not exist\'; exit 5; fi\n";
	}
	$retStr .= "####### END file location check ######\n";
	return $retStr;
}
sub adaptSDMopt(){ 
#adaptSDMopt($baseSDMopt,$globalLogDir,$samplReadLength);
	my ($baseSF,$oDir,$RL) = @_;
	my $newSDMf = $oDir."/sdmo_$RL.txt";
	my $nRL = $RL - 10 - int($RL/15);
	open I,"<$baseSF" or die "Can't open sdm opt in:\n $baseSF\n"; my $str = join("", <I>); ;close I;
	if ($RL != 0){
		#$str =~ s/minSeqLength\t\d+/minSeqLength\t$nRL/;
		if ($RL < 50){$str =~ s/maxAccumulatedError\t.*\n/maxAccumulatedError\t0.5\n/;
		} elsif ($RL < 90){$str =~ s/maxAccumulatedError\t\d+\.\d+/maxAccumulatedError\t1.2/;} #really shitty GAII platform
		if ($RL < 90){$str =~ s/TrimWindowWidth	\d+/TrimWindowWidth	8/;$str =~ s/TrimWindowThreshhold	\d+/TrimWindowThreshhold	16/;	
		$str =~ s/maxAmbiguousNT	\d+/maxAmbiguousNT	1/;}
	}
	foreach my $so (keys %sdm_opt){$str =~ s/$so	[^\n]*\n/$so	$sdm_opt{$so}\n/;}
	#die $str."\n";
	open O,">$newSDMf" or die "Can't open new sdm opt out:\n $newSDMf\n"; print O $str; close O;
	#die $newSDMf."\n";
	return $newSDMf;
}

# $sdmjN = cleanInput($cfp1ar,$cfp2ar,$sdmjN,$GlbTmpPath);
 sub cleanInput($ $ $ $){
	my ($cfp1ar,$cfp2ar,$sdmjN,$saveD) = @_;
	my @c1 = @{$cfp1ar}; my @c2 = @{$cfp2ar};
	my $cmd = "";
	for (my $i=0;$i<@c1;$i++){
		if ($c1[$i] =~ m/$saveD/){
			
			$cmd .= "rm -f $c1[$i] $c2[$i]\n";
		} else {
			#die "$c1[$i] =~ m/$saveD/\n";
		}
	}
	#die $cmd."\n";
	my $jobName = $sdmjN;
	if ($cmd ne ""){
		$jobName = "_PC$JNUM"; my $tmpCmd;
		($jobName, $tmpCmd) = qsubSystem($logDir."ClnUnzip.sh",$cmd,1,"1G",$jobName,$sdmjN,"",1,\@General_Hosts,\%QSBopt);
	}
	return $jobName;
 }
 
 
#Gap Filler to refine scaffolds
sub GapFillCtgs(){
	my($ar1,$ar2,$scaffolds,$GFdir_a,$dep,$xtrTag) = @_; #.= "_GFI1";
	my @pa1 = @{$ar1}; my @pa2 = @{$ar2};
	my @inserts;
	my $prefi = "GF";
	my $numCore = 16;
	mkdir($GFdir_a.$prefi);
	my $log = $GFdir_a.$prefi."/GapFiller.log";
	system("mkdir -p $GFdir_a$prefi");
	my $GFlib = ($GFdir_a."GFlib.opt");
	my @libFiles;
	for (my $i=0;$i<@pa1;$i++){
		push(@libFiles, ($pa1[$i].",".$pa2[$i]));
		push(@inserts,450);#default value for our hiSeq runs..
	}
	createGapFillopt($GFlib,\@libFiles,\@inserts);
#GF round 1
	my $cmd = "";
	$cmd .= "mkdir -p $GFdir_a$prefi\n";
	$cmd .= $GFbin . " -l $GFlib -s $scaffolds -m 75 -o 2 -r 0.7 -d 70 -t 10 -g 1 -T $numCore -b ".$prefi." -D $GFdir_a > $log \n";
	$cmd .= "rm -r $GFdir_a$prefi/alignoutput $GFdir_a$prefi/intermediate_results $GFdir_a$prefi/reads\n";
	my $finalGFfile = $GFdir_a."$prefi/$prefi.gapfilled.final.fa";
	#die $cmd."\n";
	if ($cmd ne ""){
		my $tmpCmd;
		($dep, $tmpCmd) = qsubSystem($logDir."GapFill_ext_$xtrTag.sh",$cmd,$numCore,int(80/$numCore)."G","_GFE$JNUM",$dep,"",1,[],\%QSBopt);
	}
	return $dep;
}


#scaffolding via mate pairs
sub scaffoldCtgs(){
	my ($ar1,$ar2,$liar, $xar1,$xar2, $refCtgs, $tmpD1,$outD,$dep,$Ncore,$smplName,$spadesRef,$xtraTag) = @_;
	my @pa1 = @{$ar1}; my @pa2 = @{$ar2}; my @libs = @{$liar};
	my @xpa1 = @{$xar1}; my @xpa2 = @{$xar2};
	return ("",$dep) if (@pa1 ==0 );
	my $tmpD = $tmpD1."/scaff/";
	my $bwtIdx = $refCtgs.".bw2"; my $cmdDB="";
	my $clnCmd = ""; my $spadesDir = ""; my $spadFakeDir = "";
	if ($spadesRef ){#&& -e ("$refCtgs/scaffolds.fasta")){#this is supposed to be the spades dir
		my $oldDir = "spades_ori";
		$spadFakeDir = $refCtgs;
		$clnCmd .= "mkdir -p $refCtgs/../$oldDir/; mv $refCtgs/* $refCtgs/../$oldDir/; mv $refCtgs/../$oldDir $refCtgs\n";
		$spadesDir = "$refCtgs$oldDir"; 
		$clnCmd .= "cp $spadesDir/smpls_used.txt $refCtgs\n";
		$refCtgs .= "$oldDir/scaffolds.fasta";
		$clnCmd .="\ngzip $spadesDir/*";
		$clnCmd .="\ngunzip $refCtgs";
		($cmdDB,$bwtIdx) = buildMapperIdx("$refCtgs",$Ncore,0,$MapperProg);#$Ncore);
	} elsif (!-e $bwtIdx){
		($cmdDB,$bwtIdx) = buildMapperIdx("$refCtgs",$Ncore,0,$MapperProg);
	}
	
	
	#my @rd1; my @rd2;
	my $algCmd = "$clnCmd\n$cmdDB\n";
	$algCmd .= "mkdir -p $tmpD\n";
	my @bams; my $cnt=0; my @insSiz; my @orientations;
	for (my $i=0;$i<@pa1;$i++){
		next unless ($libs[$i] =~ m/mate/i);
		#push @rd1,$rds[$i];push @rd2,$rds[$i+1];
		my $tmpOut = "$tmpD/tmpMateAlign$cnt.bam";
		my $tmpBAM = "$tmpD/tmpMateAlign$cnt.srt.bam";
		$algCmd .= "$bwt2Bin --no-unal --end-to-end -p $Ncore -x $bwtIdx -X $mateInsertLength -1 $pa1[$i] -2 $pa2[$i] | $smtBin view -b -F 4 - > $tmpOut\n";
		$algCmd .= "$smtBin sort -@ $Ncore -T kk -O bam -o $tmpBAM $tmpOut; $smtBin index $tmpBAM\n";
		#$algCmd .= "$novosrtBin --ram 50G -o $tmpBAM -i $tmpOut \n";
		$algCmd .= "rm $tmpOut\n\n";
		push(@bams,$tmpBAM); push(@insSiz,10000);push(@orientations,"fr");
		$cnt++;
		#print "sc  mat\n";
	}
	for (my $i=0;$i<@xpa1;$i++){#assumes paired end reads
		#push @rd1,$rds[$i];push @rd2,$rds[$i+1];
		my $tmpOut = "$tmpD/tmpMateAlign$cnt.bam";
		my $tmpBAM = "$tmpD/tmpMateAlign$cnt.srt.bam";
		$algCmd .= "$bwt2Bin --no-unal --end-to-end -p $Ncore -x $bwtIdx  -1 $xpa1[$i] -2 $xpa2[$i] | $smtBin view -b -F 4 - > $tmpOut\n";
		$algCmd .= "$smtBin sort -@ $Ncore -T kk -O bam -o $tmpBAM $tmpOut; $smtBin index $tmpBAM\n";
		#$algCmd .= "$novosrtBin --ram 50G -o $tmpBAM -i $tmpOut \n";
		$algCmd .= "rm $tmpOut\n\n";
		push(@bams,$tmpBAM);push(@insSiz,500);push(@orientations,"fr");
		$cnt++;
		#print "sc  mat\n";
	}
	
	#die "\n\n\n$algCmd\n\n\n";
	return ("",$dep) if (@bams ==0 );
	#create bowtie2 mapping

	my $zcmd = "-z 5000";  #-z 10000";
	 my $ori = "--orientation ".join(" ",@orientations);
	#for (my $i=1;$i<@bams;$i++){ $ori .= " fr"}#$zcmd .= " 10000";
	my $cmd = "";
	my $bams = join(" ",@bams);
	system "mkdir -p $outD" unless (-d $outD);
	$cmd .= "mkdir -p $outD\n";
	$cmd .= "$besstBin $zcmd $ori -f $bams -o $outD -c $refCtgs\n"; #-q $Ncore <- unstable?
	my $jobName = $dep;
	$cmd .= "rm -r $tmpD\n";
	#cleanup2, fake spades result folder
	if ($spadesDir ne ""){
		$cmd .= "mv $outD/BESST_output/pass1/Scaffolds_pass1.fa $spadFakeDir/scaffolds.fasta\n";
		$cmd .= "$renameCtgScr $spadFakeDir/scaffolds.fasta $smplName\n";
		$cmd .= "$sizFiltScr $spadFakeDir/scaffolds.fasta 400 200\n";
		$cmd .= "$assStatScr -scaff_size 500 $spadFakeDir/scaffolds.fasta > $spadFakeDir/AssemblyStats.500.txt\n";
		$cmd .= "$assStatScr $spadFakeDir/scaffolds.fasta > $spadFakeDir/AssemblyStats.ini.txt\n";
		$cmd .= "$assStatScr $spadFakeDir/scaffolds.fasta.filt > $spadFakeDir/AssemblyStats.txt\n";
		my ($cmdX,$bwtIdxX) = buildMapperIdx("$spadFakeDir/scaffolds.fasta.filt",$Ncore,0,$MapperProg);
		$cmd .= $cmdX."\n";

	}
	my $newScaffFNA = "$outD/BESST_output/pass2/Scaffolds_pass2.fa";
	$cmd .= "$renameCtgScr $newScaffFNA $smplName\n";
	$cmd .= "\ntouch $outD/scaffDone.sto\n" ;
	$cmd = "" if (-e "$outD/scaffDone.sto");
	
#die "scaff cmd $cmd\n";
	if ($cmd ne ""){
		my $tmpCmd;
		($dep, $tmpCmd) = qsubSystem($logDir."BesstScaff_ext_$xtraTag.sh",$algCmd.$cmd,$Ncore,int(50/$Ncore)."G","_BBE$JNUM$xtraTag",$dep,"",1,[],\%QSBopt);
	}
	return ($newScaffFNA,$dep);
}
#preprocess mate pairs (use nxtrim on them and communicate results)
 sub check_mates($ $ $ $ $){
	my ($ar,$ifastasPre,$mateD,$doMateCln,$dep) = @_;
	my @mat = @{$ar};
	my $cmd = "";
	my @sarPre; my $mateC=0;
	my @mates;
	foreach my $matp (@mat){
		system "mkdir -p $mateD" unless (-d $mateD);
		my @mateX = split /,/,$matp;
		push(@sarPre,"$mateD/mate${mateC}.se.fastq.gz");
		push(@mates,"$mateD/mate.${mateC}_R1.unknown.fastq.gz","$mateD/mate.${mateC}_R2.unknown.fastq.gz","$mateD/mate.${mateC}_R1.mp.fastq.gz","$mateD/mate.${mateC}_R2.mp.fastq.gz");
		next if ( -e "$mateD/matesDone.sto");
		#--rf keeps reads in rf; --joinreads joins pe
		$cmd .= "$nxtrimBin --ignorePF --separate -1 $mateX[0] -2 $mateX[1] -O $mateD/mate.$mateC\n";
		#$nxtrimBin --stdout-mp -1 $rd1 -2 $rd2 | $bwaBin mem $refCtgs -p - > out.sam
		$ifastasPre .= ";$mateD/mate.${mateC}_R1.pe.fastq.gz,$mateD/mate.${mateC}_R1.pe.fastq.gz";
		$mateC++;
	}
	#die "@mates SDS\n";
	$cmd .= "touch $mateD/matesDone.sto\n";
	#die $cmd;
	$cmd = "" if ($doMateCln == 0 || $mateC==0);
	my $jobName = $dep;
	if ($cmd ne ""){
		my $tmpCmd;
		($jobName, $tmpCmd) = qsubSystem($logDir."mateClean.sh",$cmd,1,"1G","_NX$JNUM",$dep,"",1,\@General_Hosts,\%QSBopt);
	}
	return ($ifastasPre,\@sarPre,\@mates,$jobName) ;
 }

 
 #preprocess mate pairs (use nxtrim on them and communicate results)
 sub check_matesL($ $ $ $){
	my ($pa1,$pa2,$mateD,$doMateCln) = @_;
	my %ret;
	my $cmd = "";
	$cmd .= "rm -rf $mateD; mkdir -p $mateD\n";
	 my $mateC=0;
	#foreach my $matp (@mat){
	if ($mateC > 0){die "mate pairs only supports single library\n";}
	system "mkdir -p $mateD" unless (-d $mateD);
	#my @mateX = split /,/,$matp;
	#--rf keeps reads in rf; --joinreads joins pe
	$cmd .= "$nxtrimBin --ignorePF --separate -1 $pa1 -2 $pa1 -O $mateD/mate.$mateC\n";
	$cmd .= "rm -f $pa1 $pa2\n";
	#$nxtrimBin --stdout-mp -1 $rd1 -2 $rd2 | $bwaBin mem $refCtgs -p - > out.sam
	$ret{pe1} = "$mateD/mate.${mateC}_R1.pe.fastq.gz"; $ret{pe2} = "$mateD/mate.${mateC}_R1.pe.fastq.gz";
	$ret{se} =  "$mateD/mate${mateC}.se.fastq.gz";
	$ret{un1} = "$mateD/mate.${mateC}_R1.unknown.fastq.gz"; $ret{un2} = "$mateD/mate.${mateC}_R2.unknown.fastq.gz";
	$ret{mp1} = "$mateD/mate.${mateC}_R1.mp.fastq.gz"; $ret{mp2} = "$mateD/mate.${mateC}_R2.mp.fastq.gz";
	$mateC++;
	#}
	#die "@mates SDS\n";
	my $locStone = "$mateD/matesDone.sto";
	$cmd .= "touch $locStone\n";
	#die $cmd;
	$cmd = "" if ($doMateCln == 0 || $mateC==0);
	
	return (\%ret,$cmd,$locStone) ;
 }

 
 #determines what read types are present and starts cleaning of mates, if required
 sub get_ifa_mifa($ $ $ $ $ $){
	my ($ifastas,$libInfoAr, $mateD, $doMateCln, $jdep, $singlReadMode) = @_;
	my @allFastas = split /;/,$ifastas;
	my @libInfo = @{$libInfoAr}; my @miSeqFastas = (); my $mcnt=0;
	my @mates;
	foreach (@libInfo){
		if ($_ =~ m/.*miseq.*/i) {
			push (@miSeqFastas, $allFastas[$mcnt]);
			splice(@allFastas, $mcnt, 1);
		} elsif ($_ =~ m/.*mate.*/i){ #remove from process
			push (@mates, $allFastas[$mcnt]);
			splice(@allFastas,$mcnt,1);
		}
		$mcnt++;
	}
	$ifastas = join(";",@allFastas);
	my $mi_ifastas = join(";",@miSeqFastas);
	my $singleAddAr = []; my $matRef = [];
	#moved to seedUnzip2tmp
	#($ifastas,$singleAddAr,$matRef,$jdep) = check_mates(\@mates,$ifastas,$mateD,$doMateCln,$jdep);

	return ($ifastas, $mi_ifastas, $singleAddAr, $matRef, $jdep);
 }
 sub get_sdm_outf($ $ $ $){
	my ($ifastas, $mi_ifastas,$finD,$singlReadMode) = @_;
	my @ret1;my @ret2;my @sret;
 	if ($mi_ifastas ne ""){
		push(@ret1,$finD."filtered_mi.1.fastq");push( @ret2, ($finD."filtered_mi.2.fastq"));push(@sret, ($finD."filtered_mi.singl.fastq"));
	}
	if ($ifastas ne "" && !$singlReadMode){
		push(@ret1,$finD."filtered.1.fastq");push( @ret2, ($finD."filtered.2.fastq"));push(@sret, ($finD."filtered.singl.fastq"));
	} elsif ($ifastas ne "" && $singlReadMode){
		push(@sret, $finD."filtered.s.fastq");
	}
	
	return (\@ret1,\@ret2,\@sret);
}
 sub check_sdm_loc(){
	die "defunct functon  check_sdm_loc!\n";
	my ($ar1,$ar2,$libInfoAr,$sdmO) = @_;
 	my $ifastasPre = ${$ar1}[0].",".${$ar2}[0];
	for (my $i=1;$i<@{$ar1};$i++){if ($i>0){$ifastasPre .= ";".${$ar1}[$i].",".${$ar2}[$i];}}
	my ($ifastas, $mi_ifastas, $singlAddAr, $matAr, $jdep) = get_ifa_mifa($ifastasPre,$libInfoAr,"",0,"","???");
	my ($ar1x,$ar2x,$sar) = get_sdm_outf($ifastas, $mi_ifastas,$sdmO,"???");
	my @ret1=@{$ar1x};my @ret2=@{$ar2x};my @sret=@{$sar};
	my $presence=1;
	for (my $i=0;$i<@ret1;$i++){	if (!-e $ret1[0] || -z $ret1[0]){$presence=0;}	}
	for (my $i=0;$i<@ret2;$i++){	if ( !-e $ret2[0]|| -z $ret1[0]){$presence=0;}	}
	if (!-e $sret[0]){$presence=0;}
	my $stone = "$sdmO/filterDone.stone";
	if (!-e $stone){$presence=0;}
	if (!$presence){return 0;}
	my $assInputFlaw = 0;
	if (-e $logDir."spaderun.sh.otxt"){open I,"<$logDir/spaderun.sh.otxt" or die "Can't open old assembly logfile $logDir\n"; my $str = join("", <I>); close I;
		if ($str =~ /Assertion `seq_\.size\(\) == qual_\.size\(\)' failed\./){$assInputFlaw=1;}	
		if ($str =~ /paired_readers\.hpp.*Unequal number of read-pairs detected in the following files:/){$assInputFlaw=1;}	
	}

	if ($presence==0 || $assInputFlaw==1){
		return 0;
	} else {
		return 1;
	}
}

#($mergRdsHsh,$mergJbN) = mergeReads($arp1,$arp2,$sdmjN,$GlbTmpPath."merge_clean/");
sub mergeReads(){
	my ($arp1,$arp2,$jdep,$outdir,$doMerge) = @_;
	my $numCores = 8;
	my $outT = "sdmCln";
	my %ret = (mrg => "", pair1 => "", pair2 => "");
	if (!$doMerge || !$readsRpairs){return (\%ret,"");}
	if (@{$arp1} > 1 ){die "Array with reads provided to merging routine is too large!\n@{$arp1}\n";}
	my $mergCmd = "$flashBin -M 250 -o $outT -d $outdir -t $numCores ${$arp1}[0] ${$arp2}[0]\n";
	my $stone = "$outdir/$outT.sto";
	$mergCmd .="touch $stone\n ";
	my $jobName = "";
	if (!-e $stone){
		$jobName = "_FL$JNUM"; my $tmpCmd;
		($jobName, $tmpCmd) = qsubSystem($logDir."flashMrg.sh",$mergCmd,$numCores,"3G",$jobName,$jdep,"",1,[],\%QSBopt);
	}
	
	#die "$logDir/flashMrg.sh\n$mergCmd\n";
	$ret{mrg} = "$outdir/$outT.extendedFrags.fastq";
	$ret{pair1} = "$outdir/$outT.notCombined_1.fastq";
	$ret{pair2} = "$outdir/$outT.notCombined_2.fastq";
	return (\%ret,$jobName);
}
 
sub sdmClean(){
	my ($ar1,$ar2,$ars,$libInfoAr,$tmpD,$mateD,$finD,$sdmO,$jobd) = @_;
	#create ifasta path
	#die "cllll\n";
	if (@{$ar1} == 0 && @{$ars}==0){die "Empty array to sdmClean given\n";}
	my $ifastasPre=""; my $singlReadMode=0;
	if (@{$ar1} != 0){
		$ifastasPre = ${$ar1}[0].",".${$ar2}[0];
		for (my $i=1;$i<@{$ar1};$i++){if ($i>0){$ifastasPre .= ";".${$ar1}[$i].",".${$ar2}[$i];}}
	} else {
		$ifastasPre = ${$ars}[0];
		for (my $i=1;$i<@{$ars};$i++){if ($i>0){$ifastasPre .= ";".${$ars}[$i];}}
		$singlReadMode=1;
	}
	open O,">$curOutDir/input_fil.txt"; print O $ifastasPre; close O;

	system("mkdir -p $finD");
	my $cmd = "";
	$cmd .= "mkdir -p $tmpD\n" unless ($tmpD eq "");
	$cmd .= "mkdir -p $finD\n" unless ($tmpD eq "");
	
	my ($ifastas, $mi_ifastas, $singlAddAr, $matAr, $jdep) = get_ifa_mifa($ifastasPre,$libInfoAr,$mateD,0,$jobd,$singlReadMode);
	$jobd = $jdep; #replaces old dep
	my ($ar1x,$ar2x,$sar) = get_sdm_outf($ifastas, $mi_ifastas,$finD,$singlReadMode);
	
	push(@{$sar},@{$singlAddAr});
	my @ret1=@{$ar1x};my @ret2=@{$ar2x};my @sret=@{$sar};
	#die "$mi_ifastas\n$ifastas XX\n";
	#$baseSDMoptMiSeq;
	my $ofiles = "";
	my $mi_ofiles = "";
	system("mkdir -p $logDir/sdm/");
	my $useLocalTmp = 1;
	if ($tmpD eq ""){$useLocalTmp = 0;$tmpD = $finD;}
	my $nsdmPair = 2;
	my $sdm_def= " -ignore_IO_errors 1 -i_qual_offset auto -binomialFilterBothPairs 1 ";
	if ($ifastas ne ""){
		if ($singlReadMode){
			$ofiles = $tmpD."/filtered.s.fastq"; $nsdmPair=1;
		}else {
			$ofiles = $tmpD."/filtered.1.fastq,".$tmpD."/filtered.2.fastq";
		}
		$cmd .= "touch $tmpD/filtered.1.fastq $tmpD/filtered.1.fastq.singl  $tmpD/filtered.2.fastq $tmpD/filtered.2.fastq.singl\n" if (!$singlReadMode);
		$cmd .= "$sdmBin -i \"$ifastas\" -o_fastq $ofiles -options $sdmO -paired $nsdmPair  -log $logDir/sdm/filter.log $sdm_def\n\n";
		#push(@ret1,$finD."filtered.1.fastq");push( @ret2, ($finD."filtered.2.fastq"));push(@sret, ($finD."filtered.singl.fastq"));
		$cmd .= "cat  $tmpD/filtered.1.fastq.singl  $tmpD/filtered.2.fastq.singl >  $sret[-1]\n" if (!$singlReadMode);
	}
	if ($mi_ifastas ne ""){
		$mi_ofiles =  $tmpD."/filtered_mi.1.fastq,".$tmpD."/filtered_mi.2.fastq";
		$cmd .= "touch $tmpD/filtered_mi.1.fastq.singl  $tmpD/filtered_mi.2.fastq.singl\n";
		$cmd .= " $sdmBin -i \"$mi_ifastas\" -o_fastq $mi_ofiles -options $baseSDMoptMiSeq -paired 2  -log $logDir/sdm/filter_xtra.log $sdm_def\n\n";
		#$cmd .= "/g/bork3/home/hildebra/dev/C++/sdm/./sdm -i_fastq \"$mi_ifastas\" -o_fastq $ofiles -options $baseSDMoptMiSeq -paired 2 -i_qual_offset auto -log $logDir/sdm/filter.log -binomialFilterBothPairs 1\n\n";
		#push(@ret1,$finD."filtered_mi.1.fastq");push( @ret2, ($finD."filtered_mi.2.fastq"));push(@sret, ($finD."filtered_mi.singl.fastq"));
		$cmd .= "cat  $tmpD/filtered_mi.1.fastq.singl  $tmpD/filtered_mi.2.fastq.singl > $sret[-1]\n";
	}
	#die $cmd."\n";

	my $stone = "$finD/filterDone.stone";
	if ($useLocalTmp){
		$cmd .= "mv -f ".join(" ",(split(/,/,$ofiles),split(/,/,$mi_ofiles)))." $finD\n";
		#mv -f $tmpD/filtered.2.fastq $finD\n";
		$cmd .= "rm -rf $tmpD\n";
	}
	$cmd .= "touch $stone\n";
	#die "$cmd\n";
	my $jobName = "";
	my $presence=1;
	for (my $i=0;$i<@ret1;$i++){	if (!-e $ret1[0] || -z $ret1[0]){$presence=0;}	}
	for (my $i=0;$i<@ret2;$i++){	if ( !-e $ret2[0]|| -z $ret1[0]){$presence=0;}	}
	if (!-e $sret[0]){$presence=0;}
	if (!-e $stone){$presence=0;}
	#common error: spade input failed
	my $assInputFlaw = 0;
	if (-e $logDir."spaderun.sh.otxt"){open I,"<$logDir/spaderun.sh.otxt" or die "Can't open old assembly logfile $logDir\n"; my $str = join("", <I>); close I;
		if ($str =~ /Assertion `seq_\.size\(\) == qual_\.size\(\)' failed\./){$assInputFlaw=1;}	
		if ($str =~ /paired_readers\.hpp.*Unequal number of read-pairs detected in the following files:/){$assInputFlaw=1;}	
	}
	#die("ass inp flaw: $assInputFlaw\n");

	if ($presence==0 || $assInputFlaw==1){
		$jobName = "_SDM$JNUM"; my $tmpCmd;
		($jobName, $tmpCmd) = qsubSystem($logDir."sdmReadCleaner.sh",$cmd,1,"1G",$jobName,$jobd.";".$jdep,"",1,\@General_Hosts,\%QSBopt);
	}
	
	return (\@ret1,\@ret2,\@sret,$matAr,$jobName);
}
#$importMocat==1
sub mocat_reorder(){
	my($ar1,$ar2,$singlAr,$inJob) = @_;
	
	#my @pairs = split(";",$ifastas);
	my @ret1 = @{$ar1}; my @ret2 = @{$ar2}; 
	#print "@ret1\n@ret2\n";
	#foreach (@pairs){		my @spl = split /,/;		push(@ret1,$spl[0]);push(@ret2,$spl[1]);	}
	#my $jobName = $inJob;
	return (\@ret1,\@ret2,$singlAr,$inJob);
}
sub SEEECER(){
	my ($p1ar,$p2ar,$tmpD) = @_;
	my @p1 = @{$p1ar}; my @p2 = @{$p2ar};
	if ( $p1[0] =~ m/.*\.gz/){
		system("gunzip -c ".join(" ",@p1)." > $tmpD/pair.1.fastq");	system("gunzip -c ".join(" ",@p2)." > $tmpD/pair.2.fastq");
		@p1 = ("$tmpD/pair.1.fastq"); @p2 = ("$tmpD/pair.2.fastq");
	} elsif ( @p1 > 1 ){
		system("cat ".join(" ",@p1)." > $tmpD/pair.1.fastq");	system("cat ".join(" ",@p2)." > $tmpD/pair.2.fastq");
		@p1 = ("$tmpD/pair.1.fastq"); @p2 = ("$tmpD/pair.2.fastq");
	}
	my $SEEbin = "bash /g/bork5/hildebra/bin/SEECER-0.1.3/SEECER/bin/run_seecer.sh";
	system("mkdir -p $tmpD/tmpS");
	my $cmd = $SEEbin . " -t $tmpD/tmpS $p1[0] $p2[0]";
	die("SEECER deactive\n");
	#qsubSystem($logDir."SEECERCleaner.sh",$cmd,1,"30G",1);
	my @ret1= ($tmpD."pair.1.fastq_corrected.fa");my @ret2= ($tmpD."pair.2.fastq_corrected.fa");
	return (\@ret1,\@ret2);
}
sub complexGunzCpMv($ $ $ $){
	my ($fastap,$in,$GlbTmpPath,$finDest) = @_;
	my $unzipcmd = "";
	
	
	if ($in =~ m/\.gz$/){
		$unzipcmd .= "gunzip -c $fastap$in";
		$in =~ s/\.gz$//;
		$unzipcmd .= " > $finDest$in \n";
		
	} else {
		$unzipcmd .= "cp $fastap$in $finDest\n";
	}
	if (0 && $GlbTmpPath ne $finDest){
		$unzipcmd .= "mv -f $GlbTmpPath$in $finDest\n";
		die "$GlbTmpPath ne $finDest";
	}
	my $newRDf = $finDest."$in";
	$unzipcmd .= "chmod +w $newRDf\n";
	return ($unzipcmd,$newRDf);
}

sub seedUnzip2tmp(){
	my ($fastp,$xtrMapStr,$jDepe,$tmpPath,$finDest, $WT,$mocatImport,$libNum,$calcUnzp) = @_;
	my $doMateCln = 1; #nxtrim  mate pairs ?
	my $rawReads=""; my $mmpu = "";
	my $fastp2 = ""; my $xtraRdsTech = "";
	$xtrMapStr = "" if (!defined $xtrMapStr) ;
	if ( $xtrMapStr ne ""){
		my @spl = split(/:/,   $map{$samples[$JNUM]}{"SupportReads"}   );
		$xtraRdsTech = $spl[0].$libNum;
		$fastp2 = $spl[1];
		#die $fastp2."\n".$xtraRdsTech."\n";
	}
	


	#if ($fastp ne "") {print "Looking for fq.gz in ".$fastp;} else { print "No primary dir.. "; }
	if ($fastp2 ne ""){ print " and $fastp2";}
	print "\n";
	my @pa1; my @pa2; my @pas; my @paX1; my @paX2;
	if ($fastp ne ""){
		if ($mocatImport==0){
			opendir(DIR, $fastp) or die "$!: $fastp\n";	
			if ($readsRpairs){
				@pa2 = sort ( grep { /$rawFileSrchStr2/ && -e "$fastp/$_" } readdir(DIR) );	rewinddir DIR;
				@pa1 = sort ( grep { /$rawFileSrchStr1/  && -e "$fastp/$_"} readdir(DIR) );	close(DIR);
			}
			@pas = sort ( grep { /$rawFileSrchStrSingl/  && -e "$fastp/$_"} readdir(DIR) );	close(DIR);
		} else {
			$fastp .= $map{mocatFiltPath}; 
			opendir(DIR, $fastp) or die "$!\n$fastp\n";;	
			@pa2 = sort ( grep { /.*2\.f[^\.]*q\.gz$/ && -e "$fastp/$_" } readdir(DIR) );	rewinddir DIR;
			@pa1 = sort ( grep { /.*1\.f[^\.]*q\.gz$/  && -e "$fastp/$_"} readdir(DIR) );	rewinddir DIR;
			@pas = sort ( grep { /.*single\.f[^\.]*q\.gz$/  && -e "$fastp/$_"} readdir(DIR) );	close(DIR);
		}
	}
	#print "@pa1\n";
	#import xtra long rds
	if (@pa1 != @pa2 && $readsRpairs){die "For dir $fastp, unequal fastq files exist:\nP1\n".join("\n",@pa1)."\nP2:\n".join("\n",@pa2)."\n";}
	#check for date of files (special filter)
	if ($doDateFileCheck){
		my @pa1t = @pa1; my @pa2t = @pa2; undef @pa1; undef @pa2;
		for (my $i=0; $i<@pa1t;$i++){
			my $date = POSIX::strftime( "%m", localtime( ( stat "$fastp/$pa1t[$i]" )[9] ) );
			#print $date;
			if ($date < 3){push(@pa1, $pa1t[$i]);push(@pa2, $pa2t[$i]);}
		}
		if (@pa1 != 1){die "still too many file: $fastp\n @pa1\n@pa2\n";}
	}
	my @libInfo = ("lib$libNum") x int @pa1;
	if (@pa1==0 && @pas!=0){#case of just single read in metag
		@libInfo = ("lib$libNum") x int @pas;
	}
	if ($fastp2 ne ""){
		opendir(DIR, $fastp2) or die $!;	
		@paX2 = sort ( grep { /$rawFileSrchStrXtra2/ && -e "$fastp2/$_" } readdir(DIR) );	rewinddir DIR;
		@paX1 = sort ( grep { /$rawFileSrchStrXtra1/  && -e "$fastp2/$_"} readdir(DIR) );	close(DIR);
		if (@paX2 != @paX1){die "For dir $fastp, unequal fastq files exist:P1\n".join("\n",@paX1)."\nP2:\n".join("\n",@paX2)."\n";}
		if (@paX2 == 0){die "Can't find files with pattern $rawFileSrchStrXtra2 in $fastp2\n";}
		#die $paX1[0]."\n";
		push @pa1 , @paX1; push @pa2 , @paX2; my @lib2 = $xtraRdsTech x @paX1; push @libInfo, @lib2;
	}
		#die "@libInfo\n$libInfo[1]\n";

	#die "@pa1";
	#die "For dir $fastp, unual fastq files exist:P1\n".join("\n",@pa1)."\nP2:\n".join("\n",@pa2)."\n";
	
	#check if symlink and if this is valid & create raw read link:
	for (my $i =0; $i<@pa1; $i++){
		my $pp = $fastp;
		$pp = $fastp2 if ($libInfo[$i] eq $xtraRdsTech);
		if (-l $pp.$pa1[$i]){
			if (!-f abs_path($pp.$pa1[$i])){die "File $pa1[$i] is not file.\n";}
		}
		if (-l $pp.$pa2[$i]){
			if (!-f abs_path($pp.$pa2[$i])){die "File $pa2[$i] does not exist.\n";}
		}
		if ($i==0){
			$rawReads="$pp$pa1[$i],$pp$pa2[$i]";
		} else {
			$rawReads.=";".$pp.$pa1[$i].",".$pp.$pa2[$i];
		}
		my $realP = `readlink -f $pp$pa1[$i]`;
		if ($realP =~ m/\/(MMPU[^\/]+)\//){
			$mmpu = $1;
		}
	}
	#die $rawReads."\n";
	#DEBUG fix to reduce file sizes
	#@fastap2 = ($fastap2[0]); @fastap1 = ($fastap1[0]);
	if (@pa1 == 0 && @pas ==0 ){die"Can;t find files in $fastp\n";}
	my $finishStone = "$finDest/rawRds/done.sto";
	$tmpPath.="/rawRds/";
	my $unzipcmd = "";
	$unzipcmd .= "sleep $WT\nset -e\n";
	$unzipcmd .= "rm -f $finishStone\nrm -r -f $tmpPath\nmkdir -p $tmpPath\nmkdir -p $finDest/rawRds/\n";

	#make sure input is unzipped
	#die ($tmpPath."\n");
	system("mkdir -p $finDest/rawRds/");# my $newRDf ="";
	my $testf1 = "";my $testf2 = "";
	for (my $i=0; $i<@pa1; $i++){
		#print $pa1[$i]."\n";
		my $pp = $fastp;
		$pp = $fastp2 if ($libInfo[$i] eq $xtraRdsTech);
		my ($tmpCmd,$newF) = complexGunzCpMv($pp,$pa1[$i],$tmpPath,$finDest."/rawRds/");
		$unzipcmd .= $tmpCmd."\n";
		$pa1[$i] = $newF;
		($tmpCmd,$newF) = complexGunzCpMv($pp,$pa2[$i],$tmpPath,$finDest."/rawRds/");
		$unzipcmd .= $tmpCmd."\n";
		$pa2[$i] = $newF;
	}
	for (my $i=0; $i<@pas; $i++){
		my $pp = $fastp;
		$pp = $fastp2 if ($libInfo[$i] eq $xtraRdsTech);
		my ($tmpCmd,$newF) = complexGunzCpMv($pp,$pas[$i],$tmpPath,$finDest."/rawRds/");
		$unzipcmd .= $tmpCmd."\n";
		$pas[$i] = $newF;
		if ($splitFastaInput != 0){ #in case of input assemblies (MG-RAST.. arghh!!)
			$unzipcmd .= "\n$sizSplitScr $pas[$i] $splitFastaInput\n";
		}
	}
	$unzipcmd .= "touch $finishStone\n";
	my $jobN = "";
	my $jobNUZ = $jobN;
	
	
	#die $unzipcmd."\n";
	#print "  HH ".-s $testf2 < -s $testf1." FF \n";
	my $tmpCmd;
	if ($calcUnzp){ #submit & check for files
		my $presence=1;
		for (my $i=0;$i<@pa1;$i++){	if (!-e $pa1[0] || -z $pa1[0]){$presence=0;}	}
		for (my $i=0;$i<@pa2;$i++){	if ( !-e $pa2[0] || -z $pa2[0]){$presence=0;}	}
		if (!$presence || !-e $finishStone ){
			$jobN = "_UZ$JNUM"; 
			my $tmpSHDD = $QSBopt{tmpSpace};
			$unzipcmd = "" if ($presence && -e $finishStone);
			
			$QSBopt{tmpSpace} = "150G"; #set option how much tmp space is required, and reset afterwards
			($jobN, $tmpCmd) = qsubSystem($logDir."UNZP.sh",$unzipcmd,1,"20G",$jobN,$jDepe,"",1,\@General_Hosts,\%QSBopt) ;
			#### 1 : UNZIP
			$QSBopt{tmpSpace} = $tmpSHDD;
			$WT += 30;
			#print " FDFS ";
		}
	#	die($jobN);
		$jobN = krakHSap(\@pa1,\@pa2, \@pas,$tmpPath,$jobN);
		#### 2 : remove human contamination
	}
	
		#check already here for mate pair support reads, deactivate fastp2
	my ($mateCmd,$mateSto) = ("","/sh");
#	print "@libInfo\n";
	my $ii=0; my @matePrps;
	while($ii<@libInfo){
		#print $ii." \n";
		#next;
		unless ($libInfo[$ii] =~ m/.*mate.*/i){$ii++;next;} #remove from process
		my $href;
		($href,$mateCmd,$mateSto) = check_matesL($pa1[$ii],$pa2[$ii],$finDest."mateCln/",$doMateCln);
		#remove this ori file from raw reads
		splice(@pa1,$ii,1);splice(@pa2,$ii,1);splice(@libInfo,$ii,1);
		push(@matePrps,$href);
	}
	if ( ($mateSto ne "/sh" && !-e $mateSto) && $calcUnzp){
		$mateCmd = "" if ( -e $mateSto);
		($jobN, $tmpCmd) = qsubSystem($logDir."MATE.sh",$mateCmd,1,"20G","_MT$JNUM",$jobN,"",1,\@General_Hosts,\%QSBopt) ;
		#### 3 : if mate pairs, process these now (and remove corresponding raw fiiles
	}
	foreach my $hrr (@matePrps){
		my %nateFiles = %{$hrr};
		push(@pa1,$nateFiles{pe1});push(@pa2,$nateFiles{pe2});push(@pas,$nateFiles{se}); push(@libInfo,"pe4mt$ii");
		push(@pa1,$nateFiles{mp1});push(@pa2,$nateFiles{mp2});push(@libInfo,"mate$ii");
		push(@pa1,$nateFiles{un1});push(@pa2,$nateFiles{un2});push(@libInfo,"mate_unkn$ii");
	}
#die "@pa1\n@libInfo\n";
	
	#if (!$calcUnzp ){ return ($jobN,\@pa1,\@pa2, \@pas, $WT, $rawReads, $mmpu,\@libInfo, $jobNUZ);}

#mmpu number of EMBL
	return ($jobN,\@pa1,\@pa2, \@pas, $WT, $rawReads, $mmpu,\@libInfo, $jobNUZ);
}


sub prepKraken($){
	#my ($DBdir) = @_;
	my $oriKrakDir = "/g/scb/bork/hildebra/DB/kraken/";
	my %DBname;
	$DBname{"hum1stTry"} = 1 if ($humanFilter);
	$DBname{"minikraken_2015"} = 1 if ($DO_EUK_GENE_PRED);
	$DBname{$globalKraTaxkDB} = 1 if ($globalKraTaxkDB ne "");
	my $cmd = "";my $jobN= "";
	#die "krakper\n";
	foreach my $kk (keys %DBname ){
		if (!-d "$oriKrakDir$kk"){die "can't find kraken db $oriKrakDir$kk\n";}
		if (!-d "$krakenDBDirGlobal/$kk" && !-e "$krakenDBDirGlobal/$kk/cpFin.stone" ){
			$cmd =  "cp -r $oriKrakDir$kk $krakenDBDirGlobal/\n";
			$cmd .= "touch $krakenDBDirGlobal/$kk/cpFin.stone\n";
		}
	}
	my $tmpCmd="";
	if ($cmd ne ""){
		$jobN= "_KRDB";
		($jobN, $tmpCmd) = qsubSystem("$baseOut/LOGandSUB/krkDB.sh",$cmd,1,"1G",$jobN,"","",1,\@General_Hosts,\%QSBopt) ;
		print "Fix KRDB 2416\n";
	}
	#die $cmd;
	return $jobN;
}

sub krakHSap($ $ $ $ $){
	my ($ar1, $ar2, $ars,$tmpD,$jDep) = @_;
	my @pa1 = @{$ar1}; my @pa2 = @{$ar2}; my @pas = @{$ars};
	return $jDep unless ($humanFilter);
	#my ($DBdir,$DBname) = @_;
	my $DBdir = $krakenDBDirGlobal;
	my $unsplBin = "perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/unsplit_krak.pl ";
	my $DBname = "hum1stTry"; my $numThr = 6;
	my $cmd = "\n\nmkdir -p $tmpD\n\n"; my $tmpF ="$tmpD/krak.tmp.fq";
	my $fileDir = "";
	for (my $i=0;$i<@pa1;$i++){
		my $r1 = $pa1[$i]; my $r2 = $pa2[$i];
		$pa1[$i] =~ m/(.*\/)[^\/]+$/; $fileDir= $1 if ($fileDir eq "");
		$cmd .= "$krkBin --paired --preload --threads $numThr --fastq-input --unclassified-out $tmpF --db $DBdir/$DBname  $r1 $r2 > /dev/null\n";
		#overwrites input files
		$cmd .= "$unsplBin $tmpF $r1 $r2\n";
	}
	$cmd .= "\n\nrm -f $fileDir/krak.stone\n\n";
	for (my $i=0;$i<@pas;$i++){
		my $rs = $pas[$i]; 
		$cmd .= "$krkBin --preload --threads $numThr --fastq-input --unclassified-out $tmpF --db $DBdir/$DBname  $rs > /dev/null\n";
		$cmd .= "rm -f $rs; mv $tmpF $rs\n";
		#overwrites input files
	}
	$cmd .= "\n\n";
	$cmd .= "touch $fileDir/krak.stone\n";
	my $jobN = ""; my $tmpCmd="";
	unless (-e "$fileDir/krak.stone"){
		$jobN = "_KR$JNUM";
		($jobN, $tmpCmd) = qsubSystem($logDir."KrakHS.sh",$cmd,$numThr,"10G",$jobN,$jDep.";$krakDeps","",1,\@General_Hosts,\%QSBopt) ;
	}
	return $jobN;
}

sub krakenTaxEst(){
	my ($arp1,$arp2, $ars, $outD, $tmpD,$name,$jobd) = @_;
	my @pa1 = @{$arp1}; my @pa2 = @{$arp2}; my @pas = @{$ars};
	#$outD.= "$globalKraTaxkDB/";
	my $krakStone = "$outD/krakDone.sto";
	return $jobd if (-d $outD && -e $krakStone);
	my $numCore = $krakenCores;
	my @thrs = (0.01,0.02,0.04,0.06,0.1,0.2,0.3);
	my $cmd = "mkdir -p $outD\nrm -rf $tmpD\nmkdir -p $tmpD\n";
	my $it =0;
	my $curDB = "$krakenDBDirGlobal/$globalKraTaxkDB";
	#die $curDB."\n";
	#paired read tax assign
	for (my $i=0;$i<@pa1;$i++){
		my $r1 = $pa1[$i]; my $r2 = $pa2[$i];
		$pa1[$i] =~ m/(.*\/)[^\/]+$/; 
		$cmd .= "$krkBin --paired --preload --threads $numCore --fastq-input  --db $curDB  $r1 $r2 >$tmpD/rawKrak.$it.out\n";
		for (my $j=0;$j< @thrs;$j++){
			$cmd .= "$krkBin-filter --db $curDB  --threshold $thrs[$j] $tmpD/rawKrak.$it.out | $krkBin-translate --mpa-format --db $curDB > $tmpD/krak_$thrs[$j]"."_$it.out\n";
		}
		#$cmd .= " $krakCnts1 $tmpD/krak$it.out $tmpD/krak$it.tax\n";
		$it++;
	}
	$cmd .= "\n\n";
	#single read tax assign
	for (my $i=0;$i<@pas;$i++){
		my $rs = $pas[$i]; 
		$cmd .= "$krkBin --preload --threads $numCore --fastq-input  --db $curDB  $rs >$tmpD/rawKrak.$it.out\n";
		#$cmd .= " | tee ";#$krkBin-filter --db $curDB  --threshold 0.01 | $krkBin-translate --mpa-format --db $curDB > $tmpD/krak$it.out\n";
		for (my $j=0;$j< @thrs;$j++){
			$cmd .= "$krkBin-filter --db $curDB  --threshold $thrs[$j] $tmpD/rawKrak.$it.out | $krkBin-translate --mpa-format --db $curDB > $tmpD/krak_$thrs[$j]"."_$it.out\n";
		}
		#overwrites input files
		$it++;
	}
	
	#TODO: 1: make table; 2: copy to outD
	$cmd .= "\n\n";
	for (my $j=0;$j< @thrs;$j++){
	$cmd .= "cat $tmpD/krak_$thrs[$j]"."_*.out > $tmpD/allkrak$thrs[$j].out\n";
	$cmd .= "$krakCnts1 $tmpD/allkrak$thrs[$j].out $outD/krak.$thrs[$j].cnt.tax\n";
	}

	$cmd .= "touch $krakStone\n";
	$cmd .= "rm -r $tmpD";

	my $jobName = "_KT$JNUM"; my $tmpCmd;
	if (!-d $outD || !-e $krakStone){
		($jobName,$tmpCmd) = qsubSystem($logDir."KrkTax.sh",$cmd,$numCore,"20G",$jobName,$jobd.";$krakDeps","",1,\@General_Hosts,\%QSBopt);
	}
	return $jobName;
}



sub bam2cram($ $ $ $ $ $ $){#save further space: convert the bam to cram
	my ($nxtBAM,$REF,$del,$subm,$doCram, $stone, $numCore) = @_;
	my $ret = "";
	my $nxtCRAM = $nxtBAM;	
	#my $numCore = 4;
	if (!$doCram){return ($ret,"");}
	$nxtCRAM =~ s/\.bam$/\.cram/;
	if (!-e $nxtBAM && -e $nxtCRAM){print "CRAM exists, but no stone set\n"; system "rm $nxtCRAM";return $ret;}
	#my $stone = $nxtBAM;	$stone =~ s/\.bam$/\.cram\.sto/;
	$ret.="rm -f $nxtCRAM\n" if (-e $nxtCRAM);
	#$ret.="$smtBin view -T $REF -C -o $nxtCRAM $nxtBAM\n";
	$ret.="$sambambaBin view -T $REF -t $numCore -f cram -o $nxtCRAM $nxtBAM\n";
	$ret.="rm -f $nxtBAM\n" if ($del);
	$ret .= "touch $stone\n";
	#die $ret;
	if ($subm){
		my $jobN = "_CRAM$JNUM"; my $tmpCmd;
		($jobN, $tmpCmd) = qsubSystem($logDir."BAM2CRAMxtra.sh",$ret,$numCore,"5G",$jobN,"","",1,[],\%QSBopt);
	}
	return ($ret,$nxtCRAM);
}

sub jgi_depth_cmd_depr($){
	my ($nxtBAM) = @_;
	my $covCmd = "";
	$covCmd .= "\n/g/bork3/home/hildebra/bin/metabat/./jgi_summarize_bam_contig_depths";
	$covCmd .= " --outputDepth $nxtBAM.jgi.depth.txt  --percentIdentity 97 --pairedContigs $nxtBAM.jgi.pairs.sparse $nxtBAM > $nxtBAM.jgi.cov\n";
	$covCmd .= "\ngzip $nxtBAM.jgi.*\n\n";
	if (-e "$nxtBAM.jgi.cov"){$covCmd="";}
	#$covCmd .= "gzip $nxtBAM.jgi*\n";
	return $covCmd;
}

sub mapReadsToRef{
	my ($dirsHr,$outName, $is2ndMap, $par1,$par2,$Ncore,#hm,m,x,x,x,x,x
		$REF,$jDepe, $unaligned,$immediateSubm,#m,x,x,x
		$doCram,$smpl,$libAR) = @_;#x,x,x
		
			# ($dirsHr,$outName, $is2ndMap, $par1,$par2,$Ncore,#hm,m,x,x,x,x,x
	#	$REF,$jDepe, $unaligned,$immediateSubm,#m,x,x,x,m
	#	$doCram,$smpl,$libAR) = @_;#x,x,x
		#get essential dirs...
	my $mappDir = ${$dirsHr}{glbMapDir};	my $nodeTmp = ${$dirsHr}{nodeTmp}; #m,s
	my $tmpOut = ${$dirsHr}{glbTmp};	my $finalD = ${$dirsHr}{outDir}; #s,m
	my @finalDS = split /,/,$finalD;

	my $bwtIdx = $REF.".bw2";
	my @libsOri = @{$libAR};
	#already exists
	#print "IS${unaligned}SI\n";
	my $baseN = $smpl;#$RNAME."_".$QNAME;
	
	my $bashN = "";	if ($is2ndMap){$bashN = "$outName"; $bashN =~ s/,/./;}
	my @outNms = split /,/,$outName;
	if ($outName ne "" ){
		$baseN = $outNms[0];
	} else {
		@outNms = ($baseN);
	}
	my $tmpOut22 = $tmpOut."/$baseN.iniAlignment.bam";
	#global value overwrites local value
	if ($doCram){$doCram = $doBam2Cram;}
	my @pa1 = @{$par1}; my @pa2 = @{$par2};
	my %params;
	my $bamFresh = 0; #is the bam newly being created?
	my $isSorted = 0;		$isSorted=1 if (@pa1==1 && exists($make2ndMapDecoy{Lib}));
	$params{sortedbam}=$isSorted; $params{bamIsNew} = $bamFresh; $params{is2ndMap} = $is2ndMap;
	$params{immediateSubm} =  $immediateSubm;;
	
	my $outputExists = 1;
	for (my $k=0;$k<@outNms;$k++){
		$tmpOut22 = $tmpOut."/$outNms[$k].iniAlignment.bam"; 
		#print "$tmpOut22\n";
		if ($rewrite2ndMap){$outputExists=0; system "rm -fr $tmpOut22"; next;}
		next if (-e "$finalDS[$k]/$outNms[$k]-smd.bam.coverage.gz" || -e "$mappDir/$outNms[$k]-smd.bam.coverage.gz");
		#print "-e $finalDS[$k]/$outNms[$k]-smd.bam.coverage.gz || -e $mappDir/$outNms[$k]-smd.bam.coverage.gz";
		if (!-e $tmpOut22){ $outputExists=0; last;}
	}
	#die "$outputExists\n";
	if ($outputExists){
		return ("","",\%params);
	} 
	
	
	my $nxtCRAM = "$tmpOut/$baseN-smd.cram";
	my $tmpOut21 = $nodeTmp."/$baseN.iniAlignment.bam";
	my $sortTMP = $nodeTmp."/$baseN.srt";
	#print "$nxtBAM\n";
	#die("too far\n");
	
	my $retCmds="";my $tmpCmd=""; my $xtraSamSteps1="";
	$nodeTmp.="/map/";
	#my $nxtCRAM = "$mappDir/$baseN-smd.cram";

	system("mkdir -p ".join(" ",split(/,/,$mappDir))." " .join(" ",split(/,/,$tmpOut)));
	#move ref DB & unzip
	#system("rm -r $tmpOut\n mkdir -p $tmpOut\n cp $REF $tmpOut");
	#$REF = $tmpOut.basename($REF);
	my $unzipcmd = "";
	if ($REF =~ m/\.gz$/){$unzipcmd .= "gunzip $REF\n";$REF =~ s/\.gz$//;}
	my $jobN = "";
	my $tmpUna = $tmpOut."/unalTMP/";
	#Error: No EOF block on /g/scb/bork/hildebra/SNP/MeHiAss/MH0411//tmp//mapping/Alignment.bam, possibly truncated file.
	if (-e $logDir.$bashN."bwtMap2.sh.etxt"){open I,"<$logDir/".$bashN."bwtMap2.sh.etxt" or die "Can't open old bowtie_2 logfile $logDir\n"; my $str = join("", <I>); close I;
		if ($str =~ /Error: No EOF block on (.*), possibly truncated file\./){system("rm $1 $tmpOut22");}	
		close I;	
	}
	#depending on setup, mappDir == tmpOut
	my $algCmd = "rm -rf ".join(" ",split(/,/,$mappDir))."\nmkdir -p ".join(" ",split(/,/,$mappDir))." $tmpOut $nodeTmp\n"; 
	system ("mkdir -p $logDir/mapStats/");
	my $statsF = $logDir."mapStats/".$bashN."mapStats.txt";
	
	my @regs;#subset of DB seqs to filter for 
	my @reg_lcs;
	$REF =~ m/([^\/]+)$/; my $REFnm = $1; $REFnm =~ s/\.fna//; my $decoyModeActive=0;
	if (exists($make2ndMapDecoy{Lib})){
		$bwtIdx = "$nodeTmp/$REFnm.decoyDB.fna";
		$algCmd.= "\n$decoyDBscr $REF $make2ndMapDecoy{Lib} $bwtIdx $Ncore $outName $finalD\n";
		#die "$algCmd\n";
		$bwtIdx .= ".bw2";
		$decoyModeActive=1;
		$xtraSamSteps1 = "$smtBin sort -@ $Ncore -T $sortTMP - |";
		@regs = @{$make2ndMapDecoy{regions}};
		@reg_lcs =  @{$make2ndMapDecoy{region_lcs}};
		#for (my $k=0;$k<@regs;$k++){
		#die "$bwtIdx\n$algCmd\n";
	}
	
	my @subBams; my @subBamsCnt;
	if ($MapperProg==1 || $MapperProg==2){
		
		#die "TODO bowtie 2 \n";
		my $algCmdBase = "$bwt2Bin --no-unal --end-to-end -p $Ncore -x $bwtIdx ";
		if ($unaligned ne ""){
			#die "IS${unaligned}SI\n";
			#system "mkdir -p $unaligned $tmpUna\n";
			$algCmd .= "mkdir -p $unaligned $tmpUna\n\n";
			$algCmdBase .= " --un-conc $tmpUna ";
		}
		if ($MapperProg==2){
			$algCmdBase = "$bwaBin mem -t $Ncore ";
		}
		#--al-conc #prints pairs that hit
		#my $tmpOut22 = $tmpOut."ali.2.sam";
		
		#join(",",@pa2)
		
		for (my $i=0; $i< @pa1; $i++){
			#$pa1[$i] =~ m/\/([^\/]+)\.f.*q$/;
			if ($MapperProg==1){
				my $rgID = "$smpl"; my $rgStr = "--rg SM:$smpl --rg PL:ILLUMINA "; #PU:lib1
				$rgStr .= "--rg LB:$libsOri[$i] ";
				$rgStr .= " -X $mateInsertLength " if ($libsOri[$i] =~ m/mate/);
				$algCmd .="$algCmdBase -1 ".$pa1[$i]." -2 ".$pa2[$i]." --rg-id $rgID $rgStr ";
			} elsif ($MapperProg==2){
				my $rgID = "$smpl"; my $rgStr = "'\@RG\\tID:$1\\tSM:$smpl\\tPL:ILLUMINA";
				$rgStr .= "\\tLB:lib1'";
				$algCmd .= $algCmdBase." -R $rgStr $REF " .$pa1[$i]." ".$pa2[$i];#| $smtBin view -@ $Ncore -b - > $tmpOut21.$i\n  $smtBin flagstat $tmpOut21.$i > $statsF.$i \n";
			}
			$algCmd .= "| $xtraSamSteps1  $smtBin view -b1 -@ $Ncore -F 4 - > $tmpOut21.$i\n";
			
			if ($decoyModeActive){#remove unnecessary reads
				$algCmd .= "$smtBin index $tmpOut21.$i\n";
				for (my $k=0;$k<@regs;$k++){
					$algCmd .= "$bamHdFilt_scr $tmpOut21.$i $reg_lcs[$k] 0 > $tmpOut21.$i.decoy.sam.$k\n";
					$algCmd.= "  $smtBin view -@ $Ncore $tmpOut21.$i $regs[$k] >> $tmpOut21.$i.decoy.sam.$k\n";
					$algCmd.= "  $smtBin view -b -h -@ $Ncore $tmpOut21.$i.decoy.sam.$k > $tmpOut21.$i.decoy.$k\nrm $tmpOut21.$i.decoy.sam.$k \n";
					$subBams[$k] .= " $tmpOut21.$i.decoy.$k";  
				}
				
				$algCmd.= "rm $tmpOut21.$i\n";# mv $tmpOut21.decoy.$i $tmpOut21.$i\n";
			} else {
				$subBams[0] .= " $tmpOut21.$i"; 
			}

			#| $smtBin view -1 - |tee > $tmpOut21.$i  |  $smtBin flagstat - > $statsF.$i \n\n";
			
		}
		#" -S $tmpOut21\n";
		#$algCmd .= "rm -f $bwtIdx"."*\n";
		for (my $k=0;$k<@subBams;$k++){
			$tmpOut22 = $tmpOut."/$outNms[$k].iniAlignment.bam";
#			print $tmpOut22."\n";
			#if (-e $tmpOut22){ next;}
			$bamFresh=1;
			if (@pa1 > 1){
				$algCmd .= "\n$smtBin cat ".$subBams[$k]." > $tmpOut22\n";
				$algCmd .= "\nrm ".$subBams[$k]."\n";
			} else { #nothing else, compression is now in the bowtie2 phase
				#$algCmd .= "\n$smtBin view -bS -F 4 -@ $Ncore $subBams[0] > $tmpOut22\n";
				$algCmd .= "\nmv $subBams[$k] $tmpOut22\n";
			}
		}
		#die "$algCmd\n";
		my $mapProgNm = "";
		if ($MapperProg==1){
			$mapProgNm = "bowtie2";
		}elsif ($MapperProg==2){
			$mapProgNm = "bwa";
		}
		foreach my $mpdSS (split(/,/,$mappDir)) {
			$algCmd .= "echo '$mapProgNm' > $mpdSS/map.sto\n"; last;
		}

		#system($algCmd." -U $p2 -S $tmpOut22\n");
		#die();
	}
	my $nodeCln = "\nrm -rf $nodeTmp\n";
	my $unalignCmd = "";
	if ($unaligned ne "" && $MapperProg == 1 && !-e "$unaligned/unal.1.fq.gz"){
		$unalignCmd .= "mkdir -p $unaligned\n";
		$unalignCmd .= "$pigzBin -p $Ncore -c $tmpUna/un-conc-mate.1 > $unaligned/unal.1.fq.gz;\n$pigzBin -p $Ncore -c $tmpUna/un-conc-mate.2 > $unaligned/unal.2.fq.gz\n";
		$unalignCmd .= "rm -f $tmpUna/un-conc-mate.*\n";
	}
	#if (exists($make2ndMapDecoy{Lib})){#rm newly created decoy lib again..
	#	$bwtIdx =~ s/.bw2$//;
	#	$nodeCln .= "rm $bwtIdx*\n";
	#}
	#die "$unzipcmd.$algCmd.$nodeCln";
	if ( $bamFresh   || !-e "$mappDir/map.sto"){
		$jobN = "_BT$JNUM$outNms[0]"; $bamFresh = 1; 
		($jobN,$tmpCmd) = qsubSystem($logDir.$bashN."bwtMap.sh",
				$unzipcmd.$algCmd.$unalignCmd.$nodeCln,
				$Ncore,$bwtMapMem,$jobN,$jDepe,"",$immediateSubm,\@General_Hosts,\%QSBopt) ;
		$retCmds .= $tmpCmd;
	}
	#xtraSamSteps1 = isSorted
	$params{sortedbam}=$isSorted; $params{bamIsNew} = $bamFresh; $params{is2ndMap} = $is2ndMap;
	$params{immediateSubm} = $immediateSubm;
	#die "$jobN\n$immediateSubm\n";
	return($jobN,$retCmds,\%params);
}	
	
sub bamDepth{
	my ($dirsHr,  #   $mappDir,$tmpOut,$nodeTmp,$finalD,#dirs to work in
			$outName,$REF,$jDep,$doCram,$mapparhr) = @_;
	my %params = %{$mapparhr};
	my ($isSorted , $bamFresh,$is2ndMap) =($params{sortedbam},$params{bamIsNew},$params{is2ndMap});
	#recreate base pars from map2tar sub
	
	my $immediateSubm = $params{immediateSubm} ;
	my $mappDir = ${$dirsHr}{glbMapDir};	my $nodeTmp = ${$dirsHr}{nodeTmp};
	my $tmpOut = ${$dirsHr}{glbTmp};	my $finalD = ${$dirsHr}{outDir};

	my $baseN = $outName;
	my $bashN = "";	if ($is2ndMap){$bashN = "$outName"; $bashN =~ s/,/./;}
	my $tmpOut22 = $tmpOut."/$baseN.iniAlignment.bam";
	my $cramSTO = "$mappDir/$baseN-smd.cram.sto";
	my $nxtBAM = "$tmpOut/$baseN-smd.bam";
	my $sortTMP = $nodeTmp."/$baseN.srt";
	
	#check if already done
	if ($doCram && (-e "$finalD/$baseN-smd.cram.sto" && -s "$finalD/$baseN-smd.bam.coverage.gz" && !$params{bamIsNew} ) ){
		return ("","");
	} elsif (!$doCram && -e "$finalD/$baseN-smd.bam" && -s "$finalD/$baseN-smd.bam.coverage.gz" && !$params{bamIsNew}){
		return ("","");
	}
	if ( -e "$mappDir/$baseN-smd.bam.coverage.gz" && (-e "$mappDir/$baseN-smd.bam" || -e "$mappDir/$baseN-smd.cram") && !$params{bamIsNew}){ #already stored in mapping dir, still needs to be copied
		return ("","1");
	}
	#die "assd";
	
	
	
	my $numCore = 6;#new functionality with sambamba
	#my $biobambamBin = "bamsormadup";
	#depth profile
	#my $cmd = "$smtBin faidx $REF\n $smtBin view -btS $REF.fai $tmpOut/$baseN.sam > $tmpOut/$baseN.bam\n";
	my $cmd = "sleep 3\n";
	#$cmd .= "mkdir -p $nodeTmp\n";
	$cmd .= "mkdir -p $mappDir $tmpOut $nodeTmp\n";
	#sort & remove duplicates
	if (1 || -z $nxtBAM){
		#$cmd .= "$sambambaBin sort -t $numCore -o $tmpOut2s $tmpOut22 \n$sambambaBin markdup -r -t $numCore $tmpOut2s $nxtBAM\n";
		#$cmd .= "$biobambamBin sort -t $numCore $tmpOut22 /dev/stdout | $sambambaBin markdup -r -t $numCore /dev/stdin $nxtBAM\n";
		#$cmd .= "$novosrtBin --removeDuplicates --ram 100G -o $nxtBAM -i $tmpOut22 \n\n";
		#$cmd .= "/g/bork3/home/hildebra/bin/htsbox/./htsbox pileup -f $REF -Q20 -q30 -Fs3 sorted.bam > $REF.cns\n" if ($map_DoConsensus);
		
		if ($doRmDup){
			if ($isSorted  ){#no sort required
				$cmd .= "$smtBin rmdup -s $tmpOut22 $nxtBAM;\n";
			} else {
				$cmd .= "$smtBin sort -@ $numCore -T $sortTMP $tmpOut22 | $smtBin rmdup -s - $nxtBAM;\n";
			}
			#$cmd .= "$smtBin sort -@ $numCore $tmpOut22  | $sambambaBin markdup -r -t $numCore /dev/stdin $nxtBAM\n";
		} else {
			if ($isSorted ){
				$cmd .= "mv $tmpOut22 $nxtBAM\n";
			} else {
				$cmd .= "$smtBin sort -@ $numCore -T $sortTMP $tmpOut22  > $nxtBAM\n";
			}
		}
		$cmd .= "$smtBin index $nxtBAM\n";
		$cmd .= "rm -f $tmpOut22 \n" if (!$isSorted);
	}
	#$cmd .= "$smtBin view -bS $tmpOut21 > $tmpOut/$baseN.bam\n";
	#my $mdJar = "/g/bork5/hildebra/bin/picard-tools-1.119/MarkDuplicates.jar";
	#$cmd .= "java -Xms1g -Xmx24g -XX:ParallelGCThreads=2 -XX:MaxPermSize=1g -XX:+CMSClassUnloadingEnabled ";		$cmd .= "-jar $mdJar INPUT=$mappDir/$baseN-s.bam OUTPUT=$mappDir/$baseN-smd.bam METRICS_FILE=$mappDir/$baseN-smd.metrics ";		
	#$cmd .= "AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE\n\n";
	
	#$cmd .= "$smtBin index $nxtBAM\n";#$mappDir/$baseN-smd.bam\n";
	#takes too much space, calc on the fly!
	#$cmd .= "$smtBin depth $mappDir/$baseN-s.bam > $mappDir/$baseN.depth.bam.txt\n";#depth file (different estimator)
	my $covCmd = "";
	if (!-e "$nxtBAM.coverage"){
		$covCmd = "/g/bork5/hildebra/bin/bedtools2-2.21.0/bin/genomeCoverageBed -ibam $nxtBAM -bg > $nxtBAM.coverage\n";
		$covCmd .= "rm -f $nxtBAM.coverage.gz\ngzip $nxtBAM.coverage\n";
	}
	
	#jgi depth profile - not required, different call to metabat
	$covCmd .= jgi_depth_cmd($nxtBAM,$nxtBAM,95);
	#$covCmd .= "/g/bork5/hildebra/bin/bedtools2-2.21.0/bin/genomeCoverageBed -ibam $nxtBAM -bg  | ".'awk \'BEGIN {pc=""} {	c=$1;	if (c == pc) {		cov=cov+$2*$5;	} else {		print pc,cov;		cov=$2*$5;	pc=c}';
	#$covCmd .= "} END {print pc,cov}\' $nxtBAM.coverage | tail -n +2 > $nxtBAM.coverage.percontig";
	
	my ($CRAMcmd,$CRAMf) = bam2cram($nxtBAM,$REF,1,0,$doCram,$cramSTO, $numCore);
	$CRAMcmd .= "\nmv $nxtBAM* $CRAMf $mappDir\n" unless ($mappDir eq $tmpOut);
	$CRAMcmd .= "echo \"".basename($nxtBAM)."\" > $mappDir/done.sto\n";


	my $newJobN = "_SBAM$JNUM$outName";	
	#subsequent jobs not dependent on this one
	my $jobN2 = $jDep; my $retCmds="";
	$covCmd = "" if ( -s "$nxtBAM.coverage");
	
	#my $cleaner = "mv $tmpD $fin
	if (0){
		print "B1 " if ( $doCram && !-e $cramSTO); print "B2 " if (  !-s "$nxtBAM.coverage"); print "B3 " if ( $bamFresh);
		print " ".$cramSTO."\n";
	}
	my $nodeCln = "\nrm -rf $nodeTmp\n";
	#die "($doCram && !-e $cramSTO) || ( !-s $nxtBAM.coverage) || $bamFresh\n";
	if ( ($doCram && !-e $cramSTO) || ( !-s "$nxtBAM.coverage") || $bamFresh){
		($jobN2,$retCmds) = qsubSystem($logDir.$bashN."bwtMap2.sh",
				$cmd."\n".$covCmd."\n".$CRAMcmd."\n$nodeCln\n"#.$covCmd2
				,$numCore,int(110/$numCore)."G",$newJobN,$jDep,"",$immediateSubm,\@General_Hosts,\%QSBopt);
	}

	#die();
	return($jobN2,$retCmds);
}

sub metphlanMapping{
	my ($inF1a,$inF2a,$tmpD,$finOutD,$smp,$Ncore,$deps) = @_;
	my $stone = $finOutD."$smp.MP2.sto";
	if (!$DoMetaPhlan ||  -e $stone){return;}
	my @car1 = @{$inF1a}; my @car2 = @{$inF2a};
	my $inF1 = join(",",@car1); my $inF2 = join(",",@car2); 
	system "mkdir -p $finOutD" unless (-d $finOutD);
	my $finOut = $finOutD."$smp.MP2.txt";
	my $finOut_noV = $finOutD."$smp.MP2.noV.txt";
	my $finOut_noVB = $finOutD."$smp.MP2.noV.noB.txt";
	my $finOut_Vo = $finOutD."$smp.MP2.VirusOnly.txt";
	my $sam = "$tmpD/metph2.sam";
	#path to metaphlan DB
	my $mpDB = "/g/bork3/home/hildebra/bin/metaphlan2/db_v20/mpa_v20_m200";
	my $cmd = "mkdir -p $tmpD\n$bwt2Bin --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S $sam -p $Ncore -x $mpDB -1 $inF1 -2 $inF2 \n";
	$cmd .= "$metPhl2Bin $sam --input_type sam --mpa_pkl $mpDB.pkl > $finOut\n";
	$cmd .= "$metPhl2Bin $sam --input_type sam --ignore_viruses --mpa_pkl $mpDB.pkl > $finOut_noV\n";
	$cmd .= "$metPhl2Bin $sam --input_type sam --ignore_bacteria --ignore_viruses --ignore_archaea --mpa_pkl $mpDB.pkl > $finOut_noVB\n";
	$cmd .= "$metPhl2Bin $sam --input_type sam --ignore_bacteria --ignore_eukaryotes --ignore_archaea --mpa_pkl $mpDB.pkl > $finOut_Vo\n";
	$cmd .= "rm $sam\n";
	my $mergeStr = "$metPhl2Merge *.MP2.txt > $finOutD/comb.MP2.txt";
	$cmd .= "echo \' $mergeStr \' > $stone\n";
	my $jobN = "MP2$JNUM";

	my ($jobN2,$tmpCmd) = qsubSystem($logDir."metaPhl2.sh",
			$cmd,$Ncore,"3G",$jobN,$deps,"",1,[],\%QSBopt);
	$jobN  = $jobN2;
	return $jobN;
}


sub clean_tmp{#routine moves output from temp dirs to final dirs (that are IO limited, thus single process)
	my ($clDar,$cpref,$cpnodel,$jDepe,$tag,$curJname) = @_;
	my @clDa = @{$clDar};
	my @cps = @{$cpref};
	my @cpsND = @{$cpnodel};
	my $cmd = "";
	for (my $i=0;$i<@cps;$i+=2){
		next if (length($cps[$i+1] ) < 5);
		$cmd .= "rm -r -f $cps[$i+1]\nmkdir -p $cps[$i+1]\nrsync -r --remove-source-files $cps[$i] $cps[$i+1]\n" if (@cps > 1);
	}
	for (my $i=0;$i<@cpsND;$i+=2){
		next if (length($cpsND[$i+1] ) < 5);
		$cmd .= "mkdir -p $cpsND[$i+1]\nrsync -r --remove-source-files $cpsND[$i] $cpsND[$i+1]\n"  if (@cpsND > 1);
	}
	$cmd .= "rm -f -r ".join(" ",@clDa)."\n" if (@clDa > 0  );
	$cmd .= "echo $tag >> $collectFinished\n" unless ($tag eq "");
	my $clnBash = $logDir."clean$curJname.sh";
	if ($curJname eq ""){
		$curJname = "_cln$JNUM";
	}
	#add extra job dependency on d2star
	my $xtraJDep = "";
	#$xtraJDep = $QSBoptHR->{rTag}."_d2met";
	my ($jdep,$tmpCmd) = qsubSystem($clnBash,$cmd,1,"1G",$curJname,$jDepe.";$xtraJDep","",1,[],\%QSBopt);
	return ($jdep );
}



sub remComma($){
	my ($in) = @_;
	$in =~ s/,//g;
	return $in;
}
sub smplStats(){
	my ($inD,$assDir) = @_;
	#print "STATS   \n";
	my %ret; my $outStr = ""; my $outStrDesc = "";my $outStr5 = ""; my $outStrDesc5 = "";
	my $filStats = "";
		#return ($outStrDesc,$outStr,$outStrDesc5, $outStr5);
	
	if (-e "$inD/LOGandSUB/KrakHS.sh.etxt"){
		$filStats = `cat $inD/LOGandSUB/KrakHS.sh.etxt`;
		if ($filStats =~ m/\d+ sequences classified \((\d+\.?\d*)%\)/){
			$outStr.= "$1%\t";
		} else {
			$outStr.= "?%\t";
		}
	} else {$outStr .= "\t" x 1;}
	$outStrDesc .= "HumanRdsPerc\t";
	
		#seq filter stats  /LOGandSUB/sdm/filter.log
	$filStats = "";
	$filStats = `cat $inD/LOGandSUB/sdm/filter.log` if (-e "$inD/LOGandSUB/sdm/filter.log");
	
	if (!$readsRpairs && $filStats =~ m/Reads processed: ([0-9,]+)/){#single end format
		my $totRds =  remComma($1);
		$filStats =~ m/Rejected: ([0-9,]+)\n/;
		#die $1. "   $2\n";
		my $Rejected1 =  remComma($1);my $Rejected2 =  0;
		$filStats =~ m/Accepted: ([0-9,]+) \(\d+/;
		my $Accepted1 =  remComma($1);my $Accepted2 =  0;
		my $Singl1 =  0;my $Singl2 =  0;
		$filStats =~ m/- Seq Length :\s*\d+.*\/(\d+.*)\/\d+.*\n/;
		my $AvgLen = $1;
		$filStats =~ m/- Quality :\s*\d+.*\/(\d+.*)\/\d+.*\n/;
		my $AvgQual = $1;
		$filStats =~ m/- Accum. Error (\d+\.\d+)\n/;
		my $accErr = $1;
		$outStr .= "$totRds\t$Rejected1\t$Rejected2\t$Accepted1\t$Accepted2\t$Singl1\t$Singl2\t$AvgLen\t$AvgQual\t$accErr\t";
	}elsif ($readsRpairs && $filStats =~ m/Reads processed: ([0-9,]+); ([0-9,]+) \(pa/)
			{ #only do this if newest format
		my $totRds =  remComma($1);
		$filStats =~ m/Rejected: ([0-9,]+); ([0-9,]+)\n/;
		#die $1. "   $2\n";
		my $Rejected1 =  remComma($1);my $Rejected2 =  remComma($2);
		$filStats =~ m/Accepted: ([0-9,]+); ([0-9,]+) \(\d+/;
		my $Accepted1 =  remComma($1);my $Accepted2 =  remComma($2);
		#Singletons among these: 269,516; 6,686
		$filStats =~ m/Singletons among these: ([0-9,]+); ([0-9,]+)\n/;
		my $Singl1 =  remComma($1);my $Singl2 =  remComma($2);
		
		$filStats =~ m/- Seq Length :\s*\d+.*\/(\d+.*)\/\d+.*\n/;
		my $AvgLen = $1;
		$filStats =~ m/- Quality :\s*\d+.*\/(\d+.*)\/\d+.*\n/;
		my $AvgQual = $1;
		$filStats =~ m/- Accum. Error (\d+\.\d+)\n/;
		my $accErr = $1;
		$outStr .= "$totRds\t$Rejected1\t$Rejected2\t$Accepted1\t$Accepted2\t$Singl1\t$Singl2\t$AvgLen\t$AvgQual\t$accErr\t";
		
	} else {
		$outStr .= "\t" x 10;
	}
	$outStrDesc .= "totRds\tRejected1\tRejected2\tAccepted1\tAccepted2\tSingl1\tSingl2\tAvgLen\tAvgQual\taccErr\t";
	#check for flash merged reads
	if (-e "$inD/LOGandSUB/flashMrg.sh.otxt"){
		$filStats = `cat $inD/LOGandSUB/flashMrg.sh.otxt`;
		if ($filStats =~ m/\[FLASH\]     Combined pairs:   (\d+)/){$outStr.="$1\t";} else {$outStr.="?\t";}
		if ($filStats =~ m/\[FLASH\]     Uncombined pairs: (\d+)/){$outStr.="$1\t";} else {$outStr.="?\t";}
	} else {
		$outStr .= "\t" x 2;
	}
	$outStrDesc .= "Merged\tNotMerged\t";
	
	
	#check if corrected dir still exists..
	system "rm -rf $inD/assemblies/metag/corrected" if (-d "$inD/assemblies/metag/corrected");
	#assembly stats
	my $tmpassD = "";
	if (-e "$inD/assemblies/metag/assembly.txt"){
		$tmpassD = `cat $inD/assemblies/metag/assembly.txt`; chomp $tmpassD;
	} elsif (-e "$inD/assemblies/metag/AssemblyStats.txt"){
		$tmpassD = "$inD/assemblies/metag/";
	} elsif (-e $assDir."metag/AssemblyStats.txt"){
		$tmpassD = "$assDir/metag/";
	}
#print $tmpassD."\n";
	if (-e "$tmpassD/AssemblyStats.txt"){
		my $assStats = `cat $tmpassD/AssemblyStats.txt`;
		if ($assStats =~ m/Number of scaffolds\s+(\d+)/){ $ret{NScaff} = $1;} else { die "scfcnt wrg1\n";}
		if ($assStats =~ m/Total size of scaffolds\s+(\d+)/){ $ret{ScaffSize} = $1;} else { die "scfcnt wrg2\n";}
		if ($assStats =~ m/Longest scaffold\s+(\d+)/){ $ret{ScaffMaxSize} = $1;} else { die "scfcnt wrg3\n";}
		if ($assStats =~ m/N50 scaffold length\s+(\d+)/){ $ret{ScaffN50} = $1;} else { die "scfcnt wrg4\n";}
		if ($assStats =~ m/Number of scaffolds > 1K nt\s+(\d+)/){ $ret{NScaffG1k} = $1;} else { die "scfcnt wrg5\n";}
		if ($assStats =~ m/Number of scaffolds > 10K nt\s+(\d+)/){ $ret{NScaffG10k} = $1;} else { die "scfcnt wrg6\n";}
		if ($assStats =~ m/Number of scaffolds > 100K nt\s+(\d+)/){ $ret{NScaffG100k} = $1;} else { die "scfcnt wrg7\n";}
		$outStr .= "$ret{NScaff}\t$ret{NScaffG1k}\t$ret{NScaffG10k}\t$ret{NScaffG100k}\t$ret{ScaffN50}\t$ret{ScaffMaxSize}\t$ret{ScaffSize}\t";
	} else {
		$outStr .= "\t" x 7;
	}
	$outStrDesc .="NScaff400\tNScaffG1k\tNScaffG10k\tNScaffG100k\tScaffN50\tScaffMaxSize\tScaffSize\t";
	
		#assembly stats
	my $assStats="";
	if (-s "$tmpassD/AssemblyStats.500.txt"){
		my $assStats = `cat $tmpassD/AssemblyStats.500.txt`;
	}
	if (length($assStats) > 100){
		if ($assStats =~ m/Number of scaffolds\s+(\d+)/){ $ret{NScaff} = $1;} else { die "scfcnt wrg1\n";}
		if ($assStats =~ m/Total size of scaffolds\s+(\d+)/){ $ret{ScaffSize} = $1;} else { die "scfcnt wrg2\n";}
		if ($assStats =~ m/Longest scaffold\s+(\d+)/){ $ret{ScaffMaxSize} = $1;} else { die "scfcnt wrg3\n";}
		if ($assStats =~ m/N50 scaffold length\s+(\d+)/){ $ret{ScaffN50} = $1;} else { die "scfcnt wrg4\n";}
		if ($assStats =~ m/Number of scaffolds > 1K nt\s+(\d+)/){ $ret{NScaffG1k} = $1;} else { die "scfcnt wrg5\n";}
		if ($assStats =~ m/Number of scaffolds > 10K nt\s+(\d+)/){ $ret{NScaffG10k} = $1;} else { die "scfcnt wrg6\n";}
		if ($assStats =~ m/Number of scaffolds > 100K nt\s+(\d+)/){ $ret{NScaffG100k} = $1;} else { die "scfcnt wrg7\n";}
		$outStr5 .= "$ret{NScaff}\t$ret{NScaffG1k}\t$ret{NScaffG10k}\t$ret{NScaffG100k}\t$ret{ScaffN50}\t$ret{ScaffMaxSize}\t$ret{ScaffSize}\t"; 
	} else {
		$outStr5 .= "\t" x 7;
	}
	$outStrDesc5 .="NScaff500\tNScaffG1k\tNScaffG10k\tNScaffG100k\tScaffN50\tScaffMaxSize\tScaffSize\t";

	#bowtie map
	my @spl = ();
	if (-s "$inD/LOGandSUB/bwtMap.sh.etxt"){
		my $alignStats = `cat $inD/LOGandSUB/bwtMap.sh.etxt`;
		@spl = split(/\n/,$alignStats);
	}
	my $idx =0;my $dobwtStat=0;
	if (@spl > 12){
		$dobwtStat=1;
		while($spl[$idx] !~ m/\d+ reads; of these:/){
			$idx++; 
			last if ($idx > @spl);
		}
		if ($spl[$idx+0] =~ m/(\d+) reads; of these:/){
			$ret{totReadPairs} = $1;
		} else {
			$dobwtStat=0;#die "wrong bwtOut: $inD/LOGandSUB/bwtMap.sh.etxt X $spl[2]";
		}
	}
	if ($dobwtStat){
		if ($spl[$idx+2] =~ m/(\d+) \(.+\) aligned concordantly 0 times/){ $ret{notAlign} = $1;} else { die "bwtOut wrg1\n";}
		if ($spl[$idx+3] =~ m/(\d+) \(.+\) aligned concordantly exactly 1 time/ ){ $ret{uniqAlign} = $1;} else { die "bwtOut wrg2\n";}
		if ($spl[$idx+4] =~ m/(\d+) \(.+\) aligned concordantly >1 times/ ){ $ret{multAlign} = $1;} else { die "bwtOut wrg3  $spl[6]\n";}
		if ($spl[$idx+7] =~ m/(\d+) \(.+\) aligned discordantly 1 time/ ){ $ret{DisconcAlign} = $1;} else { die "bwtOut wrg4 $spl[9]\n";}
		if ($spl[$idx+12] =~ m/(\d+) \(.+\) aligned exactly 1 time/ ){ $ret{SinglAlign} = $1;} else { die "bwtOut wrg5 $spl[14]\n";}
		if ($spl[$idx+13] =~ m/(\d+) \(.+\) aligned >1 times/ ){ $ret{SinglAlignMult} = $1;} else { die "bwtOut wrg6 $spl[15]\n";}
		if ($spl[$idx+14] =~ m/(\d+\.\d+)\% overall alignment rate/ ){ $ret{AlignmRate} = $1;} else { die "bwtOut wrg7 $spl[16]\n";}
		$outStr .= "$ret{totReadPairs}\t$ret{AlignmRate}\t".$ret{uniqAlign}/$ret{totReadPairs}*100 ."\t".$ret{multAlign}/$ret{totReadPairs}*100 ."\t".$ret{DisconcAlign}/$ret{totReadPairs}*100 
			."\t".$ret{SinglAlign}/$ret{totReadPairs}*100 ."\t".$ret{SinglAlignMult}/$ret{totReadPairs}*100 ."\t";
		$outStr5 .= "$ret{totReadPairs}\t$ret{AlignmRate}\t".$ret{uniqAlign}/$ret{totReadPairs}*100 ."\t".$ret{multAlign}/$ret{totReadPairs}*100 ."\t".$ret{DisconcAlign}/$ret{totReadPairs}*100 
			."\t".$ret{SinglAlign}/$ret{totReadPairs}*100 ."\t".$ret{SinglAlignMult}/$ret{totReadPairs}*100 ."\t";
	} else {
		$outStr .= "\t" x 7;
		$outStr5 .= "\t" x 7;
	}
	$outStrDesc .="ReadsPaired\tOverallAlignment\tUniqueAlgned\tMultAlign\tDisconcAlign\tSingleUniqAlign\tSingleMultiAlign\t";
	$outStrDesc5 .="ReadsPaired\tOverallAlignment\tUniqueAlgned\tMultAlign\tDisconcAlign\tSingleUniqAlign\tSingleMultiAlign\t";
	
	#novocraft sort
	my $doDup = 1;
	if (!-e "$inD/LOGandSUB/bwtMap2.sh.etxt"){$doDup=0;}
	if ($doDup){
		my $alignStats2 = `cat $inD/LOGandSUB/bwtMap2.sh.etxt`;
		if ($alignStats2 =~ m/Proper Pairs\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/){$ret{EstLibSize} = $4; $ret{duplOptic} = $3; $ret{duplPCR} = $2; $ret{duplPass} = $1;} else {$doDup=0;}#die "Can't find dupl stats: $alignStats2\n";}
		if ($alignStats2 =~ m/Improper Pairs\s+(\d+)\s+(\d+)\s+(\d+)/){$ret{duplOptic} += $3; $ret{duplPCR} += $2; $ret{duplPass} += $1;} else {$doDup=0;}#{die "Can't find impr. dupl stats: $alignStats2\n";}
	}
	if ($doDup){
		$outStr .= "$ret{duplOptic}\t$ret{duplPCR}\t$ret{duplPass}\t$ret{EstLibSize}\t";
		$outStr5 .= "$ret{duplOptic}\t$ret{duplPCR}\t$ret{duplPass}\t$ret{EstLibSize}\t";
	} else {
		$outStr .= "\t" x 3;
		$outStr5 .= "\t" x 3;
	}
	$outStrDesc .= "OpticalDuplicates\tPCRduplicates\tPassedMD\tBS_estLibSize\t";
	$outStrDesc5 .= "OpticalDuplicates\tPCRduplicates\tPassedMD\tBS_estLibSize\t";
	#die "$inD/assemblies/metag/ContigStats/GeneStats.txt";
if (-e "$inD/assemblies/metag/ContigStats/GeneStats.txt"){	
		open I,"<$inD/assemblies/metag/ContigStats/GeneStats.txt"; 
		#GeneNumber	AvgGeneLength	AvgComplGeneLength	BpGenes	BpNotGenesGcomplete	G5pComplete	G3pComplete	Gincomplete
		my $tmp = <I>;$outStrDesc .= $tmp;$outStrDesc5 .= $tmp;
		$tmp = <I>; $outStr .= $tmp; $outStr5 .= $tmp; 
		close I;
	} else {
		$outStrDesc .= "GeneNumber	AvgGeneLength	AvgComplGeneLength	BpGenes	BpNotGenesGcomplete	Gcomplete	G5pComplete	G3pComplete	Gincomplete";
		$outStrDesc5 .= "GeneNumber	AvgGeneLength	AvgComplGeneLength	BpGenes	BpNotGenesGcomplete	Gcomplete	G5pComplete	G3pComplete	Gincomplete";
		$outStr .= "\t" x 7;$outStr .= "\t" x 7;
	}
	chomp $outStrDesc; chomp $outStr;chomp $outStrDesc5; chomp $outStr5;
#	my @gest = split(/\t/,$gsta);
	#my $nGenes = `wc -l $inD/assemblies/metag/genePred/genes.gff | cut -f1 -d ' '`;	chomp $nGenes;
	#die ($outStr."\n$nGenes\n");
#	$outStr .= "$nGenes\t";	$outStrDesc .= "NumGenes\t";


	return ($outStrDesc,$outStr,$outStrDesc5, $outStr5);
	#die $outStr."\n";
}

sub readG2M($){
	my ($inf) = @_;
	open K,"<",$inf;
	my %COGs2motus;
	my %motus;
	while (my $li = <K>){
		chomp($li);
		my @spl = split("\t",$li);
		$motus{$spl[0]} = $spl[3];
		$COGs2motus{$spl[0]} = $spl[2];
	}
	close K;
	print "Read gene to map file\n";
}

sub ReadsFromMapping{
	my ($tmpFile,$linkFile) = @_;
	my $newOfile = ""; my $GID = "";
	open TO,">",$linkFile;
	open I,"<",$tmpFile;
	while (my $line = <I>){
		chomp($line);
		my @spl = split("\t",$line);
		my $readID = $spl[0];
		my $geneMtch = $spl[2];
		#get read and genomeMatch
		#unless (exists($motus{$geneMtch})){print "could not find gene ID $geneMtch\n"; next;}
		#my $newHd = "@".$COGs2motus{$geneMtch}."___".$motus{$geneMtch};
		$newOfile = $GID.".rdM";
		#die ($newOfile);
		
		$readID=~m/(.*)#/;
		print TO $1."\t"."\t".$newOfile."\n";#$newHd."\n".$spl[9]."\n+\n".$spl[10]."\n";
		#die;
	}
	close TO;
	close I;
}

sub RayAssembly(){
 "mpiexec -n 1 /g/bork5/hildebra/bin/Ray-2.3.1/ray-build/Ray -o test -p test/test_1.fastq test/test_2.fastq -k 31"
 }
 
 
 sub createPsAssLongReads(){
	my ($arp1,$arp2,$singAr,$jdep, $pseudoAssFile, $Fdir, $smplName) = @_;
	my @allRds = (@{$arp1},@{$arp2},@{$singAr});
	my $psDir = $pseudoAssFile; $psDir =~ s/[^\/]+$//;
	my $psFinal = $pseudoAssFile; $psFinal =~ s/^.+\//$Fdir\//;
	#die $psFinal;
	my $pseudoAssFileFlag = $pseudoAssFile.".sto";
	
	my $psFile = $pseudoAssFile;#"$finalCommAssDir/longReads.fasta.filt";
	if (-d $psFinal){system "rm -r $psFinal";}
	if (!-e $pseudoAssFile && -e $psFinal && -e $psFinal.".sto"){
		$pseudoAssFileFlag = $psFinal.".sto"; 
		$psFile = $psFinal;
	}

	#die "@allRds\n";
	my $cmd = "";
	$cmd .= "mkdir -p $psDir $Fdir\n";
	$cmd .= "$sizFiltScr ".join(",",@allRds)." 400 -1 $pseudoAssFile\n";
	$cmd .= "$renameCtgScr $psFile $smplName\n";
	$cmd .= "touch $pseudoAssFileFlag\n";
	#die $cmd;
	my ($jname,$tmpCmd) = ("","");
	if (!-e $pseudoAssFileFlag || !-e $psFinal.".sto"){
		$jname = "_PA$JNUM";;
		($jname,$tmpCmd) = qsubSystem($logDir."pseudoAssembly.sh",$cmd,1,"5G",$jname,$jdep,"",1,\@General_Hosts,\%QSBopt) ;
	} else {$cmd="";}
	my $megDir = $psFile; $megDir =~ s/[^\/]+$//;
	#die "$megDir\n";
	return ($jname,$psFile, $megDir);
 }

 
sub spadesAssembly(){
	my ($asHr,$cAssGrpX,$nodeTmp,$finalOut,$doClean,$helpAssembl,$smplName,$hostFilter,$libInfoRef,$mateFlag) = @_;
	
	my $p1ar = $AsGrps{$cAssGrpX}{FilterSeq1};
	my $p2ar = $AsGrps{$cAssGrpX}{FilterSeq2};
	my $singlAr = $AsGrps{$cAssGrpX}{FilterSeqS};
	my $jDepe = $AsGrps{$cAssGrpX}{SeqClnDeps};
	#print all samples used 
	my $nCores = $Spades_Cores;#6
	my $noTmpOnNode = 0; #prevent usage of tmp space on node
	if ($noTmpOnNode){
		$nodeTmp = $finalOut;
	}
 	my $cmd = "rm -rf $nodeTmp\nmkdir -p $nodeTmp\nmkdir -p $finalOut\n\n";
	$cmd .= "\necho '". $AsGrps{$cAssGrpX}{AssemblSmplDirs}. "' > $nodeTmp/smpls_used.txt\n\n";
	my $defTotMem = $Spades_Memory;#60;
	my $defMem = ($defTotMem/$nCores);
	$defTotMem = $defMem * $nCores; #total really available mem (in GB)

	$cmd .= $spadesBin;
	my $K = $Spades_Kmers ;
	#insert single reads
	my $errStep = "";
	$errStep = "--only-assembler " if ($doClean == 0);
	my $numInLibs = scalar @{$p1ar};
	my $sprds = inputFmtSpades($p1ar,$p2ar,$singlAr,$logDir);
	if ($numInLibs <= 1) {$cmd .= " --meta " ;} else {$cmd .= " --sc " ;} #deactivated as never done with other T2 samples..
	$cmd .= " $K $sprds -t $nCores $errStep -m $defTotMem ";#--mismatch-correction "; # --meta  --sc "; #> $log #--meta :buggy in 3.6
	$cmd .= " --mismatch-correction " if ($spadesMisMatCor);
	if ($helpAssembl ne ""){
		$cmd .= "--untrusted-contigs $helpAssembl ";
	}
	$cmd .= "-o $nodeTmp\n";
	#from here could as well be separate 1 core job
	#cleanup assembly
	$cmd .= "\nrm -f -r $nodeTmp/K* $nodeTmp/tmp $nodeTmp/mismatch_corrector/*\n";
	#dual size filter
	$cmd .= "$renameCtgScr $nodeTmp/scaffolds.fasta $smplName\n";
	$cmd .= "$sizFiltScr $nodeTmp/scaffolds.fasta 400 200\n";
	
	if ($JNUM > 1 && $helpAssembl ne ""){
		#cluster short reads using CD-HIT, replaces scaffolds.fasta.filt2 file
		my $secAss = $nodeTmp."secondary_shorts/";
		system("mkdir -p $secAss");
		$cmd .= "cat $nodeTmp/scaffolds.fasta >> $nodeTmp/scaffolds.fasta2\n";
		$cmd .= $spadesBin." -s $nodeTmp/scaffolds.fasta2 --only-assembler -m 1000 -t $nCores $K -o $secAss\n";
		#cleanup
		$cmd .= "cp  $secAss/contigs.fasta $nodeTmp/scaffolds.fasta2\nrm -f -r $secAss\n";
	}
	$cmd .= "$assStatScr -scaff_size 500 $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.500.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.ini.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta.filt > $nodeTmp/AssemblyStats.txt\n";

	my ($cmdDB,$bwtIdx) = buildMapperIdx("$nodeTmp/scaffolds.fasta.filt",$nCores,0,$MapperProg);#$nCores);
	$cmd .= $cmdDB unless($mateFlag); #doesn't need bowtie index
	
	#clean up
	$cmd .= "\ngzip -r $nodeTmp/misc/\ngzip -r $nodeTmp/before_rr* $nodeTmp/*contigs.fa* $nodeTmp/*.fastg\n";
	$cmd .=  "rm -rf $nodeTmp/corrected\n";
	unless ($noTmpOnNode){
		$cmd .=  "mkdir -p $finalOut\ncp -r $nodeTmp/* $finalOut\n" ;
		$cmd .= "rm -rf $nodeTmp\n";
	}
	my $jname = "";
	
	if (-e $logDir."spaderun.sh.otxt"){	#check for out of mem
		open I,"<$logDir/spaderun.sh.otxt" or die "Can't open old assembly logfile $logDir\n"; my $str = join("", <I>); close I;
		if ($str =~ / Error in malloc(): out of memory/ ||$str =~ m/TERM_MEMLIMIT: job killed after reaching LSF memory usage limit/){ #memory error for real
			my $replMem  = "";
			if ($str =~ /\n    Max Memory :     (\d+) MB\n/){	$replMem = int($1*1000/$nCores*1.7);
			} elsif ($str =~ /\nMAX MEM (\d+)G\n/){	$replMem = int($1/$nCores*1.7);}
			unless ($replMem eq ""){
				if (($replMem *$nCores)< 50){$replMem = 12;} 
				$defMem = $replMem;
				print $defMem."G: new MEM\n"; #die $defMem."\n";
			}
		}
	}
	$cmd .= "echo \"MAX MEM ".$defTotMem."G\"";
	#print "in Assembly\n$jDepe\n";
	#print "$finalOut/scaffolds.fasta.filt\n";
	unless (-e "$finalOut/scaffolds.fasta.filt" && !-z "$finalOut/scaffolds.fasta.filt" && !-z "$finalOut/AssemblyStats.txt"){
		#my $size_in_mb = (-s $fh) / (1024 * 1024);
		my $tmpCmd="";
		$jname = "_ASS$JNUM";#$givenJName;
		if ($hostFilter || $SpadesAlwaysHDDnode){
			my $tmpSHDD = $QSBopt{tmpSpace};
			$QSBopt{tmpSpace} = $Spades_HDspace."G" unless ($Spades_HDspace =~ m/G$/); #set option how much tmp space is required, and reset afterwards
			($jname,$tmpCmd) = qsubSystem($logDir."spaderun.sh",$cmd,(int($nCores/2)+1),int($defMem*2)."G",$jname,$jDepe,"",1,\@Spades_Hosts,\%QSBopt) ;
			$QSBopt{tmpSpace} = $tmpSHDD;
		} else {
			($jname,$tmpCmd) = qsubSystem($logDir."spaderun.sh",$cmd,(int($nCores/2)+1),int($defMem*2)."G",$jname,$jDepe,"",1,\@General_Hosts,\%QSBopt) ;
		}
	} else {
		print "Assembly still on tmp dir\n";
	}
	#die("SPADE\n");
	return ($jname);
}
#GenePrediction
sub run_prodigal($ $ $ $) {
	my ($inputScaff, $outDir, $jobDepend,$finDir,$specJname) = @_;
	system("mkdir -p $outDir");
	my $prodigalBin = "/g/bork5/hildebra/bin/Prodigal-2.6.1/prodigal";
	my $output_format_prodigal = "gff"; 
	my $expectedD = "$finDir/assemblies/metag/genePred";
	my $prodigal_cmd = ("$prodigalBin -i $inputScaff -o $outDir/genes.gff -a $outDir/proteins.faa -d $outDir/genes.fna -f $output_format_prodigal -p meta\n");
	$prodigal_cmd .= "sleep 2;cut -f1 -d \" \" $outDir/proteins.faa > $outDir/proteins.shrtHD.faa\n" ;
	$prodigal_cmd .= "cut -f1 -d \" \" $outDir/genes.fna > $outDir/genes.shrtHD.fna\n" ;
	$prodigal_cmd .= "rm -f $outDir/proteins.faa $outDir/genes.fna";
	
	#if (!-e "$inputScaff") {die "Input scaffold $inputScaff does not exists.\n";}
	#print $prodigal_cmd."\n";
	#if (system $prodigal_cmd){die "Finished prodigal with errors";}
	my $jname = "";
	#print "$expectedD/proteins.shrtHD.faa\n$expectedD/genes.gff\n";
	if ( (!-s "$expectedD/proteins.shrtHD.faa" || !-s "$expectedD/genes.gff")){
		if (!-s "$outDir/proteins.shrtHD.faa" || !-s "$outDir/genes.gff"){
			$jname = "_GP$JNUM"; my $tmpCmd="";
			$jname = $specJname.$JNUM if ($specJname ne "");
			#print "genesubm\n";
			($jname,$tmpCmd) = qsubSystem($logDir."prodigalrun.sh",$prodigal_cmd,1,"1G",$jname,$jobDepend,"",1,\@General_Hosts,\%QSBopt); 
		}
	}
	return $jname;
}
sub run_prodigal_augustus($ $ $ $) {
	my ($inputScaff, $outDir, $jobDepend,$finDir,$specJname,$GlbTmpPath) = @_;
	my $tmpGene = $outDir."/tmpCalls/";
	system("mkdir -p $outDir");
	my $scrDir = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/euk_gene_caller/bin";
	my $output_format_prodigal = "gff"; 
	my $expectedD = "$finDir/genePred";
	my $augCmd  = ""; my $cmpCmd = "";
	my $tmpCmd=""; #placeholder for qsub return
	my $splitDep = $jobDepend;
	my $inputBac = $inputScaff; my $inputEuk = "$GlbTmpPath/euk.kraken.fasta";
	my $bacmark="";
	if ($DO_EUK_GENE_PRED ){ $bacmark = ".bac";	}
	
	#print "$outDir/proteins$bacmark.shrtHD.faa\n";
	
	if ( (-s "$expectedD/proteins$bacmark.shrtHD.faa" || -s "$expectedD/genes$bacmark.gff") &&
			(-s "$outDir/proteins$bacmark.shrtHD.faa" || -s "$outDir/genes$bacmark.gff") ){
		return "";
	}
	#use kraken to classify contigs..
	if ($DO_EUK_GENE_PRED ){ #stupid way of doing things..
		#first split euk vs bac contigs
		my $splCores = 10;
		my ($refDB,$shrtDB,$clnCmd) = prepDiamondDB("NOG",$globaldDiaDBdir,$splCores);
		my $scmd = "$splitKgdContig $inputScaff $GlbTmpPath $krakenDBDirGlobal $globaldDiaDBdir $splCores\n";
		$scmd .= "touch $GlbTmpPath/bac.euk.split.sto\n";
		#die $scmd."\n";
		$inputBac = "$GlbTmpPath/bact.kraken.fasta";
		#die "$inputEuk || !-e $inputBac\n";
		if (!-e $inputEuk || !-e $inputBac || !-e "$GlbTmpPath/bac.euk.split.sto"){
			($splitDep,$tmpCmd) = qsubSystem($logDir."splitAssembly.sh",$scmd,$splCores,int(50/$splCores)."G","_SC$JNUM","$jobDepend;$krakDeps","",1,\@General_Hosts,\%QSBopt); 
		}
	}

#die "toofar";
	my $prodigal_cmd = "";
	$prodigal_cmd .= "mkdir -p $tmpGene\n";
		
		$prodigal_cmd.="if [ ! -s $inputBac ];then\n";
		$prodigal_cmd.="touch $tmpGene/proteins$bacmark.shrtHD.faa $tmpGene/genes$bacmark.per.ctg $tmpGene/genes$bacmark.shrtHD.fna $tmpGene/genes$bacmark.gff\n";
		$prodigal_cmd.="else\n";
	
	$prodigal_cmd .= ("$prodigalBin -i $inputBac -o $tmpGene/genes$bacmark.gff -a $tmpGene/proteins.faa -d $tmpGene/genes.fna -f $output_format_prodigal -p meta\n");
	$prodigal_cmd .= "sleep 2;cut -f1 -d \" \" $tmpGene/proteins.faa > $tmpGene/proteins$bacmark.shrtHD.faa\n" ;
	$prodigal_cmd .= "cut -f1 -d \" \" $tmpGene/genes.fna > $tmpGene/genes$bacmark.shrtHD.fna\n" ;
	$prodigal_cmd .= "cut -f1 $tmpGene/genes$bacmark.gff | sort | uniq -c | grep -v '#' |  awk -v OFS='\\t' {'print \$2, \$1'} > $tmpGene/genes$bacmark.per.ctg\n";
	$prodigal_cmd .= "rm -f $tmpGene/proteins.faa $tmpGene/genes.fna\n";	
	#copy everything to the right place (for now, might have to change this later to take care of Eukarya)
	$prodigal_cmd .= "mkdir -p $expectedD\n";
		
		$prodigal_cmd.="fi\n";
	
	
	
	if ($DO_EUK_GENE_PRED ){ #stupid way of doing things..
		#1st set mark that prodigal genes are on bacteria only
		$prodigal_cmd .= "touch $tmpGene/genes.prodigal.bact.only\n";
		$augCmd .= "mkdir -p $tmpGene\n";
		my $eukmark = ".euk";

		$augCmd.="if [ ! -s $inputEuk ];then\n";
		$augCmd.="touch $tmpGene/genes$eukmark.gff $tmpGene/genes$eukmark.codingseq $tmpGene/genes$eukmark.fna\n";
		$augCmd.="else\n";
		$augCmd .= "$augustusBin --species=ustilago_maydis $inputEuk --protein=off --codingseq=on > $tmpGene/augustus.gff;\n";
		$augCmd .= "perl $scrDir/getAnnoFast.pl --seqfile=$inputEuk $tmpGene/augustus.gff;";
		$augCmd.="mv $tmpGene/augustus.gff $tmpGene/genes$eukmark.gff\n mv $tmpGene/augustus.codingseq $tmpGene/genes$eukmark.codingseq\nmv $tmpGene/augustus.cdsexons $tmpGene/genes$eukmark.fna\n";
		$augCmd .= "cut -f1 -d \" \" $tmpGene/genes$eukmark.fna > $tmpGene/genes$eukmark.shrtHD.fna\n" ;
		$augCmd .= "rm $tmpGene/genes$eukmark.fna\n";
		$augCmd.="fi\n";

		
		my $axl = "$tmpGene/augustus.list"; my $pxl = "$tmpGene/genes.list";
		my $cmpCmd = "grep -v \"^#\" $tmpGene/genes.gff > $pxl;\n";
		$cmpCmd .= "grep -w \"transcript\" $tmpGene/augustus.gff > $axl;\n";
		$cmpCmd .= "grep \">\" $inputEuk > $tmpGene.scaff.list;\n";
		$cmpCmd .= "awk \'!/^>/ { printf \"%s\", $0; n = \"\n\" } /^>/ { print n $0; n = \"\" } END { printf \"%s\", n }\' $tmpGene/augustus.codingseq > $axl.nolines.fa;\n";
		$cmpCmd .= "awk \'!/^>/ { printf \"%s\", $0; n = \"\n\" } /^>/ { print n $0; n = \"\" } END { printf \"%s\", n }\' $tmpGene/genes.shrtHD.fna > $pxl.nolines.fa;\n";
		$cmpCmd .= "perl $scrDir/process_nodes.pl $pxl $axl;\n";
		$cmpCmd .= "perl $scrDir/match_nodes_euk.pl $axl.output.txt $axl.nolines.fa;\n";
		$cmpCmd .= "perl $scrDir/match_nodes_bac.pl $pxl.output.txt $pxl.nolines.fa;\n";
		$cmpCmd .= "perl $scrDir/overlap_nodes.pl $tmpGene.scaff.list $axl.output.txt $pxl.list.output.txt ;\n"; #$f.blast.matched.txt
		$cmpCmd .= "perl $scrDir/finalgeneparsing.pl $tmpGene.scaff.list.ALLmatched.txt\n";
		#now, let's pull the appropriated gene. for each contig, need to retrieve ALL genes for that contig.
		$cmpCmd .= "cut -f 4 $tmpGene.scaff.list.ALLmatched.txt.decision.txt | sort |uniq> $tmpGene/formatching.txt\n";
		#$cmpCmd .= "perl ./bin/getgenes.pl $tmpGene/formatching.txt $f.aug.nolines.fa.renamed.fa \n";
		#$cmpCmd .= "perl ./bin/getgenes.pl $tmpGene/formatching.txt $f.mgm.nolines.fa.renamed.fa\n";
		#$cmpCmd .= "cat $f.aug.nolines.fa.renamed.fa.output.txt $f.mgm.nolines.fa.renamed.fa.output.txt > $f.genecalled.fa\n";
		#die $augCmd;
	}
	
	#finish up file transfer etc
	$prodigal_cmd .= "mv $tmpGene/* $outDir\nrm -r $tmpGene\n";

	#if (!-e "$inputScaff") {die "Input scaffold $inputScaff does not exists.\n";}
	#print $prodigal_cmd."\n";
	#if (system $prodigal_cmd){die "Finished prodigal with errors";}
	my $jname = "";
	#print "$expectedD/proteins.shrtHD.faa\n$expectedD/genes.gff\n";
	if ( (!-s "$expectedD/proteins$bacmark.shrtHD.faa" || !-s "$expectedD/genes$bacmark.gff") &&
			(!-s "$outDir/proteins$bacmark.shrtHD.faa" || !-s "$outDir/genes$bacmark.gff") ){
		$jname = "_GP$JNUM"; 
		$jname = $specJname.$JNUM if ($specJname ne "");
		#print "genesubm\n";
		($jname,$tmpCmd) = qsubSystem($logDir."prodigalrun.sh",$augCmd.$prodigal_cmd,1,"1G",$jname,$splitDep,"",1,\@General_Hosts,\%QSBopt); 

	}
	return $jname;
}

sub filterSizeFastaCmd(){
my ($fasta,$siz,$siz2) = @_;
my $cmdd= ("/g/bork3/home/hildebra/dev/Perl/assemblies/./sizeFilterFas.pl $fasta $siz $siz2");
return ($fasta.".filt",$cmdd);
}








