#!/usr/bin/env perl
#uses a multi sample assembly to bin contigs with metabat
# ./compoundBinning.pl /g/scb/bork/hildebra/SNP/GNMassSimuC
# ./compoundBinning.pl /g/scb/bork/hildebra/SNP/GNMass3/alien-11-374-0
use warnings;
use strict;
use Data::Dumper;
use Mods::GenoMetaAss qw(readMap qsubSystem emptyQsubOpt  unzipFileARezip);
use Mods::IO_Tamoc_progs qw(jgi_depth_cmd );

#add this? https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4748697/figure/fig-1/

sub runMetaBat;

my $inD = $ARGV[0];
$inD.="/" unless($inD =~ m/\/$/);
my $doSubmit = 1;

#set up basic structures
my $QSBoptHR = emptyQsubOpt($doSubmit,"");
my $tmpD = $inD."/tmp/";
my $outD = $inD."/Binning/MetaBat/";
my $logDir = $outD."LOGandSUB/";
my $singleSample = 0;
if (-e "$inD/LOGandSUB/inmap.txt"){ #this is the outdir of a whole TAMOC run
	$singleSample =0;
	print "Compound Assembly MetaBatting..\n";
} else { #single sample output
	$singleSample=1;
	print "Single Sample MetaBatting..\n";
}

system "mkdir -p $tmpD $outD $logDir";
if ($singleSample){
	
	#get map
	my $SmplNm = `cat $inD/mapping/done.sto`;
	$SmplNm =~ s/-smd.bam\n?//;
	#my $inBAM = $inD."mapping/$SmplNm-smd.bam";
	my $jgiCov = $inD."mapping/$SmplNm-smd.bam.jgi.cov";
	#print "$jgiCov\n";
	my $jgiConn = $inD."mapping/$SmplNm-smd.bam.jgi.pairs.sparse";
	my ($befZ,$aftZ) = unzipFileARezip( [$jgiCov,$jgiConn] );
	my $metaGD = `cat $inD/assemblies/metag/assembly.txt`; chomp $metaGD;

	my $scaffs = $metaGD."/scaffolds.fasta.filt";
	
	die runMetaBat($befZ,$aftZ,$inD."mapping/$SmplNm-smd.bam.jgi",$outD,$SmplNm,$scaffs);
	my $ess100 = $metaGD."/ContigStats/ess100genes/ess100.id.txt";

	#die $befZ.$aftZ;

} else {
	#die Dumper($hrm);
	my ($hrm,$asGrpObj) = readMap($inD."LOGandSUB/inmap.txt");
	my %map = %{$hrm};
	my %assGrp = %{$asGrpObj};
	#infer Assembly dirs & corrsponding bams with several Samples (compound assemblies)
	my @smpls = @{$map{smpl_order}};
	my %DOs;
	#figure out which compound assemblies there are..
	foreach my $smpl (@smpls){
		my $cntAim = $assGrp{ $map{$smpl}{MapGroup} }{CntAimMap};
		#die "XX".$cntAim."\n";
		if ( $cntAim > 1){
			my $dir2rd = $inD.$map{$smpl}{dir};
			my $cAssGrp = $map{$smpl}{AssGroup};
			my $tar = "AssmblGrp_$cAssGrp";
			if (!-e "$inD/$tar"){
				print "Can't read $inD/$tar\nSkipping..\n";
				next;
			}
			#(-s "$curOutDir/mapping/Align_ment-smd.bam" || -s "$curOutDir/mapping/Align_ment-smd.cram")
			#my $bam = $inD."mapping/Align_ment-smd.bam";

			if ( !exists($DOs{$tar}{dir}) ){
				$DOs{$tar}{dir} = $dir2rd;
				$DOs{$tar}{SmplID} = $map{$smpl}{SmplID};
			} else {
				$DOs{$tar}{dir} .= ",".$dir2rd;
				$DOs{$tar}{SmplID} .= ",".$map{$smpl}{SmplID};
			}
		}
	}
	#run metabat on each assembly group
	foreach my $Doo (keys %DOs){
		my $bef = "";
		$bef = jgi_depth_cmd($DOs{$Doo}{dir},$tmpD."/depth",95) unless (-e );
		system runMetaBat($bef,"","$tmpD/depth.jgi",$outD,"test","$inD/$Doo/metag/scaffolds.fasta.filt");
	}
}



sub runMetaBat($ $ $ $ $ $){
	my ($before,$after,$jgO,$outD,$nm, $fna) = @_;
	print "Running MetaBat..\n";
	my $numCore = 6;
	my $mbBin = "/g/bork3/home/hildebra/bin/metabat/./metabat";
	my $outCtgNms = "$outD/$nm.ctgs.txt";
	my $outFna = "$outD/$nm.fasta.fna";
	my $outMat = "$outD/$nm.mat.txt";
	#my $jgO = "$tmp/depth.jgi";
	if (!-e $fna){die "Can't find requried scaffold file in metabat routine: $fna\n"; }
	#start metabat
	my $cmd = $before."\n";
	$cmd .= "$mbBin -i $fna -a $jgO.depth.txt -o $outFna -p $jgO.pairs.sparse -l $outCtgNms --minCVSum  10 -t $numCore --saveCls $outMat -v\n";
	$cmd .= $after;
	$cmd .= "echo \"$nm\" > $outD/MeBa.sto";
	return $cmd;
	#die $cmd."\n";
	#my $jobName = "MBat";
	#my ($jobNameX, $tmpCmd) = qsubSystem($logDir."metaBat.sh",$cmd,$numCore,"60G",$jobName,"","",1,[],$QSBoptHR);

}
