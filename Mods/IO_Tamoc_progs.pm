package Mods::IO_Tamoc_progs;
use warnings;
use Cwd 'abs_path';
use strict;
#TAMOC programs related to IO to other programs, program paths .. not real subroutines that do anything

use Exporter qw(import);
our @EXPORT_OK = qw(inputFmtSpades inputFmtMegahit jgi_depth_cmd createGapFillopt getProgPaths buildMapperIdx);


sub getProgPaths{
	my @var = @_;
	my $srchVar = $var[0] ;
	my @multVars = ();my @retA;
	if (ref $srchVar eq 'ARRAY') {
		print "ARRAY\n";
		@multVars = @{$srchVar};
	}
	my $required=1;
	#print $_ . " => " . $INC{$_} . "\n" for keys %INC;
	
	my $modDir = $INC{"Mods/IO_Tamoc_progs.pm"};
	$modDir =~ s/IO_Tamoc_progs.pm//;
	if (@var > 1){$required = $var[1];}
	#die "$srchVar\n$required\n@multVars\n";
	#my $optF = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/Mods/MATAFILERcfg.txt";
	my $optF = "$modDir/MATAFILERcfg.txt";
	open I,"<$optF" or die "Can't open $optF\n";
	my $TMCpath = "";my $Tset=0;
	while (<I>){
		chomp;
		next if (m/^#/);
		if (!$Tset && m/^MFLRDir\t([^#]+)/){
			$Tset=1;$TMCpath = $1;
#			next;
		}
		if (@multVars > 0){
			for (my $k=0;$k<@multVars;$k++){
				if (m/^$multVars[$k]\t([^#]+)/){
					my $reV = $1;
					$reV =~s/\[MFLRDir\]/$TMCpath/;
					$retA[$k] = $reV;
				}
			}
		} elsif (m/^$srchVar\t([^#]+)/){
			my $reV = $1;
			$reV =~s/\[MFLRDir\]/$TMCpath/;
			close I;return $reV;
		}
	}
	close I;
	if (@multVars > 0){
		return \@retA;
	}
	die "Can't find program $srchVar\n" if ($required);
	return "";
}

sub buildMapperIdx($ $ $ $){
	my ($REF,$ncore,$lrgDB,$MapperProg) = @_;
	my $bwt2Bin = getProgPaths("bwt2");
	my $bwaBin  = getProgPaths("bwa");
	my $bwtIdx = $REF.".bw2";
	my $dbCmd =$bwt2Bin."-build ";
	$dbCmd .= " --large-index "if ($lrgDB);
	$dbCmd .= "-q $REF --threads $ncore $bwtIdx\n";
	if (-s $REF.".bw2.1.bt2" || -s $REF.".bw2.1.bt2l"){$dbCmd = "";} #deactivate if already built
	if($MapperProg==2) { 
		$dbCmd = $bwaBin." index $REF";
		if (-s $REF.".pac"){$dbCmd = "";} 
		#die $dbCmd."\n";
	}
	#my $jobN = "_ASDB$JNUM"; my $tmpCmd;
	#($jobN, $tmpCmd) = qsubSystem($logDir."BAM2CRAMxtra.sh",$dbCmd,1,"10G",$jobN,"","",1,[],\%QSBopt);
	return ($dbCmd,$bwtIdx);
}


sub inputFmtSpades($ $ $ $){
	my ($p1ar,$p2ar,$singlAr,$logDir) = @_;
	my @p1 = @{$p1ar}; my @p2 = @{$p2ar};
	my @singl = @{$singlAr};
	my $doYAML = 0;
	if (@p1 > 2){$doYAML=1;}
	if (@p1 != @p2){print "Unequal paired read array lengths arrays for Spades\n"; exit(2);}
	my $sprds = "";
	if ($doYAML==0){
		for (my $i =0; $i<@p1;$i++){
			$sprds .= " --pe".($i+1) ."-1 $p1[$i] --pe".($i+1) ."-2 $p2[$i]";
		}
		for (my $i=0;$i<@singl;$i++){
			$sprds .= " --pe".($i+1) ."-s $singl[$i]";
		}
	} else {
		open O,">$logDir/spadesInput.yaml";
		print O "  [\n      {\n        orientation: \"fr\",\n        type: \"paired-end\",\n        left reads: [\n";
		for (my $i =0; $i<@p2;$i++){
			if (($i+1) == @p2){	print O "          \"$p1[$i]\"\n";	} else { print O "          \"$p1[$i]\",\n";}
		}
		print O "     ],\n        right reads: [\n";
		for (my $i =0; $i<@p1;$i++){
			if (($i+1) == @p1){	print O "          \"$p2[$i]\"\n";	} else { print O "          \"$p2[$i]\",\n";}
			#print O "          \"$p2[$i]\"\n";
		}
		print O "       ]\n      }";
		if (@singl > 0){
			print O ",\n      {\n        type: \"single\",\n        single reads: [\n";
			for (my $i=0;$i<@singl;$i++){
				if (($i+1) == @singl){	print O "          \"$singl[$i]\"\n";	} else { print O "          \"$singl[$i]\",\n";}
				#print O  "          \"$singl[$i]\"\n";
			}
			print O "       ]\n      }";
		}
		#end of yaml file
		print O "\n     ]\n";
		close O;
		$sprds = " --dataset $logDir/spadesInput.yaml";
	}
	return $sprds;
 }

 sub inputFmtMegahit($ $ $ $){
	my ($p1ar,$p2ar,$singlAr,$logDir) = @_;
	my @p1 = @{$p1ar}; my @p2 = @{$p2ar};
	my @singl = @{$singlAr};
	if (@p1 != @p2){print "Unequal paired read array lengths arrays for Spades\n"; exit(2);}
	my $sprds = "";
	if (@p1 > 0){ 
		$sprds .= "-1 ".join(",",@p1) . " -2 ".join(",",@p2)." ";
	}
	if (@singl > 0){
		$sprds .= "-r ".join(",",@singl);
	}
	return $sprds;
 }


sub jgi_depth_cmd($ $ $){
	my ($dirs,$out,$perID) = @_;
	#die $dirs."\n";
	my @dirSS = split(',',$dirs);
	#go through each dir and find sample name
	my $comBAM = "";
	foreach my $DDI (@dirSS){
		if ( $DDI =~ m/\/$/ ||  $DDI !~ m/bam$/ ){
			my $SmplNm = `cat $DDI/mapping/done.sto`;#$SmplNm =~ s/-smd.bam\n?//;
			chomp $SmplNm;
			$comBAM .= "$DDI/mapping/$SmplNm ";
		} else {
			$comBAM .= "$DDI ";
		}
	}
	#my $comBAM = join("/mapping/Align_ment-smd.bam ",@dirSS);
	my $covCmd = "";
	$covCmd .= "rm -f $out.jgi.*\n";
	$covCmd .= "\n/g/bork3/home/hildebra/bin/metabat/./jgi_summarize_bam_contig_depths";
	$covCmd .= " --outputDepth $out.jgi.depth.txt  --percentIdentity $perID --pairedContigs $out.jgi.pairs.sparse $comBAM > $out.jgi.cov\n";
	$covCmd .= "gzip $out.jgi*\n";
	if (-e "$out.jgi.cov"){$covCmd="";}

	#$covCmd .= "gzip $nxtBAM.jgi*\n";
	return $covCmd;
}

#2nd: arrray of files, paired sep by ","
sub createGapFillopt($ $ $){
 my ($ofile, $arFiles, $insertSizAr) = @_;
 my @Files = @{$arFiles};
 my @insertSiz = @{$insertSizAr};
 my $opt = "";
 for (my $cnt=0;$cnt<@Files; $cnt++){
	my $line = "LIB$cnt bwa ";
	my $lineMode = 2; 
	if ($lineMode==2){ #PE; only available option at the moment
		my @curFils = split(",",$Files[$cnt]);
		die ("Only paired reads accepatble to GapFiler.\n") if (@curFils != 2);
		$line .= $curFils[0]." ".$curFils[1];
	}
	$line .= " $insertSiz[$cnt] 0.3 FR\n";
	$opt .= $line;
 }
 #print $opt."\n";
 open O,">",$ofile; print O $opt; close O;
}
