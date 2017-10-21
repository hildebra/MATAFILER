package Mods::GenoMetaAss;
use warnings;
use Cwd 'abs_path';
use strict;
#use List::MoreUtils 'first_index'; 
use Mods::IO_Tamoc_progs qw(getProgPaths);

use Exporter qw(import);
our @EXPORT_OK = qw(convertMSA2NXS gzipwrite renameFastaCnts renameFastqCnts readNCBItax gzipopen readMap 
		readMapS renameFastHD findQsubSys emptyQsubOpt qsubSystem
		readClstrRev  unzipFileARezip systemW is_integer readGFF reverse_complement reverse_complement_IUPAC
		readFasta writeFasta readFastHD readTabByKey convertNT2AA prefix_find runDiamond median deNovo16S);





sub first_index (&@) {
    my $f = shift;
    for my $i (0 .. $#_) {
	local *_ = \$_[$i];	
	return $i if $f->();
    }
    return -1;
}

sub readFasta{
	my $fil = $_[0];
	my $cutHd=0;
	my $sepChr= "\\s";
	$cutHd = $_[1] if (@_ > 1);
	$sepChr = $_[2] if (@_ > 2);
	my %Hseq;
	if (-z $fil){ return \%Hseq;}
	my $FAS;
	open($FAS,"<","$fil") || die ("Couldn't open FASTA file $fil\n");
	#my $FAS = gzipopen($fil,"fasta file to readFasta",1);
	my $temp; my $line; 
	my $trHe =<$FAS>; chomp ($trHe);  $trHe = substr($trHe,1);
	if ($cutHd) {$trHe =~ s/$sepChr.*//;}
	while($line = <$FAS>){
		chomp($line);
		if ($line =~ m/^>/){
			#finish old fas`
			$Hseq{$trHe} = $temp; 
			#prep new entry
			
			$trHe = substr($line,1);# $trHe =~ s/;size=\d+;.*//; 
			if ($cutHd) {$trHe =~ s/$sepChr.*//;}
			$trHe = "".$trHe; #just to ensure it's a string
			$temp = "";
			#die $trHe."\n$sepChr\n";
			next;
		}
		$temp .= ($line);
	}
	$Hseq{$trHe} = $temp;
	close ($FAS);
	return \%Hseq;
}



sub convertNT2AA($){
 my ($text) = @_;
 # translate a DNA 3-character codon to an amino acid. We take three letter groups from
# the incoming string of C A T and G and translate them via a Hash.  It's listed vertically
# so we can label each of the amino acids as we set it up.

#$text = "aaatgaccgatcagctacgatcagctataaaaaccccggagctacgatcatcg";

$text =~ s/[N]*$//i;
if (length($text) % 3 != 0){
	my $ltr = length($text) % 3;
	my $lt = length($text);
	$text = substr $text,0,($lt -$ltr);
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

my %convertor = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

# We don't actually know where the groups of 3 will start in our sample piece of DNA, so we've got
# three ways of doing the coding ... here's a loop to work out each of the possibilities in turn,
# leaving the odd extra letters on the beginning or end.
#for ($s=0; $s<3; $s++) {
#last; #bs, just need on code
#        $scrap = substr($text,0,$s);
#        $main = substr($text,$s);
#        $main =~ s/(...)/"$convertor{uc $1}" || "?"/eg;
#        print "$scrap$main\n";
#        }
$text =~ s/(...)/"$convertor{uc $1}" || "?"/eg;
$text =~ s/[^ACDEFGHIKLMNPQRSTVWY*]/X/g;
#die $text;
return $text;
}



 sub deNovo16S{
	#returns 1)hash with rDNAs and 2) name of best seq
	my ($fasFile,$outfile) = @_;
	die "Fasta $fasFile doesn't exist\n" unless (-e $fasFile);
	my %ret; my $refSize = 1580;
	print "Detecting de novo 16S\n";
	
	my $newGff = "$outfile.gff";
	system "rm -f $newGff" if (-e$newGff);
	my $rnaBin = getProgPaths("rnammer");
	my $cmd = "$rnaBin -S bac -m ssu -gff $newGff < $fasFile";
	system $cmd."\n";
	#FP929038        RNAmmer-1.2     rRNA    507757  508732  1202.9  -       .       16s_rRNA
	my $genoR = readFasta($fasFile,1); my %geno = %{$genoR};
	#my @k = keys(%geno); die @k.join(" ",@k)."\n";
	open II,"<",$newGff; my $gffcnt=0;
	while(<II>){next if (/^#/ || length($_) < 1);my @spl = split(/\s+/);
		#print $_;
		my $newS = substr($geno{$spl[0]},$spl[3],$spl[4]-$spl[3]);
		if ($spl[6] eq "-"){$newS = reverse_complement_IUPAC($newS);}
		#die $newS;
		$ret{$spl[0]."_rrn_".$gffcnt} = $newS;
		$gffcnt++;
	}
	close II;
	#check for new better 16S
	my $bestV=100000;
	my @head = keys(%ret);	my $best = 0; 	my $cnt =0;
	if (@head == 0){
		#die "Empty array\n$newGff\n\n$cmd\n";
		print "Empty 16S array\n$newGff\n";
		return (\%ret,"");
	} else {
		foreach my $hd (@head){
			my $tmp = abs($refSize - length($ret{$hd}));
			if ($tmp < $bestV){$bestV = $tmp; $best = $cnt;}
			$cnt++;
			#print $tmp."\n";
		}
	}
	#print "16S ".$bestV."\n";
	my %ret3 = ($head[$best],$ret{$head[$best]});
	return (\%ret,$head[$best]);
}



sub prefix_find($){
	my ($ar) = @_;
	my @rds = @{$ar}; my @newRds = @rds;
	my $first = $rds[0]; 
	#die "$first FFF\n";
	for (my $i=0;$i<@rds; $i++){
		my $second = $rds[$i];
		my @matches;
		for my $start (0..length $first) {
			for my $len ($start+1 .. length $first) {
				my $substr = substr($first, $start, $len);
				push @matches, $second =~ m[($substr)]g;
			}
		}
		my ($len, $longest) = 0;
		length > $len and ($longest, $len) = ($_, length) for @matches;
		$first = $longest;
		#print "$longest $i\n";
		#"$first\0$second" =~ m/^(.*)\0\1/s;
		#$first = $1;
	}
	
	#print "$first common prefix\n";
	#for (my $i=0; $i<@newRds; $i++){
	#	$newRds[$i] =~ s/^$first//;
	#	$newRds[$i] .= "/" if ($newRds[$i] !~ m/\/$/ && length($newRds[$i] ) != 0 );
	#}
	return ($first);
}

 
sub is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}


sub gzipwrite{
	my ($outF,$descr) = @_;
	$outF .= ".gz" if ( $outF !~ m/\.gz$/);
	open (my $O, "| gzip -c > $outF") or die "error starting gzip pipe $outF\n$!";
	return $O;
}
sub gzipopen{
	my ($inF,$descr) = @_;
	my $dodie = 1;
	if (@_ > 2){$dodie = $_[2];}
	$inF .= ".gz" if (!-e $inF && -e $inF.".gz");
	if ($inF =~ m/\.gz$/){
		my $inFwo = $inF; $inFwo =~ s/\.gz$//;
		$inF = $inFwo if (-e $inFwo && !-e $inF);
	} else {
		$inF .= ".gz" if (!-e $inF && -e $inF.".gz");
	}
	my $I; my $OK = 1;
	if($inF =~ m/\.gz$/ ){
		my $msg = "Can't open a pipe to $descr file $inF\n";
		if (!open($I, "gunzip -c $inF |")) {if ($dodie){die $msg;} else { $OK=0;print $msg;}}
	} else{
		my $msg = "Can't open $descr file $inF\n";
		if (!open($I, "<", "$inF") ) {if ($dodie){die $msg;} else {$OK=0; print $msg;}}
	}
	return ($I,$OK);
}

sub readNCBItax($){
	#reads in Jaime's ete tax file
	my ($tIn)  = @_;
	open I,"<$tIn";
	while (my $l = <I>){
		chomp $l;
		my @spl = split (/\t/,$l);
		my $txID = $spl[0];
		my $rnksLvl = $spl[2];
		my $rnks = $spl[3];
	}
	close I;
}

sub readFastHD($){#only reads headers
	my ($inF) = @_;
	my @ret;
	open I,"<$inF" or die "can t open $inF\n";
	while (my $l = <I>){
		if ($l =~ m/^>(\S+)/){
			push @ret,$1;
		}
	}
	close I;
	return \@ret;
}

sub renameFastHD($ $ $){ #set a new name for headers in fasta files
	my ($inF,$hr,$extr) = @_;
	my %COGid = %{$hr};
	my $oS = "";
	my %CidCnt;
	open I,"<$inF" or die "can t open renameFastHD  $inF\n";
	while (my $l = <I>){
		if ($l =~ m/^>/){
			chomp $l;
			unless (exists $COGid{$l}){die "Can't find $l in COGs\n";}
			if (exists $CidCnt{$COGid{$l}}){$CidCnt{$COGid{$l}}++;
				$oS .= ">".$extr."_".$COGid{$l}."..".$CidCnt{$COGid{$l}}."\n";
			} else {$CidCnt{$COGid{$l}} = 0;
				$oS .= ">".$extr."_".$COGid{$l}."\n";
			}
			

		} else {
			$oS .= $l;
		}
	}
	close I;
	open O,">$inF" or die "Can't open rename Fasta HD out file $inF\n";
	print O $oS;
	close O;
}
sub renameFastqCnts($ $){ #set a new name for headers in fastq files, using a simple scheme of prefix and just counts afterwards
	my ($inF,$prefix) = @_;
	open I,"<$inF" or die "can t open $inF\n";
	my $cnt  = 0; my $cnt2 = 0;
	open O,">$inF.tmp"; my $plusSeen = 1;
	while (my $l = <I>){
		if ($l =~ m/^@/ && $plusSeen && $cnt2 >= 4 ){ #$l =~ m/^@/ && 
			print O "@".$prefix."_".$cnt."\n";
			$cnt ++; $plusSeen=0;$cnt2=0;
#			print L ">".$prefix."_".$cnt."\t$l";
		} else {
			print O $l;
			$plusSeen=1 if ($l =~ m/^\+\n$/);
		}
		$cnt2 ++;
	}
	close I; close O;
	system "rm $inF;mv $inF.tmp $inF";
}

sub renameFastaCnts($ $ $){ #set a new name for headers in fasta files, using a simple scheme of prefix and just counts afterwards
	my ($inF,$prefix,$logF) = @_;
	open I,"<$inF" or die "can t open $inF\n";
	my $cnt  = 0;
	open O,">$inF.tmp";
	open L,">$logF";
	while (my $l = <I>){
		if ($l =~ m/^>/){
			print O ">".$prefix."_".$cnt."\n";
			print L ">".$prefix."_".$cnt."\t$l";
		} else {
			print O $l;
		}
		$cnt ++;
	}
	close I; close O; close L;
	system "rm $inF;mv $inF.tmp $inF";
}


sub readClstrRev($){ #gets the exact assembled genes clustered in GC genes #version for my shortened index file
	my ($inF) = @_;
	my %retR; my %retF;
	print "Reading clstr $inF.. \n";
	open I,"<$inF" or die "Can't find rev clustering file $inF\n"; my $curCl="";
	while (<I>){
		chomp; 
		next if (m/^#/ || length($_) < 5);
		my @arr = split /\t/;
		#my $pos = index($_, "\t");
		$curCl = $arr[0];#substr($_,0,$pos);
		#my $rem = $arr[1]; #substr($_,$pos+1);
		#print $curCl."XX$rem\n";
		#foreach (split /,/,$rem){$retR{$_} = $curCl;}
		my @tmpArr = split /,/,$arr[1];
		@retR{@tmpArr} = ($curCl) x @tmpArr;
		$retF{$curCl} = $arr[1];

	}
	close I;
	print "done\n";
	return (\%retR,\%retF);
}




sub unzipFileARezip($){
	my ($inFar) = @_;
	my @inFs = @{$inFar}; my $totCnt=0;
	my $bef = "gunzip "; my $aft="gzip ";
	#print "$inF.gz\n";
	foreach my $inF (@inFs){
		if (-e "$inF.gz"){
			$bef .= " $inF.gz"; $aft .= " $inF"; $totCnt++;
		}
	}
	if ($totCnt==0) {
		return ("","");
	} else {
		return $bef."\n",$aft."\n";
	}
}
sub reverse_complement_IUPAC ($) {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}


sub readGFF($){
	my ($inF) =@_;
	my %ret;my $sbcnt=1;
	open I,"<$inF";
	while (<I>){
		if (m/^#/){$sbcnt=1;next;}
		chomp;
		my @spl = split(/\t/);
		$spl[8] =~ m/^ID=\d+_(\d+)/;
		my $k = ">".$spl[0]."_$1";
		#print $k;
		$ret{$k}=$_;
	}
	#die;
	close I;
	return \%ret;
}


sub emptyAssGrpsObj($){
my ($hr) = @_;
my %AsGrps = %{$hr};
foreach my $k (keys %AsGrps){
	$AsGrps{$k}{CntAss} = 0;
	$AsGrps{$k}{CntMap} = 0;
	#dependencies
	$AsGrps{$k}{ITSDeps} = "";
	$AsGrps{$k}{MapDeps} = "";
	$AsGrps{$k}{DiamDeps} = "";
	$AsGrps{$k}{SeqClnDeps} = "";
	#print $AsGrps{$k}{CntAimAss} ."\n";
	#remove tmp dirs
	$AsGrps{$k}{ClSeqsRm} = "";
	#copy to final folder
	@{$AsGrps{$k}{MapCopies}} = ();
	@{$AsGrps{$k}{AssCopies}} = ();
	$AsGrps{$k}{MapCopiesNoDel} = [];
	#filteredSequenceFiles for assembly
	@{$AsGrps{$k}{FilterSeq1}} = (); @{$AsGrps{$k}{FilterSeq2}} = (); @{$AsGrps{$k}{FilterSeqS}} = ();
	@{$AsGrps{$k}{RawSeq1}} = (); @{$AsGrps{$k}{RawSeq2}} = (); @{$AsGrps{$k}{Libs}} = ();
	$AsGrps{$k}{PostAssemblCmd} = "";
	$AsGrps{$k}{PostClnCmd} = "";
	$AsGrps{$k}{AssemblSmplDirs} = "";
	$AsGrps{$k}{scndMapping} = "";
	$AsGrps{$k}{ClSeqsRm} = "";
	#print $k."\n";;
}
$AsGrps{global}{DiamCln} = "";
$AsGrps{global}{DiamDeps} = "";
#die;
return \%AsGrps;
}

sub readMapS($ $){
	my ($inF,$folderStrClassical) = @_;
	my @spl = split /,/,$inF;
	my %ret; my %agbp;
	my @outDirs;
	my $hr1 = \%ret;my $hr2 = \%agbp;  my $cnt = -1;
	foreach my $map (@spl){
		($hr1,$hr2) = readMap($map,$cnt,$hr1,$hr2,$folderStrClassical);
		%ret = %{$hr1};
		$cnt = $ret{totSmpls};
		push(@outDirs,$ret{outDir});
		#print $cnt."\n";
	}
	#%ret = %{$hr1};
	#print keys %ret;
	$ret{outDir} = join(",",@outDirs);
	return (\%ret,$hr2);
}
sub readMap{
	my $inF = $_[0];
	my $Scnt = defined $_[1] ? $_[1] : 0;
	my %ret = defined $_[2] ? %{$_[2]} : (); 
	my %agBP = defined $_[3] ? %{$_[3]} : ();
	my $folderStrClassical = defined $_[4] ? $_[4] : 0;
#die "$folderStrClassical\n";
	my @order = exists $ret{smpl_order} ? @{$ret{smpl_order}} : ();
	my $dirCol = -1; my $smplCol = 0; my $rLenCol = -1; my $SeqTech = -1; my $SmplPrefixCol = -1;
	my $AssGroupCol = -1; my $EstCovCol = -1; my $MapGroupCol = -1; my $SupRdsCol = -1;my $ExcludeAssemble = -1;
	#some global params
	my $dir2dirs = ""; #dir on file system, where all dirs specified in map can be found (enables different indirs with different mapping files)
	my $dir2out = "";
	my $baseID = ""; my $mocatFiltPath = "";
	my $inDirSet = 0;
	my $cnt = -1;
	my @dir2dirsA;
	my %trackMGs; my %trackAGs; #hashes to track the last (final) sample in each mapgroup.. important to know this to check if assembly / mapping is done 
	my %memberMGs ; my %memberAGs ;
	#print $inF."\n";
	open I,"<$inF" or die "Can't open $inF\n";
	#AssGrps
	while (<I>){
		$cnt++; chomp;
		if (m/^#/ && $cnt > 0 ){#check for ssome global parameters
			if (m/^#DirPath\s+(\S+).*$/){$dir2dirs = $1; $dir2dirs.="/" unless ($dir2dirs=~m/\/$/); push(@dir2dirsA,$dir2dirs);}
			if (m/^#OutPath\s(\S+)/){$dir2out = $1; $inDirSet=0;$dir2out .= "/" if ($dir2out !~ m/\/$/);}
			if (m/^#RunID\s(\S+)/){$baseID = $1; $inDirSet=0;}
			if (m/^#mocatFiltPath\s(\S+)/){$mocatFiltPath = $1;}
			if (!$inDirSet && $dir2out ne "" && $baseID ne ""){$dir2out.=$baseID unless ($dir2out =~ m/$baseID[\/]*$/); $inDirSet =1;}
			$dir2out .= "/" if ($dir2out !~ m/\/$/);
			next;
		}
		my @spl = split(/\t/);
		if ($cnt == 0){
			#die "@spl\n";
			$dirCol = first_index { /Path/i } @spl;
			$smplCol = first_index { /#SmplID/ } @spl;
			$SeqTech = first_index { /SeqTech/ } @spl;
			$SmplPrefixCol = first_index { /SmplPrefix/i } @spl;
			$rLenCol = first_index { /ReadLength/ } @spl;
			$AssGroupCol = first_index { /AssmblGrps/ } @spl;
			$EstCovCol = first_index { /EstCoverage/ } @spl;
			$MapGroupCol = first_index { /MapGrps/ } @spl;
			$SupRdsCol = first_index { /SupportReads/ } @spl;
			$ExcludeAssemble = first_index { /ExcludeAssembly/ } @spl;
			
			die "Only \"Path\" or \"SmplPrefix\" can be defined in mapping file. Both is not supported.\n" if ($dirCol != -1 && $SmplPrefixCol != -1);
			next;} #maybe later check for col labels etc
		$Scnt++;
		#die $spl[0]." ".$spl[1]."\n";
		die "inPath has to be set in mapping file!\n" if ($dir2dirs eq "");
		my $curSmp = $spl[$smplCol];
		my $altCurSmp = "";
		#print $curSmp." ";
		die "Double sample ID $curSmp\n" if (exists $ret{$curSmp});
		my $cdir = ""; 
		$cdir = $spl[$dirCol] if ($dirCol >= 0); 
		my $cdir2= $cdir;
		$ret{$curSmp}{dir} = $cdir;#this one should stay without a tag
		$ret{$curSmp}{rddir} = $dir2dirs.$cdir;
		$ret{$curSmp}{rddir} .="/" unless ($ret{$curSmp}{rddir} =~ m/\/$/);
		#die "$ret{$curSmp}{rddir} $dirCol $cdir $curSmp\n $smplCol $dirCol\n";
		if ($SmplPrefixCol>=0){$cdir2 = $spl[$SmplPrefixCol];; $ret{$curSmp}{prefix} = $cdir2;} else {$ret{$curSmp}{prefix} = "";}
		$cdir2.="/" unless ($cdir2 =~ m/\/$/);
		if ($folderStrClassical){
			$ret{$curSmp}{wrdir} = $dir2out.$cdir2;
		} else {
			$ret{$curSmp}{wrdir} = $dir2out.$curSmp."/";
		}
		#die "$ret{$curSmp}{wrdir}\n";
		$ret{$curSmp}{SmplID} = $curSmp;
		$ret{$curSmp}{mapFinSmpl} = $curSmp;
		$ret{$curSmp}{assFinSmpl} = $curSmp;
		if ($SeqTech >= 0) {$ret{$curSmp}{SeqTech} = $spl[$SeqTech]; } else {$ret{$curSmp}{SeqTech} = "";}
		#die "$SupRdsCol\t@spl\n$spl[$SupRdsCol]\n";
		if ($SupRdsCol >= 0 && $SupRdsCol < @spl) { if(length($spl[$SupRdsCol]) > 0 && $spl[$SupRdsCol] !~ m/\/$/) {$spl[$SupRdsCol].="/";} $ret{$curSmp}{SupportReads} = $spl[$SupRdsCol]; } else {$ret{$curSmp}{SupportReads} = "";}
		if ($rLenCol >= 0){$ret{$curSmp}{readLength} = $spl[$rLenCol];} else {$ret{$curSmp}{readLength} = 0;}
		if ($ExcludeAssemble >= 0){$ret{$curSmp}{ExcludeAssem} = $spl[$ExcludeAssemble];} else {$ret{$curSmp}{ExcludeAssem} = 0;}
		
		if ($EstCovCol >= 0){$ret{$curSmp}{DoEstCoverage} = $spl[$EstCovCol];} else {$ret{$curSmp}{DoEstCoverage} = 0;}
		if ($AssGroupCol >= 0){
			my $curAG = $spl[$AssGroupCol];
			$ret{$curSmp}{AssGroup} = $curAG ;
			if (!exists($agBP{$curAG}{CntAimAss})){$agBP{$curAG}{CntAimAss}=0;}
			$agBP{$curAG}{CntAimAss}++;
			$altCurSmp = $curSmp."M".$agBP{$curAG}{CntAimAss};
			$memberAGs{$curAG} = [] unless (exists $memberAGs{$curAG});
			$trackAGs{$curAG} = $curSmp;  push(@{$memberAGs{$curAG}},$curSmp);
			$agBP{$curAG}{prodRun} = "";
			#print $agBP{$spl[$AssGroupCol]}{CntAimAss}. " :$spl[$AssGroupCol]\n" ;
		} else {$ret{$curSmp}{AssGroup} = $Scnt; $agBP{$Scnt}{CntAimAss}=0;$agBP{$Scnt}{prodRun} = "";}
		if ($MapGroupCol >= 0){
			my $curM = "M_".$spl[$MapGroupCol];
			$ret{$curSmp}{MapGroup} = $curM;
			if (!exists($agBP{$curM}{CntAimMap})){$agBP{$curM}{CntAimMap}=0;}
			$agBP{$curM}{CntAimMap}++;
			$memberMGs{$curM} = [] unless (exists $memberMGs{$curM});
			$trackMGs{$curM} = $curSmp; push(@{$memberMGs{$curM}},$curSmp);
			#print $agBP{$spl[$MapGroupCol]}{CntAimMap}. " :$spl[$MapGroupCol]\n" ;
		} else {$ret{$curSmp}{MapGroup} = $Scnt; $agBP{$Scnt}{CntAimMap}=0;}
		if ($altCurSmp ne ""){$ret{altNms}{$altCurSmp} = $curSmp;}
		push(@order,$curSmp);
		#print $spl[0]."\n";
	}
	
	#insert final sample destination for all AssGroups and MapGroups
	foreach my $k (keys %memberMGs){
		foreach (@{$memberMGs{$k}}){  $ret{$_}{mapFinSmpl} = $trackMGs{$k};  }
	}
	foreach my $k (keys %memberAGs){
		foreach (@{$memberAGs{$k}}){  $ret{$_}{assFinSmpl} = $trackAGs{$k};  }
	}
	
	#die();
	close I;
	$ret{totSmpls} = $Scnt;
	$ret{smpl_order} = \@order;
	$ret{inDir} = join(",",@dir2dirsA) if ($dir2dirs ne "");
	$ret{outDir} = $dir2out if ($dir2out ne "");
	$ret{baseID} = $baseID if ($baseID ne "");
	$ret{mocatFiltPath} = $mocatFiltPath if ($mocatFiltPath ne "");
	#@order = keys %agBP;die "@order\n";
	my $asGrpHr = emptyAssGrpsObj(\%agBP);
	return (\%ret,$asGrpHr);
}

sub randStr($){
	my ($len) = @_;
	my @letters=('A'..'Z','a'..'z',1..9);
	my @letters2=('A'..'Z','a'..'z');
	my $total=scalar(@letters);
	my $newletter ="";
	$newletter = $letters2[rand scalar(@letters)];
	for (my $i=1;$i<$len;$i++){
		$newletter .= $letters[rand $total];
	}
	return $newletter;
}



sub findQsubSys($){
	my ($iniVal) = @_;
	#my $iniVal = "lsf";
	if ($iniVal ne ""){
		$iniVal = lc $iniVal; $iniVal = "lsf" if ($iniVal eq "bsub");
		$iniVal = "sge" if ($iniVal eq "qsub");
		$iniVal = "slurm" if ($iniVal eq "sbatch");
	} else {
		$iniVal = "lsf";
		my $bpath = `which bsub  2>/dev/null`;chomp $bpath;my $bpresent=0; 
		$bpresent=1 if ($bpath !~ m/\n/ && -e $bpath);
		my $spath = `which sbatch  2>/dev/null`;chomp $spath;my $spresent=0; 
		$spresent=1 if ($spath !~ m/\n/ && -e $spath);
		my $qpath = `which qsub  2>/dev/null`; chomp $qpath;
		my $qpresent=0; $qpresent=1 if ($qpath !~ m/\n/ && -e $qpath);
		#print "$qpath\n";
		if ($spresent ){#slurm gets preference
			$iniVal="slurm";
		}elsif (!$bpresent && $qpresent){
			$iniVal = "sge";
		}elsif (!$qpresent && !$bpresent && !$spresent){
			print "Warning: No queing system found (sbatch / qsub / bsub command)\nUsing LSF (bsub), though this will likely cause errors\n";
			
		}
	}
	print "Using qsubsystem: $iniVal\n";
	#die;
	return $iniVal;
}
sub emptyQsubOpt{
	my ($doSubm) = $_[0];
	my $locChkStr = $_[1];
	my $qmode;
	
	if (@_ > 2){$qmode = $_[2];}
	else {$qmode = findQsubSys("");}
	die "qsub system mode has to be \'lsf\', \'slurm\' or \'sge\'!\n" if ($qmode ne "lsf" &&$qmode ne "slurm" && $qmode ne "sge");
	my $MFdir = getProgPaths("MFLRDir");
	my $xtraNodeCmds = getProgPaths("subXtraCmd",0);
	$xtraNodeCmds = "" unless (defined $xtraNodeCmds);
	my $longQ = "medium_priority"; my $shortQ = "medium_priority";
	if ($qmode eq "slurm"){$shortQ = "1day"; $longQ="1month";}
	my %ret = (
		rTag => randStr(3),
		doSubmit => $doSubm,
		LocationCheckStrg => $locChkStr,
		doSync => 0,
		longQueue => $longQ,
		shortQueue => $shortQ,
		useLongQueue => 0,
		qsubPEenv => getProgPaths("qsubPEenv"),
		perl5lib => "$MFdir:\$PERL5LIB",
		cpplib => "/g/bork3/home/hildebra/env/env/miniconda/bin/gcc/",
		tmpSpace => "30G",
		xtraNodeCmds => $xtraNodeCmds,
		qmode => $qmode,
		#LOG => undef,
	);
	#die "$MFdir\n";
	return \%ret;
}
sub qsubSystem($ $ $ $ $ $ $ $ $ $){
	#args: 1[file to save bash & error & output] 2[actual bash cmd] 3[cores reserved for job]
	# 4["1G": Ram usage per core in GB] 5[0/1: synchronous execution] 6[name of job] 
	# 7[name of job dependencies, separated by ";"]
	# 8[0/1: excute in cwd?] 9[0/1: return qsub cmd or submit job to cluster]
	# Falk Hildebrand, may 2015
	my ($tmpsh,$cmd,$ncores,$mem,$jname,$waitJID,$cwd,$immSubm, $restrHostsAR, $optHR) = @_;
	#$doSync, 5th arg
	#14,12G
	#die $tmpsh."\n";
	#my $jname = $tmpsh;
	#$jname =~ s/.*\///g;$jname =~ s/\.sh$//g;
	#\n#\$ -N $tmpsh
	return("") if ($cmd eq "");
	my $LSF = 0;
	my $qbin = "qsub";
	my $xtra = "";
	my $memory = $mem;
	my $rTag = $optHR->{rTag};
	my $qmode = $optHR->{qmode};
	my $xtraNodeCmds = $optHR->{xtraNodeCmds};
	my @restrHosts = @{$restrHostsAR};
	#different format for bsub and slurm
	if ($memory =~ s/G$//){$memory *= 1024 * $ncores;};
	my $tmpSpace = $optHR->{tmpSpace};
	if ($tmpSpace =~ s/G$//){$tmpSpace *= 1024 ;};
	#die ($memory."\n");
	my $queues = "\"".$optHR->{shortQueue}."\"";#"\"medium_priority\"";
	$queues = "\"".$optHR->{longQueue}."\"" if ($optHR->{useLongQueue} ==1);#"\"medium_priority\"";
	if ($memory > 250001){$queues = "\"scb\"";}
	open O,">",$tmpsh or die "Can't open qsub bash script $tmpsh\n";

	#die "$cmd\n";
	#print "$memory   $queues\n";
	#if (`hostname` !~ m/submaster/){
	if ($qmode eq "slurm"){$LSF = 2;$qbin="sbatch";
		if ($memory > 250001){$queues = "\"bigmem\"";}
		system "rm -f $tmpsh.otxt $tmpsh.etxt";
		print O "#!/bin/bash\n#SBATCH -N 1\n#SBATCH --cpus-per-task  $ncores\n#SBATCH -o $tmpsh.otxt\n"; #\n#SBATCH -n  $ncores
		print O "#SBATCH --tmp=$tmpSpace\n#SBATCH -e $tmpsh.etxt\n#SBATCH --mem $memory\n#\$ --export ALL\n";
		print O "#SBATCH -p $queues\n";
		print O "#\$ -S /bin/bash\n#\$ -v LD_LIBRARY_PATH=".$optHR->{cpplib}."\n#\$ -v TMPDIR=/dev/shm\n";
		print O "#\$ -v PERL5LIB=".$optHR->{perl5lib}."\n";
	} elsif ($qmode eq "sge"){
		system "rm -f $tmpsh.otxt $tmpsh.etxt";
		print O "#!/bin/bash\n#\$ -S /bin/bash\n#\$ -cwd\n#\$ -pe ".$optHR->{qsubPEenv}." $ncores\n#\$ -o $tmpsh.otxt\n#\$ -e $tmpsh.etxt\n#\$ -l h_vmem=$mem\n";
		print O "#\$ -v LD_LIBRARY_PATH=".$optHR->{cpplib}."\n#\$ -v TMPDIR=/dev/shm\n";
		print O "#\$ -v PERL5LIB=".$optHR->{perl5lib}."\n";
	} else {$LSF = 1;$qbin="bsub";
		print O "#!/bin/bash\n";
		print O "export LD_LIBRARY_PATH=/g/bork3/home/hildebra/env/env/miniconda/lib/:/g/bork3/home/hildebra/env/zlib-1.2.8/:/g/bork8/costea/boost_1_53_0/:/shared/ibm/platform_lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/g/bork3/x86_64/lib64:/g/bork3/x86_64/lib:\${LD_LIBRARY_PATH}\n\n";
		#print O "export LD_LIBRARY_PATH=/g/bork3/home/hildebra/env/zlib-1.2.8:/g/bork3/x86_64/lib64:/lib:/lib64:/usr/lib64:\${LD_LIBRARY_PATH}\n\n";#:/g/software/linux/pack/python-2.7/lib/\nexport PATH=/g/bork3/home/zeller/py-virtualenvs/py2.7_bio1/bin/:\${PATH}\n\n";
		##BSUB -n $ncores\n#BSUB -o $tmpsh.otxt\n#BSUB -e $tmpsh.etxt\n#BSUB -M $mem\n#\$ -v LD_LIBRARY_PATH=/g/bork3/x86_64/lib64:/lib:/lib64:/usr/lib64\n#\$ -v TMPDIR=/dev/shm\n#BSUB -q medium_priority\n";
		if ( @restrHosts > 0){
			$xtra .= " -m \"".join(" ",@restrHosts)."\" ";
			$queues = "\"medium_priority scb\"";
		}
		$xtra .= "-n $ncores -oo $tmpsh.otxt -eo $tmpsh.etxt -q $queues -M $memory -R \"select[(mem>=$memory)] rusage[tmp=$tmpSpace] span[hosts=1]\" -R \"rusage[mem=$memory]\" "; #
	}
	#set abortion on program fails
	print O "set -e\n";
	print O "ulimit -c 0;\n";
	#any xtra commands (like module load perl?)
	print O "$xtraNodeCmds\n";
	#prevent core dump files
	#file location check availability
	print O $optHR->{LocationCheckStrg};

	print O $cmd."\n";
	close O;
	my $depSet=0;
	my @jspl = split(";",$waitJID); @jspl = grep /\S/, @jspl;
	if ($LSF==2){#slurm
		if ($optHR->{doSync} == 1){$qbin = "srun";}
		if (length($waitJID) >3) {
			#my $idsstr= "squeue --name=".join(":",@jspl)." | awk '{if(NR>1)print}' | sed 's/\\s\\+/ /g' | cut -f2 -d' ' ";
			#$idsstr = `$idsstr`;
			#die $idsstr."\n";
			#@jspl = split(/\n/,$idsstr);
			#die "@jspl\n";
			#$depSet=1;
			for (@jspl) {s/$rTag//;}
			$xtra .= "--dependency=afterok:".join(":",@jspl)." " if (@jspl > 0);
			#print $xtra;
		}
		if ($jname ne ""){$xtra.="-J $rTag$jname ";}
		#finish up
		$xtra .= " ";
	} elsif ($LSF==1){ #bsub #-M memLimit; -q queueName;  -m "host_name[@cluster_name]; -n minProcessors; 
		if ($optHR->{doSync} == 1){$xtra.="-K ";}
		if ($jname ne ""){$xtra.="-J $rTag$jname ";}
		if (length($waitJID) >3) {
			my @jspl = split(";",$waitJID);
			#remove empty elements
			@jspl = grep /\S/, @jspl;
			if (@jspl > 0 ){
				$waitJID = join(") && done(",@jspl);
				$xtra.="-w \"done($waitJID)\" ";
			}
		}
		$tmpsh = " < ".$tmpsh;
	} else{ #qsub
		if ($optHR->{doSync} == 1){$xtra.="-sync y ";}
		if ($jname ne ""){$xtra.="-N $rTag$jname ";}
		if (length($waitJID) >3) {
			if (@jspl > 0 ){$xtra.="-hold_jid ".join(",",@jspl) ." ";}
		}
			#$waitJID =~ s/;/,/g;$xtra.="-hold_jid $waitJID ";}
	}
	if ($cwd ne ""){if ($LSF==2){$xtra .= "--workdir $cwd";} elsif ($LSF==1) {$xtra.="-cwd $cwd"; }  else {$xtra.="-wd $cwd";} }
	my $qcm = "$qbin $xtra $tmpsh \n";
	my $LOGhandle = "";
	if (exists $optHR->{LOG}){ $LOGhandle = $optHR->{LOG};}
	#print("$qcm\n\n");
	#if (@restrHosts > 0){die $qcm;}
	if ($optHR->{doSubmit} != 0 && $immSubm){
		print $LOGhandle $qcm."\n" unless ($LOGhandle eq "" || !defined($LOGhandle) );
		print "SUB:$jname\t";
		my $ret = `$qcm`; 
		if ($LSF == 2){#slurm get jobid
			chomp $ret; $ret =~ m/(\d+)$/; #$ret = $1;
			$jname=$1;
		}
		#die "$qcm\n$cmd\n$ret XX\n";
		#if ($depSet){die "$qcm";}
	}

	if ($immSubm==0){
		return "$rTag$jname",$qcm;
	}
	return "$rTag$jname",$qcm;
}

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub systemW($){
	my ($cmddd) = @_;
	if ($cmddd eq ""){return;}
	my $stat= system($cmddd);
	if ($stat){
		print "system call \n$cmddd\nfailed.\n";
		exit(33);
	}
}

sub readTabByKey{
	my ($inF) = @_;
	my %ret;
	my ($I,$OK) = gzipopen($inF,"tab file");
	my $maxTabs = 0;
	while (my $l = <$I>){
		chomp $l; my @tmp = split /\t/,$l;
		$ret {$tmp[0]} = $tmp[1];
		if (@tmp > $maxTabs){$maxTabs = @tmp;}
	}
	close $I;
	if ($maxTabs > 2){print "Warning in Tab reader: more than 2 columns were present ($maxTabs)\n";}
	return %ret;
}

sub writeFasta{
	my ($hr,$of) = ($_[0],$_[1]);
	my $maxFs = -1;
	$maxFs = $_[2] if (@_ > 2);
	my %FA = %{$hr};
	my $cnt=0;
	open O,">$of" or die "can't open out fasta $of\n";
	foreach my $k (keys %FA){
		$cnt++; 
		if ($k =~ m/^>/){
			print O $k."\n".$FA{$k}."\n";
		} else {
			print O ">".$k."\n".$FA{$k}."\n";
		}
		last if ($maxFs > 0 && $cnt > $maxFs);
	}
	close O;
}

sub convertMSA2NXS{
	my $filename = $_[0];
	my $outF = "";
	$outF = $_[1] if (@_>1);
	my $hr = readFasta($filename);
	my %FNAs = %{$hr};
	my @kk = keys %FNAs;
	my $ostr="";
	my $numtaxa = scalar(@kk);
	my $maxlength = length($FNAs{$kk[0]}); #all seqs should be same length in MSA format Format datatype=dna missing=? gap=-;
	$ostr = "#NEXUS\nBegin data;\nDimensions ntax=$numtaxa nchar=$maxlength;\nFormat datatype=dna missing=? gap=-;\nMatrix\n";

	foreach my $k (@kk) {
		my $len=length$FNAs{$k};
		die "Error nexus format conversion: $len != $maxlength in $k\n" if ($len != $maxlength);
		#if ($len<$maxlength) { my $add=$maxlength-$len; for (my $j=0; $j<$add; $j++) {$seqs[$i]=$seqs[$i].'-';}}
		$ostr.= "\n$k\t$FNAs{$k}";
	}

	$ostr .= "\n;\nend;";
	if ($outF ne ""){
		open O,">$outF" or die "Can't open nxs out file $outF";
		print O $ostr;
		close O;
	}
	return $ostr;
}


 