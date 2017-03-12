#!/usr/bin/env perl
#parseBlastFunct.pl /g/bork5/hildebra/data/Ach_proteins/res/DiaAssignment.txt [DB] [mode]
#mode == 0: just parse; ==2 check for all files present
use warnings;
use strict;
use FileHandle;
use Mods::TamocFunc qw(sortgzblast uniq);
use Mods::GenoMetaAss qw(gzipwrite gzipopen);

sub main;
sub readGene2COG; sub readNogKingdom;
sub readCOGdef;
sub readGene2KO;sub readKeggTax;
sub check_files; sub remove_file;
sub writeAllTable;
sub readMohTax;
sub combineBlasts;
sub bestBlHit;



die "no input args!\n" if (@ARGV == 0 );
my $blInf = $ARGV[0];
my $mode = 0;
my $DBmode = "NOG";
if (@ARGV > 1){
	$DBmode = $ARGV[1]; #0=normal 2=file check 4=remove outputfile
}

my $quCovFrac = 0; #how much of the query needs to be covered?
my $noHardCatCheck = 0; #select for the hit with KO assignment rather than the real best hit (w/o KO assignment)

#$DBmode = uc $DBmode;
my $minBLE= 1e-7;my $minScore=0;
my $minAlLen = 20; #my $minPidGlobal = 40;


if (@ARGV > 2){
	$minBLE = $ARGV[2];
}
if (@ARGV > 3){$mode = $ARGV[3];}
#die $blInf;

my $tmpD="";
if (@ARGV >= 6){$tmpD = $ARGV[6];}
#die "$tmpD\n";
#work out if sorted etc
$blInf = sortgzblast($blInf,$tmpD) unless ($mode == 2 || $mode == 4);

my %DBlen;
if (@ARGV > 4 ){#read the length of DB proteins
	my $lengthF = $ARGV[4];
	my @spl = split /,/,$lengthF;
	foreach my $lF2 (@spl){
		open I,"<$lF2" or die "length file $lF2 no present!\n";
		while (my $l = <I>){
			chomp $l;
			my @spl = split /\t/,$l;
			$DBlen{$spl[0]} = int $spl[1];
		}
		close I;
	}
	print "read length DB\n";
}


#die "$blInf";

die "Too few args\n" if (@ARGV<5 && $mode > 1);
my $DButil ;
my $bl2dbF ;
my $cogDefF ;
my $NOGtaxf ;
my $KEGGlink ;
my $KEGGtaxDb ;
if (@ARGV >= 5){
	$DButil = $ARGV[5];
	 $bl2dbF = "$DButil/NOG.members.tsv";
	 $cogDefF = "$DButil/NOG.annotations.tsv";
	 $NOGtaxf = "$DButil/all_species_data.txt";
	 $KEGGlink = "$DButil/genes_ko.list";
	$KEGGtaxDb = "$DButil/kegg.tax.list";

}

$blInf =~ m/(.*)\/([^\/]+$)/;
my $inP = $1; my $inF = $2;
#die "$inP\t$inF\n";

my @blRes;
my @kgdOpts = qw (0 1 2 3);
my @kgdName = ("Bacteria","Archaea","Eukaryota","Fungi");
my @kgdNameShrt = ("BAC","ARC","EUK","FNG");
#which best blast decide algo to use..
my $scomode = 0; my $emode = 1;
my @normMethods = ("cnt","GLN");
my $czyEuk =2; my $czyFungi=5;
my %czyTax;
my @aminBLE = split /,/,$minBLE;
#my @aminPID = (40,45,50,55,60,65,70,75,80,85,90,95);
#@aminBLE = ("1e-9") x scalar(@aminPID);
my @aminPID = ("30") x scalar(@aminBLE);

my %NOGkingd ; my %KEGGtax;
my %COGdef ;my %g2COG ; my %c2CAT ;
my $tabCats = 0; my $cazyDB = 0; my $ACLdb = 0; my $KGMmode=0;
if ($mode == 2 || $mode == 4){ #scan for all reads finished
} elsif ( $DBmode eq "KGB" || $DBmode eq "KGE" || $DBmode eq "KGM"){
	print "Reading KEGG DBs..\n";
	my ($hr1) = readGene2KO($KEGGlink);
	%c2CAT = %{$hr1};
	#my @tmp = keys %c2CAT; die "$tmp[0] $tmp[1] $tmp[2]\n";
	@kgdOpts = qw (3); $tabCats = 2;
	if ($DBmode eq "KGB"){ @kgdName = ("","","","Bacteria");@kgdNameShrt = ("","","","BAC") ;}
	if ($DBmode eq "KGE"){ @kgdName = ("","","","Eukaryotes");@kgdNameShrt = ("","","","EUK") ;}
	if ($DBmode eq "KGM"){ @kgdName = ("Prokayotes","Eukaryotes","","Unclassified");@kgdNameShrt = ("BAC","EUK","","UNC") ; 
		@kgdOpts = qw (0 1 3); $KGMmode=1; 
	}
	$hr1 = readKeggTax($KEGGtaxDb);
	%KEGGtax = %{$hr1};
} elsif ( $DBmode eq "NOG"){
	print "Reading NOG DBs..\n";
	my ($hr1,$hr2) = readGene2COG($bl2dbF);
	%g2COG = %{$hr1}; %c2CAT = %{$hr2};
	$hr1 = readCOGdef($cogDefF);
	%COGdef = %{$hr1}; $tabCats = 1;
	%NOGkingd = readNogKingdom($NOGtaxf);
} elsif ($DBmode eq "CZy"){
	@kgdOpts = qw (0 1 2 3 4 5);
	my $hr = readMohTax($DButil."/MohCzy.tax");
	%czyTax = %{$hr};
	$czyEuk =2; $czyFungi=5;
	@kgdName = ("Bacteria","Archaea","Eukaryota","unclassified","CBM","Fungi");
	@kgdNameShrt = ("BAC","ARC","EUK","UNC","CBM","FNG");
	$cazyDB = 1;
	print "CAZy database\n";
} elsif ($DBmode eq "ACL"){
	@kgdOpts = qw (3);
	$czyEuk =2; $czyFungi=5;
	@kgdName = ("Plasmid","Prophage","Virus","unclassified");
	@kgdNameShrt = ("PLA","PRO","VIR","UNC");
	$ACLdb = 1;
	print "aclame database\n";
} else{
	@kgdOpts = qw (3);
	@kgdName = ("","","","All");
	@kgdNameShrt = ("","","","ALL");
	print "Non - NOG DB\n";
}


#read DB file
my $readLim = 1000000; my $lcnt=0;my $lastQ=""; my $stopInMiddle=0; my $reportGeneCat=0;
my %COGhits; my %CAThits; my %GENEhits; my $O2;
for (my $i=0; $i<@aminBLE ; $i++){$COGhits{$i} = {}; $CAThits{$i} = {}; $GENEhits{$i} = {};}
if ($mode == 0 || $mode==1){ #mode1 = write gene assignment
	#die "$blInf\n";
	$reportGeneCat = 1 if ($mode==1);
	#insert here KEGG BAC ? EUK fix
	if ($KGMmode && !-e $blInf){
		#zcat dia.KGB.blast.gz dia.KGE.blast.gz | sort | gzip > dia.KGM.blast.gz
		my $KGBf = $blInf; my $KGEf = $blInf; 
		$KGBf =~ s/KGM/KGB/; $KGEf =~ s/KGM/KGE/;
		my $cmd = "zcat $KGBf $KGEf | sort | gzip > $blInf\n";
		system $cmd."\n";
	}
	
	my ($I,$OK) = gzipopen($blInf,"diamond output file",1); 
	my $OK2;
	($O2,$OK2) = gzipwrite($blInf,"diamond output file",1); 
	my @splSave;
	
	while (1){
		my %wordv1; my %wordv2;
		while (my $line = <$I>){
			chomp $line; 
			my @spl = split (/\t/,$line);
			my $query = $spl [0];
			$query =~ s/\/\d$//;
			
#			$lcnt++;
			
			if ($lastQ eq ""){$lastQ = $query;
			} elsif ($lastQ ne $query){
				$stopInMiddle=1;@splSave = @spl;last;
			}
			if (@spl > 10 ){ if ( $spl[0] =~ m/2$/){
				$wordv2{$spl [1]} = \@spl;
				} else {$wordv1{ $spl [1]} = \@spl;}
			}

#			push(@blRes,\@spl);
		}
		#1st combine scores
		my $whX = combineBlasts(\%wordv1,\%wordv2);
		#2nd: sort these, find best hit
		#my $arBhit = bestBlHit($whX); #doesn't need this, this is done in the main routine
		@blRes = values %{$whX};

		#die @blRes."\n";
		#print "Read Assignments..\n";
		#die "@aminBLE\n";
		for (my $i=0; $i<@aminBLE ; $i++){
			#my ($hr1,$hr2,$hr3) = 
			main($aminBLE[$i],$aminPID[$i],$i);#$st1{$i}, $CAThits{$i},$GENEhits{$i});
			#$st1{$i} = $hr1; $CAThits{$i} = $hr2; $st3{$i} = $hr3;
			#print "..";
		}
		undef @blRes;
		push(@blRes,\@splSave);
		if ($stopInMiddle==0){last;}
		$lcnt=0;$lastQ="";$stopInMiddle=0;
	}
	close $I;
	close $O2 if ($reportGeneCat);
	print"writing Tables\n";
	
	
	for (my $i=0; $i<@aminBLE ; $i++){
		my $pathXtra = "/CNT_".$aminBLE[$i]."_".$aminPID[$i]."/";
		my $outPath = $inP.$pathXtra;
		writeAllTable($DBmode."parse",$outPath,$aminBLE[$i],$aminPID[$i],$i);#$st1{$i}, $CAThits{$i},$st3{$i});
	}
	system "touch $blInf.stone";
	print "all done\n";
	exit(0);

} elsif ($mode==2) {
	for (my $i=0; $i<@aminBLE ; $i++){
		my $pathXtra = "/CNT_".$aminBLE[$i]."/";
		unless(check_files($DBmode."parse",$inP.$pathXtra,$aminBLE[$i])){exit(3);}
	}
} elsif ($mode==4) {

	for (my $i=0; $i<@aminBLE ; $i++){
		my $pathXtra = "/CNT_".$aminBLE[$i]."/";
		remove_file($DBmode."parse",$inP.$pathXtra,$aminBLE[$i]);
	}
} else {die"unkown run mode!!!\n";}

if ($mode ==2){
}

print "all done\n";
exit(0);


sub remove_file{
	my ($outF,$outD,$minBLE) = @_;
	my $out = $outD."/".$outF;
	system "rm -f $out*" if (-d $outD);
}
sub check_files{
	my ($outF,$outD,$minBLE) = @_;
	my $out = $outD."/".$outF;
	my @srchFls = (#"$out.$DBmode.GENE2$DBmode", #"$out.GENE2$DBmode"."aCAT"
	"$out.$DBmode.$kgdNameShrt[0].gene.cnts","$out.$DBmode.$kgdNameShrt[0].cat.cnts");
	push (@srchFls,"$out.$kgdNameShrt[0].CATcnts") if ( $DBmode eq "NOG");
	foreach (@srchFls) { unless (-e $_.".gz"){print "Fail $_\n";return 0;}}
	return 1;
}

sub writeAllTable(){
	my ($outF,$outD,$minBLE,$minPID,$pos) = @_;#,$hr1,$hr2,$hr3) = @_;
	system "mkdir -p $outD" ;#or die "Failed to create out dir $outD\n";
	my $out = $outD."/".$outF;
	system "rm -f $out*";
	my %COGabundance=%{$COGhits{$pos}}; my %CATabundance=%{$CAThits{$pos}}; my %funHit=%{$GENEhits{$pos}};
	#my %COGabundance=%{$hr1}; my %CATabundance=%{$hr2}; my %funHit=%{$hr3};
		foreach my $normMethod (@normMethods){

		foreach my $y (@kgdOpts){
			my @kk; my $O = FileHandle->new;
			if ($tabCats==0){ #MOG CZy etc
				#don't write this for NOGs any longer..
				@kk = sort { $funHit{$normMethod}{$y}{$b} <=> $funHit{$normMethod}{$y}{$a} } keys(%{$funHit{$normMethod}{$y}});
				$O = gzipwrite("$out.$DBmode.$kgdNameShrt[$y].$normMethod.gene.cnts","Dia gene $y Counts");
				#open O,">$out.$DBmode.gene.cnts";
				foreach my $k (@kk){print $O $k."\t".$funHit{$normMethod}{$y}{$k}."\n";}
				close $O;
			}
				#print COGs
			if ($tabCats != 2){
				@kk = sort { $COGabundance{$normMethod}{$y}{$b} <=> $COGabundance{$normMethod}{$y}{$a} } keys(%{$COGabundance{$normMethod}{$y}});
				$O = gzipwrite("$out.$DBmode.$kgdNameShrt[$y].$normMethod.cat.cnts","Dia cat $y Counts");
				#open O,">$out.$DBmode.cat.cnts"; 
				my $cogCnts=0;
				foreach my $k (@kk){print $O $k."\t".$COGabundance{$normMethod}{$y}{$k}."\n";$cogCnts++;}
				close $O;
				print "Total of $cogCnts categories in $kgdName[$y].\n";
			}
			if ($tabCats){ #NOG || KEGG
				#print CATs
				@kk = sort { $CATabundance{$normMethod}{$y}{$b} <=> $CATabundance{$normMethod}{$y}{$a} } keys(%{$CATabundance{$normMethod}{$y}});
				#print @kk."\n";
				$O = gzipwrite("$out.$kgdNameShrt[$y].$normMethod.CATcnts","Dia NOG cat $y Counts");
				#open O,">$out.CATcnts";
				my $CATcnt =0;
				foreach my $k (@kk){print $O $k."\t".$CATabundance{$normMethod}{$y}{$k}."\n"; $CATcnt ++;  }#else {print "X";} }
				close $O;
				#print "Total of $CATcnt categories: failed $CATfail of ".($CATfail+$CATexist)." assignments\n";
			}
		}
	}
	#sleep(5);
	#is being counted double due to normethod counts

}

sub main(){
	my ($minBLE,$minPid,$pos)=@_;#$hr1,$hr2,$hr3) = @_;
	my %COGabundance=%{$COGhits{$pos}}; my %CATabundance=%{$CAThits{$pos}}; my %funHit=%{$GENEhits{$pos}};

#	my %COGabundance=%{$hr1}; my %CATabundance=%{$hr2}; my %funHit=%{$hr3};
	if (!$emode && !$scomode){die "No scoring mode selected!!\n";}
	#print "Scanning $DBmode with an eval of $minBLE\n";
	my $NOGtreat = 0; $NOGtreat = 1 if ($DBmode eq "NOG");
	#my %kgd = %{$kgdHR};
	
	my $bestSbj="";my $bestID=0; my $bestAlLen=0;my $bestQuery = "";
	my $COGexists=0;
	my $bestScore = 0; my $CBMmode = 0; my $bestE=1000; 
	
	my @tmp = ("");
	@tmp = @{$blRes[0]} if (@blRes > 0);
	my $qold=$tmp[0];
	my $COGfail=0; my $CATfail=0; my $totalCOG=0; my $CATexist=0;  my $ii=0;
	#this routine simply counts up number of hits to COGXX
	
	#while (1){ #don't need this, should only go once over this..
	if (1){
		#if ( $ii+1 >= @blRes){last;}#print "OUT   ";
		my $fndCat=0;
		
		for( ;$ii<@blRes;$ii++){
			#print "@{$blRes[$ii]}[0]\n";
			my ($Query,$Subject,$id,$AlLen,$mistmatches,$gapOpe,$qstart,$qend,$sstart,$send,$eval,$bitSc) = @{$blRes[$ii]};
			#print $eval."\n";
			#sort by eval #changed from bestE -> bestScore
			if ( ( ($emode && $bestE > $eval) || ($scomode && $bestScore< ($id * $AlLen)) ) 
						&& ($eval <= $minBLE && $bitSc >= $minScore)
						&& ($AlLen >= $minAlLen)
						&& ($id >= $minPid)
						&& ($quCovFrac == 0 || $AlLen > $DBlen{$Subject}*$quCovFrac) 
						&& (!$fndCat || $noHardCatCheck || exists $c2CAT{$Subject}) 
						) {
						#print "Y";
				$fndCat =1;
				$bestScore = $id * $AlLen;$bestSbj =$Subject; 
				$bestAlLen=$AlLen;$bestE = $eval; $bestQuery = $Query;
			}
			my @tmp;
			#if($ii+1 < @blRes){
			#	@tmp = @{$blRes[$ii+1]};
			#	if ($tmp[0] ne $qold){$qold=$tmp[0];last;}
			#} else {$qold="";last;}
		}
		#print @blRes ."  $bestSbj $bestQuery\n";
		#print "$ii ". ($ii+1 >= @blRes)." XX\n";
		if ($bestSbj ne ""){#finalize assignment to read
			my $curCOG="-";my $curCat = ""; my $curDef = "";
			my $curKgd = 3; #default always 3.. historic
			my @splC;
			if ($ACLdb || $KGMmode){
				@splC = split /:/,$bestSbj ; $splC[0] = $splC[1] if ($ACLdb);
			} elsif ($cazyDB || $tabCats==0) {
				@splC = split /\|/,$bestSbj;
			} 

			if ($cazyDB){
				if ($bestSbj =~ m/bacteria/){$curKgd =0;
				}elsif ($bestSbj =~ m/archaea/) {$curKgd =1;
				}elsif ($bestSbj =~ m/eukaryota/){ $curKgd =2 ;
					if (exists($czyTax{$splC[$#splC]})){$curKgd = $czyTax{$splC[$#splC]} ;}
				}else{ #take Mohs long list
					$curKgd = $czyTax{$splC[$#splC]} if (exists($czyTax{$splC[$#splC]}));
				}
				#else {print $bestSbj." ";}
				if ($bestSbj =~ m/\|CBM\d+/) {$CBMmode = 1; $curKgd = 4;}
			} elsif ( 0 && $ACLdb){ #cats are too primitive, doesn't need to be split up any further
				#67936 >protein:plasmid  25941 >protein:proph  28277 >protein:vir
				if ($bestSbj =~ m/plasmid/){$curKgd =0;
				}elsif ($bestSbj =~ m/proph/){$curKgd =1;
				}elsif ($bestSbj =~ m/vir/){$curKgd =2;
				}
			} elsif ($KGMmode){
				if (exists($KEGGtax{$splC[0]})){
					$curKgd = $KEGGtax{$splC[0]};
				}  #else {print "X";}
			}
			#print "XX$bestSbj"."XX\n";
			foreach my $normMethod (@normMethods){
				$totalCOG++; 
				my $score = 1;
				die "can't find length entry for $bestSbj\n" unless (exists $DBlen{$bestSbj});
				#print $DBlen{$bestSbj}." $normMethod ";
				$score = $bestAlLen / $DBlen{$bestSbj} if ($normMethod eq "GLN"); #score for this hit, normed by prot length
				if ($tabCats == 1){ #NOG
					$bestSbj =~ m/^(\d+)\./; my $taxid = $1;
					if (!exists $NOGkingd{$taxid}){print "Can't find $taxid in ref tax\n";}
					$curKgd = $NOGkingd{$taxid} if (!$CBMmode);
					#die "$bestSbj  $taxid $NOGkingd{$taxid}\n"; 
					unless (exists( $g2COG{$bestSbj} )){
						#print "can't find $bestSbj in NOG ref\n" unless ($COGfail>10);
						#print $O "$Query\t$Subject\t$bestID\t$bestE\t$bestAlLen\t\t\t\n";
						$COGfail++;
						#print "Too many COG(t/f)ails..\n" if  ($COGfail==11);
						next;
					}
					$curCOG = $g2COG{$bestSbj};
					$curDef = $COGdef{$curCOG} if (exists $COGdef{$curCOG});
					#hash{$_} = $valA for qw(a b c);
					#if (!exists($COGabundance{$normMethod}{$curCOG})){$COGabundance{$normMethod}{$curCOG}{$_}=0 foreach (@kgdOpts);}
					$COGabundance{$normMethod}{$curKgd}{$curCOG}+= $score;
					if (exists($c2CAT{$curCOG})){
						$curCat = $c2CAT{$curCOG};$CATexist++;
						#if (!exists($CATabundance{$normMethod}{$curCat})){$CATabundance{$normMethod}{$curCat}{$_}=0 foreach (@kgdOpts);}
						$CATabundance{$normMethod}{$curKgd}{$curCat} += $score;
					} else {
						$CATfail++;
						#print "Cat unknw for: $curCOG\n" ;
					}
					
					print $O2 "$bestQuery\t$curCOG\n" if ($reportGeneCat);
					#print $O "$Query\t$Subject\t$bestID\t$bestE\t$bestAlLen\t$curCOG\t$curCat\t$curDef\n";
				
				} elsif ($tabCats == 2){ #KEGG
					if (exists $c2CAT{$bestSbj} ){
						$curCat = $c2CAT{$bestSbj}; $CATexist++;
						#if (!exists($CATabundance{$normMethod}{$curCat})){  $CATabundance{$normMethod}{$curCat}{$_}=0 foreach (@kgdOpts);  }
						$CATabundance{$normMethod}{$curKgd}{$curCat} += $score;
					} else {$CATfail++; }#print " can't find kegg cat $bestSbj\n" unless ($CATfail>2);}
					print $O2 "$bestQuery\t$curCat\n" if ($reportGeneCat);
				} else { #Moh / CAZY / ACL
					$curCOG = $splC[0];
					#die $curCOG."\n";
					#$curCat = $curCOG;
					#if (!exists($COGabundance{$normMethod}{$curCOG})){$COGabundance{$normMethod}{$curCOG}{$_}=0 foreach (@kgdOpts);}
					$COGabundance{$normMethod}{$curKgd}{$curCOG}+= $score;
					#print $COGabundance{$normMethod}{$curKgd}{$curCOG}."\n";
					print $O2 "$bestQuery\t$curCOG\n" if ($reportGeneCat);
				}
				#if (!exists($funHit{$normMethod}{$bestSbj})){$funHit{$normMethod}{$bestSbj}{$_}=0 foreach (@kgdOpts); } 
				if ($tabCats == 0){ #don't need this for KEGG and NOG
					$funHit{$normMethod}{$curKgd}{$bestSbj}+= $score;
				}
			}
			$bestSbj="";$bestE=1000; $bestScore=0;$CBMmode=0;
			#$COGexists=0; 
		}
		$bestSbj="";$bestE=10; $bestScore=0;$CBMmode=0;
	}
	$totalCOG /= 2;$COGfail /= 2; $CATfail /= 2;
	#if (0&& $DBmode eq "NOG"){
	#	print "Total entries:$totalCOG\nCOG not found: $COGfail, CAT not found: $CATfail\n";
	#} 
	#die;
	#return (\%COGabundance, \%CATabundance, \%funHit);
	$COGhits{$pos} = \%COGabundance; $CAThits{$pos} = \%CATabundance;  $GENEhits{$pos} = \%funHit;

}


sub readCOGdef(){
	my ($inF) = @_;
	open I,"<$inF";
	my %ret;
	while (my $line = <I>){
		chomp $line;
		my ($txx,$COG,$n1,$n2,$cat,$def) = split(/\t/,$line);
		$ret{$COG} = $def;
	}
	close I;
	return \%ret;
}
sub readKeggTax{
	my ($inF) = @_;
	my %ret; my $eukCnt=0; my $proCnt=0;
	open I,"<$inF";
	while (<I>){
		chomp; my @spl = split /\t/;
		my $taxI = 0; 
		if ($spl[1] eq "euka"){$taxI = 1 ; $eukCnt++} else {$proCnt++;}
		$ret{$spl[0]} = $taxI;
	}
	close I;
	print "Euk in KEGG DB = $eukCnt; Pro = $proCnt\n";
	return \%ret;
}
sub readGene2KO(){
	my ($inF) = @_;
	my %c2cat;
	open I,"<$inF";
	print "reading gene -> KO file.. " or die "gene2ko file not present\n";
	while (my $line=<I>){
		chomp $line;
		my ($gn,$KO) = split(/\t/,$line);
		$KO =~ s/^ko://;
		$c2cat{$gn} = $KO;
	}
	print "Done\n";
	close I;
	return (\%c2cat);
}
sub readGene2COG(){
	my ($inF) = @_;
	my %g2c; my %c2cat;
	open I,"<$inF";
	print "reading gene -> COG file.. ";
	while (my $line=<I>){
		chomp $line;
		#NOG     COG5157 253     225     K       9796.ENSECAP00000013089,
		my ($cls,$COG,$p1,$p2,$cat,$genes) = split(/\t/,$line);
		$c2cat{$COG} = $cat;
		my @all = split(/,/,$genes);
		foreach (@all){
			#die "$_ $COG $cat\n";
			$g2c{$_} = $COG;
		}
	}
	print "Done\n";
	close I;
	return (\%g2c,\%c2cat);
}

sub readMohTax(){
	my ($inf) = @_;
	my %ret;
	open I,"<$inf" or die "Could not find Moh tax file at $inf\n";
	while (my $line = <I>){
		chomp $line;
		my @spl = split /\t/,$line;
		
		my $taxI =$czyEuk;#euk code
		$taxI = $czyFungi if ($spl[1] =~ m/fungi/i);
		$ret{$spl[0]} = $taxI;
	}
	close I;
	return \%ret;
}

sub readNogKingdom($){
	my ($inT) = @_;
	my %ret;
	open I,"<$inT" or die "Can't open NOG tax file\n";;
	while (my $line = <I>){
		my $cls = 0; #bacteria
		$cls = 1 if ($line =~ m/Archaea/);
		if ($line =~ m/Eukaryota/){
			$cls = 2;
			$cls =3 if ($line =~ m/Fungi/);
		}
		$line =~ m/^(\d+)\s/;
		$ret {$1} = $cls;
	}
	close I;
	return %ret;
}

sub combineBlasts($ $){
	my ($wh1,$wh2) = @_;
	my %bl1 = %{$wh1}; my %bl2 = %{$wh2};
	my %ret;
	
	my @allKs = uniq ( keys %bl1, keys %bl2); #
	#die "@allKs\n";
	foreach my $k (@allKs){
		my $ex1 = exists ($bl1{$k});
		unless ($ex1 && exists ($bl2{$k}) ){
			if ($ex1){$ret{$k} = $bl1{$k};
			} else { $ret{$k} = $bl2{$k};}
			next;
		}
		#pair
		$k =~ s/\/\d$/\/12/;
		$ret{$k} = $bl1{$k};
		my @hit1 = @{$bl1{$k}};
		my @hit2 = @{$bl2{$k}};
		
		my @sbss1 = sort{ $a <=> $b }(($hit1[8],$hit1[9]));
		my @sbss2 = sort{ $a <=> $b }($hit2[8],$hit2[9]);
		#sort 
		#print "$sbss1[0] > $sbss2[0]\n";
		if ($sbss1[0] > $sbss2[0]){
		#print"X";
			my @tmp = @hit1; @hit1 = @hit2; @hit2 = @tmp;
			@tmp = @sbss1; @sbss1 = @sbss2; @sbss2 = @tmp;
		}
		my @quss1 = sort { $a <=> $b }($hit1[6],$hit1[7]);my @quss2 = sort { $a <=> $b }($hit2[6],$hit2[7]);
		#overlap?
		my $overlap = 0;
		if ($sbss1[1] > $sbss2[0]){ $overlap= ( $sbss1[1] - $sbss2[0]); }#die "${$ret{$k}}[3] = $hit1[3] + $hit2[3] - ( $sbss1[1] - $sbss2[0])\n";}
		${$ret{$k}}[3] = $hit1[3] + $hit2[3] - $overlap; #ALlength
		#print "$sbss1[1] > $sbss2[0] $sbss1[0]  ${$ret{$k}}[3] $overlap\n";
		${$ret{$k}}[2] = ($hit1[2] + $hit2[2] ) /2;#%id
		${$ret{$k}}[11] = ($hit1[11] + $hit2[11]) * (1- $overlap/($hit1[3] + $hit2[3] ) );#bitscore
		${$ret{$k}}[10] = ($hit1[10], $hit2[10])[$hit1[10] > $hit2[10]];  # min($hit1[10] + $hit2[10]);
		#die "@{$ret{$k}}\n";
		#print "$k \n@sbss1 @sbss2\n@hit1\n";
	}
	
	return \%ret;
}

sub bestBlHit($){
	my ($hr) = @_;
	#my $emode=0; my $scomode = 1;
	my %blasts = %{$hr}; my $bestBit=0; my $bestk=""; my $bestScore=0;my $bestID=0; my $bestLen=0; my $bestIDever=0;
	my $k = "";
	foreach$k(keys %blasts){
		#my ($Query,$Subject,$id,$AlLen,$mistmatches,$gapOpe,$qstart,$qend,$sstart,$send,$eval,$bitSc) = ${$blRes{$k}}[11];
		#print $eval."\n";
		#sort by eval #changed from bestE -> bestScore
		#if ( ( ($emode && $bestE > $eval) || ($scomode && $bestScore< ($id * $AlLen)) ) 
		#			#&& ($eval <= $minBLE || $bitSc >= $minScore)
		#			&& ($quCovFrac == 0 || $AlLen > $DBlen{$Subject}*$quCovFrac) 
		#			&& (($fndCat || $noHardCatCheck) || exists $c2CAT{$Subject}) ) {
		#			#print "Y";
		
		#my $curScore = ${$blasts{$k}}[2] * ${$blasts{$k}}[3];
		#if ($bestScore < $curScore){$curScore = $bestScore;$bestk = $k;	}
		
		#	$bestAlLen=$AlLen;$bestE = $eval;#$bestQuery = $Query;
		
		#just sort by bitscore
		#print "@{$blasts{$k}}\n";
		#if ($bestBit < ${$blasts{$k}}[11]){
		#	$bestBit = ${$blasts{$k}}[11];  $bestk = $k;
		#}
		my $cID = ${$blasts{$k}}[2]; my $cLe = ${$blasts{$k}}[3]; my $cSc=${$blasts{$k}}[11];
		if (($bestID -5)< $cID){
			if ( ( $cID >= ($bestID *0.97) && $cID >= $bestIDever*0.9   && ( $cLe >= $bestLen * 1.15 ) )  ||  #length is just better (15%+)
					(   ($cID >= $bestID *1.03) && ( $cLe >= $bestLen * 0.9) ) ||#id is just better, while length is not too much off
					$cSc > $bestBit*0.8){   #convincing score
				$bestID = $cID;  $bestk = $k; $bestLen = $cLe; $bestBit = $cSc;
				if ($bestID > $bestIDever){$bestIDever=$bestID;}
			}
		}
	}
	if ($bestk eq ""){die "something went wrong with ABR blast:\n$k:$blasts{$k}\n";}
	return $blasts{$bestk};
}
