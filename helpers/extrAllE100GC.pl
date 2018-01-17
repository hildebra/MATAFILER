#!/usr/bin/env perl
#extracts abundance of all essential marker genes, as well as FMGs
#./extrAllE100GC.pl /g/scb/bork/hildebra/SNP/GCs/GNM3_ABR
#./extrAllE100GC.pl /g/scb/bork/hildebra/SNP/GCs/
use warnings;
use strict;
sub processSubGenes;
sub getGeneSeqsSubGenes;

use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(readMapS systemW readClstrRev);
my $smtBin = getProgPaths("samtools");#/g/bork5/hildebra/bin/samtools-1.2/samtools";
my $rarBin = getProgPaths("rare");#"/g/bork5/hildebra/dev/C++/rare/rare";

my $GCd = $ARGV[0];
my $oldNameFolders = -1;#$ARGV[1];
my $mapF = `cat $GCd/LOGandSUB/GCmaps.inf`;
my ($hrm,$hr2X) = readMapS($mapF,$oldNameFolders);
my %map = %{$hrm};
my @samples = @{$map{smpl_order}};

#my $GCd = "$inD/GeneCatalog/";
#getGeneSeqsSubGenes("FMG");die();

my ($hr1,$hr2) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx");
my %gene2cl = %{$hr1};
my %genesE1h; my %genesFMG;
#$genesE1h{1}{gg} = "falk";
#die $genesE1h{1}{gg};
foreach my $smpl(@samples){
	my $inD = $map{$smpl}{wrdir};
	my $SmplName = $map{$smpl}{SmplID};
	my $metaGD = `cat $inD/assemblies/metag/assembly.txt`; chomp $metaGD;

	#print $SmplName."\n";
	my $e1f = "$metaGD/ContigStats/ess100genes/ess100.id.txt";
	my $FMGf = "$metaGD/ContigStats/FMG/FMGids.txt";
#read in ess 100 genes in assembly
	if (-f $e1f){
		open my $I,"<$e1f" or die "Can't open e100 file $e1f\n";
		while (my $l = <$I>){
			chomp $l;next if (length($l) < 5);
			my @spl = split /\s/,$l;
			$genesE1h{$spl[1]}{$spl[0]}=1; 
		}	close $I;
	} else {print "Can't find e100 file $e1f\n";}
#read in FMG genes in assembly
	if (-f $FMGf){
		open my $I,"<$FMGf" or die "Can't open FMG file $FMGf\n";
		while (my $l = <$I>){
			chomp $l;next if (length($l) < 5);
			my @spl = split /\s/,$l;
			$genesFMG{$spl[1]}{$spl[0]}=1; 
		}	close $I;
	} else {print "Can't find FMG file $FMGf\n";}
}
print "Read ref dataset\n";

#print "@e1cat\n";
print "Creating FMG gene matrix\n";
processSubGenes(\%genesFMG,"FMG");
print "extracting FNA & FAA's of FMG genes\n";
getGeneSeqsSubGenes("FMG");


die "Done FMG\n";
print "Creating e100 gene matrix\n";
processSubGenes(\%genesE1h,"e100");
print "Finished\n";
exit(0);


sub getGeneSeqsSubGenes(){
	my ($tag) = @_;
	my $subF = "$GCd/$tag.subset.cats";
	my $fmgOD = "$GCd/$tag/";
	#die "TODO getGeneSeqsSubGenes\n";
	system "mkdir -p $fmgOD"; 
	
	open I,"<$subF" or die "can't open $subF\n"; 
	while (my $line=<I>){
		chomp $line;
		my @spl = split /\t/,$line;
		my @spl2 = split /,/,$spl[2];
		#die "\n@spl2\n";
		my $ofile = $fmgOD."/$spl[0]";
		system "$smtBin faidx $GCd/compl.incompl.95.fna ". join (" ", @spl2) . " > $ofile.fna";
		system "$smtBin faidx $GCd/compl.incompl.95.prot.faa ". join (" ", @spl2) . " > $ofile.faa";
	} 
	close I;
}

sub processSubGenes(){
	my ($ghr,$tag) = @_;
	my %genes = %{$ghr};
	my $subF = "$GCd/$tag.subset.cats";
	my @e1cat = keys %genes;

	my %selC; open O ,">$subF" or die "Can't open $subF";
	foreach my $e1c (@e1cat){
		my @spG = keys %{$genes{$e1c}};
		my %selD;   #print "\n";
		foreach my $gene (@spG){
			die $gene. " ".$gene2cl{">".$gene}."\n" unless exists ($gene2cl{">".$gene});
			#print $gene."\n";
			$selD{$gene2cl{">".$gene}} = 1;
		}
		my @keysHds = keys %selD;
		my $size = @keysHds;
		if ($size == 0){die "@spG"."\n";}
		my $addD = $keysHds[0];
		if ($size > 1 ){$addD = join(",",sort(@keysHds));}
		#print $addD  ."    @keysHds\n";
		print O $e1c."\t$size\t".$addD."\n";
		my %tmp = (%selC, %selD); %selC = %tmp;
	#die ( @spG."\n");
	}
	close O;
	#my $sedStr = join("p;",@rows);
	#system O "sed -n '1p;$sedStr"."p' $GCd/Matrix.mat > $GCd/e100subset.mat";
	my @rows = keys %selC;
	push(@rows,1);
	print "Selected ".@rows." rows.. writing to $GCd/$tag.subset.mat\n";
	open O,">$GCd/$tag.lines" or die "Can't open output lines file\n";;print O join("\n",@rows);close O;
	
	my $cmdX= "$rarBin lineExtr -i $GCd/Matrix.mat -o $GCd/$tag.subset.mat -reference $GCd/$tag.lines\n rm $GCd/$tag.lines";
	print $cmdX."\n";
	system $cmdX;
	print "Done $tag\n";
}
