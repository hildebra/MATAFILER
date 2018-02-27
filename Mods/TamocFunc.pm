package Mods::TamocFunc;
use warnings;
#use Cwd 'abs_path';
use strict;
#use List::MoreUtils 'first_index'; 

#use Mods::GenoMetaAss qw(qsubSystem);

use Exporter qw(import);
our @EXPORT_OK = qw( sortgzblast attachProteins attachProteins2 uniq getFMG readTabbed readTabbed2 readTabbed3 getSpecificDBpaths renameFMGs);
use Mods::GenoMetaAss qw(systemW readFasta renameFastHD gzipwrite gzipopen);
use Mods::IO_Tamoc_progs qw(getProgPaths);


#gets required FA genes via samtools faidx and saves them to file
#Usage: [text file with target genes] [out faa] [in faa(search genes in this fna)] [new names for output]


sub attachProteins2{
	my ($inT,$prF,$protIn) = ($_[0],$_[1],$_[2]);
	
	my %gene2num; my $doRename=0;
	if (@_ > 3){
		$doRename=1;
		my $hr = $_[3];
		%gene2num = %{$hr};
	}
	die "Protein file does not exist: $protIn\n" unless (-e $protIn);
	my $hr = readFasta($protIn);
	my %fas = %{$hr};
	#print "$protIn\n";
	#my $protStore = `cat $inT | xargs $samBin faidx $protIn`;
	#die length($protStore)."\n";
	my @prots = @{$inT};
	open O,">>$prF" or die "Can't open $prF\n";
	foreach my $pr (@prots){
		$pr =~ s/^>//;
		#1 identify protein
		my $seq = $fas{$pr};
		if ($doRename){
			unless(exists($gene2num{$pr})){die "can not identify $pr gene in index file while rewritign prot names\n$protIn\n";}
			#print "$gene2num{$spl[0]}\n $pr\n";
			$pr = $gene2num{$pr};
		} 
		print O ">".$pr."\n".$seq."\n";
	}
	close O;
}

sub attachProteins{#"$basD/tmp.txt",$protF){
	my ($inT,$prF,$protIn) = ($_[0],$_[1],$_[2]);
	
	my %gene2num; my $doRename=0;
	if (@_ > 3){
		$doRename=1;
		my $hr = $_[3];
		%gene2num = %{$hr};
	}
	my $samBin = getProgPaths("samtools");
	#the old way
	#systemW("cat $basD/tmp.txt | xargs samtools faidx $protIn  >> $prF");
	die "Protein file does not exist: $protIn\n" unless (-e $protIn);
	#print "$protIn\n";
	my $protStore = `cat $inT | xargs $samBin faidx $protIn`;
	#die length($protStore)."\n";
	my @prots = split />/,$protStore;
	open O,">>$prF" or die "Can't open $prF\n";
	foreach my $pr (@prots){
		#1 identify protein
		next if (length($pr) <= 1);
		my @spl = split /\n/,$pr;
		if ($doRename){
			unless(exists($gene2num{$spl[0]})){die "can not identify $spl[0] gene in index file while rewritign prot names\n$protIn\n";}
			#print "$gene2num{$spl[0]}\n $pr\n";
			$spl[0] = ">".$gene2num{$spl[0]};
		} else {
			$spl[0] = ">$spl[0]";
		}
		print O join("\n",@spl)."\n";
	}
	close O;
	#die "cat $inT | xargs $samBin faidx $protIn\n";
}

sub readTabbed($){
	my ($inF) = @_;
	my %ret;
	open It,"<$inF" or die "Cant open tabbed infile $inF\n";
	while (<It>){
		chomp;
		my @spl = split /\t/;
		$ret{$spl[0]} = $spl[1];
	}
	close It;
	return \%ret;
}
sub readTabbed3($ $){
	my ($inF,$col) = @_;
	my %ret;
	my ($It,$OK) = gzipopen($inF,"Tabbed file",1);
	#open It,"<$inF" or die "Cant open tabbed infile $inF\n";
	while (<$It>){
		chomp;
		my @spl = split /\t/;
		if (@spl < $col){die"Requested column $col, but file $inF has only ".@spl." columns\n";}
		$ret{$spl[0]} = $spl[$col];
	}
	close $It;
	return \%ret;
}
sub readTabbed2{
	my $inF = $_[0];
	my $useLast = 0;
	$useLast = $_[1] if (@_>1);#not used currently
	
	my %ret;my $maxDepth=0;
	open It,"<$inF" or die "Cant open tabbed infile $inF\n";
	
	while (<It>){
		chomp;
		my @spl = split /\t/;
		my $ke;
		if ($useLast==1){
			$ke	= pop @spl;
		} else {
			$ke	= shift @spl;
		}
		#print "$ke\n";
		#print $ke." ";
		$ret{$ke} = \@spl;
		#die "@spl\n";
		if (scalar(@spl) > $maxDepth){$maxDepth=scalar(@spl);}
	}
	close It;
	return (\%ret,$maxDepth);
}


sub getSpecificDBpaths($ $){
	my ($curDB,$checkDBpreped) = @_;
	if ($curDB eq "mp3"){return "","",$curDB;}
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
	elsif ($curDB eq "ABRc"){$DBpath = getProgPaths("ABRcard_path_DB"); $refDB = "card.parsed.f11.faa";$shrtDB = $curDB; }
	elsif ($curDB eq "KGE"){$DBpath = getProgPaths("KEGG_path_DB"); $refDB = "genus_eukaryotes.pep";$shrtDB = $curDB; }
	elsif ($curDB eq "KGB"){$DBpath = getProgPaths("KEGG_path_DB"); $refDB = "species_prokaryotes.pep";$shrtDB = $curDB; }
	elsif ($curDB eq "ACL"){$DBpath = getProgPaths("ACL_path_DB");$refDB = "aclame_proteins_all_0.4.fasta";$shrtDB = $curDB; }
	elsif ($curDB eq "KGM"){$DBpath = getProgPaths("KEGG_path_DB"); $refDB = "euk_pro.pep";$shrtDB = $curDB; }
	elsif ($curDB eq "TCDB"){$DBpath = getProgPaths("TCDB_path_DB"); $refDB = "tcdb.faa";$shrtDB = $curDB; }
	elsif ($curDB eq "PTV"){$DBpath = getProgPaths("PATRIC_VIR_path_DB"); $refDB = "PATRIC_VF.faa";$shrtDB = $curDB; }
	elsif ($curDB eq "PAB"){$DBpath = getProgPaths("ABprod_path_DB"); $refDB = "dedup_best_prod_predictions.faa";$shrtDB = $curDB; }
	else {die"Unknown DB for Diamond: $curDB\n";}
	
	if ($checkDBpreped){
		die "Can't find prepared diamond database at:\n$DBpath$refDB.db.dmnd" unless (-e "$DBpath$refDB.db.dmnd");
		die "Can't find length file at:\n$DBpath$refDB.length" unless (-e "$DBpath$refDB.length");
	}
	return ($DBpath ,$refDB ,$shrtDB );
}


sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
sub sortgzblast{ #function that checks if the diamond output was already sorted (required for paired end stuff with reads)
	my ($input) = @_;
	#print "$input\n";
	if ( $input =~ m/\.srt\.gz$/ ) { #redo srt in case there's a $trial file
		my $trial = $input; $trial =~ s/\.srt//; my $trialuse=0;
		if (!-e $input && !-e $trial){die "$input doesn't exist!\n";}
		if (-e $trial && -e $input && (-s $trial > -s $input)){$input = $trial; $trialuse=1; }#print "trial\n";
		if (-e $input && !$trialuse){return $input; }
		if (-e $trial && !-e $input){$input = $trial;}
		die "something went wrong in gzip sort\n$input\n" unless (-e $input);
	}
	my $tmpd=""; 
	my $cmd = "";
	my @chars = ("A".."Z", "a".."z");my $randstring;
	$randstring .= $chars[rand @chars] for 1..8;
	my $tmpDset=0;
	if (@_ >= 2){
		$tmpd = $_[1];$tmpDset=1;
	} 
	if ($tmpd eq ""){
		$input =~ m/^(.*\/)[^\/]+$/;$tmpd = $1;
	}
	my $input2=$input;
	$input2 =~ s/\.gz$//;
	my $input3=$input;
	$input3 =~ s/\.srt\.gz$//;
	if (!-e $input){ #maybe already something done here..
		if (!-e "$input2.srt.gz" && -e "$input.srt.gz"){system "mv $input.srt.gz $input2.srt.gz";}
		if (-e "$input2.srt.gz"){$input = "$input2.srt.gz";
		} elsif ( -e $input3 ){$input = $input3;
		}
	}
	#print $input."\n";
	unless ($input =~ m/\.srt\.gz$/){ #do sort (and maybe gz)
		if ($input =~ m/\.srt$/){
			$cmd = "gzip $input"; $input .= ".gz";
		} elsif ($input =~ m/\.gz$/) { #not sorted, but gz
			system "mkdir -p $tmpd" unless (-d "$tmpd");
			my $tmpf = "$tmpd/rawBLast$randstring.bla";
			$cmd = "zcat $input > $tmpf; sort $tmpf | gzip > $input2.srt.gz; rm -f $input $tmpf; ";
			if (!-e $input){die "Wrong file as input provided: $input\n";}
		} else { #not gz, not sort
			$cmd = "sort $input > $input.srt; gzip $input.srt; rm $input;";
		}
	}
	#die $cmd."\n$input\n";
	unless ($cmd eq ""){
		if (system $cmd) { die "$cmd \nfailed\n"; }
	}
	$input = "$input2.srt.gz";
	die "Something went wrong in sortgzblast\n" if (!-e $input);
	return $input;
}










