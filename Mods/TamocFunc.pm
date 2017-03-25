package Mods::TamocFunc;
use warnings;
#use Cwd 'abs_path';
use strict;
#use List::MoreUtils 'first_index'; 

#use Mods::GenoMetaAss qw(qsubSystem);

use Exporter qw(import);
our @EXPORT_OK = qw(sortgzblast uniq getE100);
use Mods::GenoMetaAss qw(systemW );
use Mods::IO_Tamoc_progs qw(getProgPaths);

sub getE100($ $ $){
	my ($oDess,$proteins,$genesNT) = @_;
	my $hmmbin = getProgPaths("hmmsearch");#"/g/bork5/hildebra/bin/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmsearch";
	my $essDB = getProgPaths("essentialHMM");#"/g/bork5/hildebra/bin/multi-metagenome-master/R.data.generation/essential.hmm";
	my $essEukDB = getProgPaths("essentialEUK");#"/g/bork3/home/hildebra/DB/HMMs/eukCore/eukCore.hmm"; #TODO
	systemW("mkdir -p $oDess");
	my $cmd = "";
	$cmd .= "$hmmbin --domtblout $oDess/assembly.hmm.orfs.txt -o $oDess/assembly.hmmsout.txt --cut_tc --cpu 1 --notextw $essDB $proteins\n";
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

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
sub sortgzblast{ #function that checks if the diamond output was already sorted (required for paired end stuff with reads)
	my ($input) = @_;
	#print "$input\n";
	if ( $input =~ m/\.srt\.gz$/ ) {
		my $trial = $input; $trial =~ s/\.srt//; my $trialuse=0;
		if (-e $trial && -e $input && (-s $trial > -s $input)){$input = $trial; $trialuse=1; }#print "trial\n";
		if (-e $input && !$trialuse){return $input; }
		if (-e $trial && !-e $input){$input = $trial;}
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










