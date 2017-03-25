#!/usr/bin/perl -w
use strict;
use LWP::Simple;
use FindBin qw($Bin);

require("$Bin\/_getmarker.pl");
require("$Bin\/_getabund.pl");
require("$Bin\/_sepReads.pl");

my $MAX_MARKER = 0;
my $LOGOUT;

my $SETTING_FILE = "setting";
my $BOWTIE2BUILD = "bowtie2-build";
my $BOWTIE2 = "bowtie2";
my $HMMSEARCH = "hmmsearch";
my $RUNFRAG = "run_FragGeneScan.pl";
my $VELVETH = "velveth";
my $VELVETG = "velvetg";

#don't need this
#checkProgram();

my $RSCRIPT = "Rscript";
my $HEATMAP_R = "$Bin\/heatmap.r";
my $MARKERHMM = "$Bin/marker.hmm";
my $MARKERNUM = 107;
my $MAXBIN = "$Bin\/src\/MaxBin";
my $ABUND_OUTPUT = "";
my $THREADNUM = 1;
my $REASSEMBLY = 0;
my $RPLOT = 0;
my $KMERLEN = 55;
my $MIN_SEQ_LENGTH = 1000;
my $MIN_BIN_SIZE = 100000;
my $REASSEM_MIN_LEN = 200;
my $geneAAfile = "";
my $essHMMfile = "";

my $USAGE = qq(MaxBin - a metagenomics binning software.
Usage:
  run_MaxBin.pl
    -contig (contig file)
    -out (output file)
    [-reads (reads file)]
    [-abund (abundance file)]
    [-min_contig_length (minimum contig length. Default 1000)]
    [-max_iteration (maximum Expectation-Maximization algorithm iteration number. Default 50)]
    [-thread (thread num; default 1)]
    [-plotmarker]
    [-reassembly]
    [-reassembly_kmer (kmer for reassembly; default 55 -- For larger kmers please reconfigure velvet package in [auxiliary] or other user directory)]
    [-markerset (marker gene sets, 107 (default) or 40.  See README for more information.)]

  for debug purpose:
    [-verbose]
    [-preserve_intermediate]

  Please specify either -reads or -abund information.
  Please read README file for more details.\n);


main();

sub main
{
	my $contig_f = "";
	my $contig_name = "";
	my $reads_f = "";
	my $abund_f = "";
	my $out_f = "";
	my $outname = "";
	my $verbose = 0;
	my $preserve = 0;
	my $maxem = -1;
	my $remodel_contig = 0;
	my $new_contig;
	my $old_contig_name;

	my $starttime = time();
	my $endtime;
	
	my $ARGC = @ARGV;
	my $i;
	my $j;
	my $k;
	my $cmd;
	my $line;
	my $param;

	# Create a temporary log file
	$k = 0;
	while ($k == 0)
	{
		$k = int(rand(100000000));
		if (-e "$k.log")
		{
			$k = 0;
		}
	}
	open(TMPLOG, ">$k.log");

	# Check if MaxBin program exist. If not then build everything, including 3rd-party software
	if (!(-e $MAXBIN))
	{
		print "Program not built. Please run \"make\" under src directory to build MaxBin program.\n";
		close(TMPLOG);
		unlink("$k.log");
		exit(-1);
	}
	
	for ($i = 0; $i < $ARGC; $i++)
	{
		if ($ARGV[$i] eq "-contig")
		{
			$i++;
			$contig_f = $ARGV[$i];
			$contig_name = $ARGV[$i];
			print "Input contig: $contig_f\n";
			print TMPLOG "Input contig: $contig_f\n";
		}
		elsif ($ARGV[$i] eq "-reads")
		{
			$i++;
			if (-e $ARGV[$i])
			{
				$reads_f = $ARGV[$i];
				print "Located reads file [$ARGV[$i]]\n";
				print TMPLOG "Located reads file [$ARGV[$i]]\n";
			}
			else
			{
				print "Cannot find reads file [$ARGV[$i]]\n";
				close(TMPLOG);
				unlink("$k.log");
				exit;
			}
		}
		elsif ($ARGV[$i] eq "-abund")
		{
			$i++;
			if (-e $ARGV[$i])
			{
				$abund_f = $ARGV[$i];
			}
			else
			{
				print "Cannot find abundance file [$ARGV[$i]]\n";
				close(TMPLOG);
				unlink("$k.log");
				exit;
			}
		}
		elsif ($ARGV[$i] eq "-out")
		{
			$i++;
			$out_f = $ARGV[$i];
			$j = rindex($out_f, "/");
			$outname = substr($out_f, $j + 1);
			print "out header: $ARGV[$i]\n";
			print TMPLOG "out header: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-thread")
		{
			$i++;
			$THREADNUM = $ARGV[$i];
			print "Thread: $ARGV[$i]\n";
			print TMPLOG "Thread: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-HMM"){
			$i++; $essHMMfile = $ARGV[$i];
			print TMPLOG "HMMfile: $ARGV[$i]\n";
		}
		
		elsif ($ARGV[$i] eq "-AA"){
			$i++; $geneAAfile = $ARGV[$i];
			print TMPLOG "AAfile: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-max_iteration")
		{
			$i++;
			$maxem = $ARGV[$i];
			print "Max iteration: $ARGV[$i]\n";
			print TMPLOG "Max iteration: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-min_contig_length")
		{
			$i++;
			$MIN_SEQ_LENGTH = $ARGV[$i];
			print "Min contig length: $ARGV[$i]\n";
			print TMPLOG "Min contig length: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-reassembly") # Reassembly currently still does not support FASTQ reads
		{
			$REASSEMBLY = 1;
			print "Reassembly: 1\n";
			print TMPLOG "Reassembly: 1\n";
		}
		elsif ($ARGV[$i] eq "-reassembly_kmer") # Reassembly kmer
		{
			$i++;
			$KMERLEN = $ARGV[$i];
			print "Reassembly kmer: $KMERLEN\n";
			print TMPLOG "Reassembly kmer: $KMERLEN\n";
		}
		elsif ($ARGV[$i] eq "-markerset") 
		{
			$i++;
			if ($ARGV[$i] == 40)
			{
				print "Switch to 40 marker genes universal for bacteria and archaea.\n";
				print TMPLOG "Switch to 40 marker genes universal for bacteria and archaea.\n";
				$MARKERHMM = "$Bin/bacar_marker.hmm";
				$MARKERNUM = 40;
			}
		}
		elsif ($ARGV[$i] eq "-verbose")
		{
			$verbose = 1;
		}
		elsif ($ARGV[$i] eq "-plotmarker")
		{
			$RPLOT = 1;
		}
		elsif ($ARGV[$i] eq "-preserve_intermediate")
		{
			$preserve = 1;
		}
		else
		{
			print "Unrecognized token \[$ARGV[$i]\]\n";
			print $USAGE;
			close(TMPLOG);
			unlink("$k.log");
			exit;
		}
	}
	
	if ($contig_f eq "")
	{
		print "No Contig file. Please specify contig file by -contig\n";
		print $USAGE;
		close(TMPLOG);
		unlink("$k.log");
		exit;
	}
	if ($out_f eq "")
	{
		print "Please specify output file by -out.\n";
		print $USAGE;
		close(TMPLOG);
		unlink("$k.log");
		exit;
	}
	if (-e "$out_f.log")
	{
		unlink "$out_f.log";
	}
	
	if ($reads_f eq "" && $abund_f eq "")
	{
		print "Please input either reads file or abundance file.\n";
		print $USAGE;
		close(TMPLOG);
		unlink("$k.log");
		exit;
	}
	if ($REASSEMBLY == 1 && $reads_f eq "")
	{
		print "Please specify \"interleaved paired-end fasta file\" for reassembly purpose.\n";
		close(TMPLOG);
		unlink("$k.log");
		exit;
	}

	openLOG("$out_f.log");
	close(TMPLOG);
	open(TMPLOG, "<$k.log");
	while(defined($line = <TMPLOG>))
	{
		print $LOGOUT $line;
	}
	close(TMPLOG);
	unlink("$k.log");
	$new_contig = checkContig($contig_f, $MIN_SEQ_LENGTH, "$out_f.tooshort");
	if ($new_contig ne "")
	{
		$old_contig_name = $contig_f;
		$contig_f = $new_contig;
		$remodel_contig = 1;
	}
	
	$param = "";
	if ($abund_f eq "")
	{
		writeLOG("Running Bowtie2...this may take a while...\n");
		# Run bowtie2 to find abundance information
		runBowtie2($old_contig_name, $reads_f, "$out_f.sam");
		getsam("$out_f.sam", "$contig_f.abund");
	}
	else
	{
		checkAbundFile($abund_f, "$contig_f.abund");
	}
	$param = $param . " -abund " . $contig_f . ".abund";

	my @binarr;
	my $currbin;
	my $currout;
	my $currnum;
	my $maxnum;
	my @summarr;
	my @reclassifyarr;
	my @noclassarr;
	push(@binarr, $contig_f);
	# Push some dummy into result array for the original dataset
	push(@summarr, "dummy");
	push(@reclassifyarr, 1);
	$currnum = 0;
	$maxnum = 1;

	# Run HMMER3 to identify seed contigs
	writeLOG("Searching against $MARKERNUM marker genes to find starting seed contigs for [$contig_name]...\n");
	if (-e "$contig_f.frag.faa")
	{
		unlink("$contig_f.frag.faa");
	}
	getHMM($contig_f, "$contig_f.hmmout",$geneAAfile) unless (-e $essHMMfile);#$essHMMfile

	while ($currnum < $maxnum)
	{
		$currbin = $binarr[$currnum];
		if ($currbin eq $contig_f)
		{
			$currout = $out_f;
		}
		else
		{
			$currout = substr($currbin, 0, length($currbin) - 6) . ".out";
		}

		if ($MAX_MARKER == 1)
		{
			$i = gethmmmarker("$essHMMfile", $currbin, $MIN_SEQ_LENGTH, "$currout.seed", 1);
		}
		else
		{
			$i = gethmmmarker("$essHMMfile", $currbin, $MIN_SEQ_LENGTH, "$currout.seed");
		}
		if ($currbin eq $contig_f && $i == -1)
		{
			writeLOG("Try harder to dig out marker genes from contigs.\n");
			$i = gethmmmarker("$essHMMfile", $currbin, $MIN_SEQ_LENGTH, "$currout.seed", 1);
			if ($i == -1)
			{
				writeLOG("Marker gene search reveals that the dataset cannot be binned (the medium of marker gene number <= 1). Program stop.\n");
				closeLOG("$out_f.log");
				exit(-1);
			}
		}
		elsif ($currbin ne $contig_f && $i == -1)
		{
			$currnum++;
			next;
		}
		elsif ($i != -1)
		{
			$reclassifyarr[$currnum] = 1;
		}

		# Running MaxBin
		writeLOG("Done data collection. Running MaxBin...\n");
		if ($verbose == 0)
		{
			if ($maxem == -1)
			{
				$cmd = "$MAXBIN -fasta $currbin $param -seed $currout.seed -out $currout -min_contig_length $MIN_SEQ_LENGTH";
			}
			else
			{
				$cmd = "$MAXBIN -fasta $currbin $param -seed $currout.seed -out $currout -min_contig_length $MIN_SEQ_LENGTH -max_run $maxem";
			}
		}
		else
		{
			if ($maxem == -1)
			{
				$cmd = "$MAXBIN -fasta $currbin $param -seed $currout.seed -out $currout -min_contig_length $MIN_SEQ_LENGTH -verbose";
			}
			else
			{
				$cmd = "$MAXBIN -fasta $currbin $param -seed $currout.seed -out $currout -min_contig_length $MIN_SEQ_LENGTH -verbose -max_run $maxem";
			}
		}
		writeLOG("Command: $cmd\n");
		system($cmd);
		if (checkResult("$currout.summary") == -1)
		{
			if ($currbin eq $contig_f)
			{
				writeLOG("Error encountered while running core MaxBin program. Error recorded in $currout.log.\nProgram Stop.\n");
				exit(-1);
			}
		}
		push(@noclassarr, $currout);

		# Read in summary file and putting bins into stack
		open(FILE, "<$currout.summary");
		while(defined($line = <FILE>))
		{
			chomp($line);
			if ($line =~ /^Bin \[([A-Za-z0-9._\\\/\@\!\|\#\$\%\^\?\<\>\[\]\{\}\(\)\+\-]+)\] \(([0-9.]+)\)/)
			{
				push(@binarr, $1);
				push(@summarr, $2);
				push(@reclassifyarr, 0);
				$maxnum++;
			}
		}
		close(FILE);

		$currnum++;
	}

	# Put everything in order

	# Remove re-classified bins
	my @tmpgenomesize;
	my @tmpgc;
	my @genomesize;
	my @gc;
	my @bin_to_noclass;
	for ($k = 0; $k < $maxnum; $k++)
	{
		($tmpgenomesize[$k], $tmpgc[$k]) = getBinInfo($binarr[$k]);
	}
	for ($k = 0; $k < $maxnum; $k++)
	{
		if ($reclassifyarr[$k] == 1 || $tmpgenomesize[$k] < $MIN_BIN_SIZE)
		{
			if ($reclassifyarr[$k] == 1 && $binarr[$k] ne $contig_f)
			{
				unlink($binarr[$k]);
			}
			elsif ($reclassifyarr[$k] == 0)
			{
				push(@bin_to_noclass, $binarr[$k]);
			}
			splice(@reclassifyarr, $k, 1);
			splice(@binarr, $k, 1);
			splice(@summarr, $k, 1);
			splice(@tmpgenomesize, $k, 1);
			splice(@tmpgc, $k, 1);
			$k--;
			$maxnum--;
		}
	}

	# First sort @summarr and store the index in another array
	my @sortarr = sort {$b <=> $a} @summarr;
	my %indexhash;
	@indexhash{@sortarr} = (0..$#summarr);
	@sortarr = ();

	my @filearr;
	for ($k = 0; $k < $maxnum; $k++)
	{
		$i = $indexhash{$summarr[$k]};
		$j = $i + 1;
		# rename original fasta file
		rename ($binarr[$k], "$out_f.$j.fasta.tmp");
		$filearr[$i] = "$out_f.$j.fasta.tmp";
		$genomesize[$i] = $tmpgenomesize[$k];
		$gc[$i] = $tmpgc[$k];
		$sortarr[$i] = $summarr[$k];
	}
	# Write summary
	open(TMPSUM, ">$out_f.tmp.summary");
	for ($i = 0; $i < $maxnum; $i++)
	{
		$j = $i + 1;
		print TMPSUM "Bin \[$out_f.$j.fasta\] \($sortarr[$i]\)\n";
	}
	close(TMPSUM);

	# Collect all no-classes and log file
	open(OUT, ">$out_f.tmp.noclass");
	$j = @noclassarr;
	for ($i = 0; $i < $j; $i++)
	{
		open(FILE, "<$noclassarr[$i].noclass");
		while(defined($line = <FILE>))
		{
			print OUT $line;
		}
		close(FILE);

		open(FILE, "<$noclassarr[$i].log");
		while(defined($line = <FILE>))
		{
			writeLOG($line);
		}
		writeLOG("\n");
		close(FILE);
	}
	foreach $currbin (@bin_to_noclass)
	{
		if (-e $currbin)
		{
			open(FILE, "<$currbin");
			while(defined($line = <FILE>))
			{
				print OUT $line;
			}
			close(FILE);
			unlink($currbin);
		}
		else
		{
			writeLOG("File $currbin not found.\n");
		}
	}
=cut
	if ((-e "$out_f.tooshort") && (-s "$out_f.tooshort" > 0))
	{
		open(FILE, "<$out_f.tooshort");
		while(defined($line = <FILE>))
		{
			print OUT $line;
		}
		close(FILE);
	}
	unlink("$out_f.tooshort");
=cut
	close(OUT);

	# delete all normal files
	#for ($i = 0; $i < $maxnum; $i++)
	#{
	#	if ($binarr[$i] ne $contig_f && (-e $binarr[$i]))
	#	{
	#		unlink($binarr[$i]);
	#	}
	#}
	$j = @noclassarr;
	for ($i = 0; $i < $j; $i++)
	{
		unlink("$noclassarr[$i].log");
		unlink("$noclassarr[$i].summary");
		unlink("$noclassarr[$i].noclass");
		unlink("$noclassarr[$i].prob_dist");
		unlink("$noclassarr[$i].dist");
		unlink("$noclassarr[$i].seed");
	}

	# rename all tmp files to normal files
	rename("$out_f.tmp.summary", "$out_f.summary");
	rename("$out_f.tmp.noclass", "$out_f.noclass");
	for ($i = 0; $i < $maxnum; $i++)
	{
		$currbin = substr($filearr[$i], 0, length($filearr[$i]) - 4);
		rename($filearr[$i], $currbin);
		$filearr[$i] = $currbin;
	}
	
	# Separate reads for reassembly
	my %lenhash;
	if ($REASSEMBLY == 1)
	{
		if ($reads_f ne "")
		{
			writeLOG("Performing reassembly. Reads file found.\n");
			writeLOG("Separating reads according to the bins...\n");
			if ($abund_f ne "")
			{
				runBowtie2($contig_f, $reads_f, "$out_f.sam");
			}
			sepReads($reads_f, $out_f, "$out_f.sam");
			mkdir("$out_f.tmp");
			open(FILE, "<$out_f.summary");
			while(defined($line = <FILE>))
			{
				if ($line =~ /Bin \[$out_f.([0-9]+).fasta\] \(([0-9.]+)\)/)
				{
					$j = $1;
					writeLOG("Reassembling bin $j\n");
					$cmd = "$VELVETH $out_f.tmp $KMERLEN -short -fasta -interleaved $out_f.reads.$j 1>/dev/null 2>/dev/null";
					writeLOG("Command: $cmd\n");
					system($cmd);
					$i = $2 * 0.1;
					$cmd = "$VELVETG $out_f.tmp -exp_cov $2 -cov_cutoff $i 1>/dev/null 2>/dev/null";
					writeLOG("Command: $cmd\n");
					system($cmd);
					if (checkResult("$out_f.tmp\/contigs.fa") == -1)
					{
						writeLOG("Error occurs while reassembling Bin $out_f.$j.fasta.\nProgram stop.\n");
						closeLOG("$out_f.log");
						exit;
					}

					writeLOG("Replacing original contigs\n");
					open(REASSEM, "$out_f.tmp\/contigs.fa");
					while(defined($line = <REASSEM>))
					{
						chomp($line);
						if ($line =~ /^>/)
						{
							$currbin = $line;
							$lenhash{$currbin} = 0;
						}
						else
						{
							$lenhash{$currbin} += length($line);
						}
					}

					seek(REASSEM, 0, 0);
					open(REASSEM_TARGET, ">$out_f.$j.fasta");
					$k = 0;
					while(defined($line = <REASSEM>))
					{
						if ($line =~ /^>/)
						{
							chomp($line);
							if ($lenhash{$line} > $REASSEM_MIN_LEN)
							{
								$k = 1;
								print REASSEM_TARGET "$line\n";
							}
							else
							{
								$k = 0;
							}
						}
						elsif ($k == 1)
						{
							print REASSEM_TARGET $line;
						}
					}
					close(REASSEM);
					close(REASSEM_TARGET);

					#unlink("$out_f.reads.$1");
				}
			}
			close(FILE);
		}
		else
		{
			writeLOG("Cannot perform reassembly since reads file is not provided.\n");
		}
	}

	# Count marker genes for bins
	my $completearr;
	if ($REASSEMBLY == 1)
	{
		$cmd = "cat $out_f.*.fasta > $out_f.all.fasta";
		print "$cmd\n";
		system($cmd);
		getHMM("$out_f.all.fasta", "$contig_f.hmmout", $geneAAfile, "cuttc");
	}
	else
	{
		#getHMM($contig_f, "$contig_f.hmmout", "cuttc");
	}
	#"$contig_f.hmmout"
	$completearr = countmarker($essHMMfile, $out_f, $outname, $MARKERHMM, "$out_f.marker");
	$i = 0;
=cut
	foreach $currbin (@filearr)
	{
		($genomesize[$i], $gc[$i]) = getBinInfo($currbin);
		$i++;
	}
=cut
	# Re-read summary file and write other info: completeness, Genome size, GC
	open(FILE, "<$out_f.summary");
	open(SUMOUT, ">$out_f.summary.tmp");
	print SUMOUT "Bin name\tAbundance\tCompleteness\tGenome size\tGC content\n";
	$i = 0;
	while(defined($line = <FILE>))
	{
		if ($line =~ /^Bin \[$out_f.([0-9]+).fasta\] \(([0-9.]+)\)/)
		{
			printf SUMOUT "$outname.$1.fasta\t%0.2f\t%0.1f%%\t$genomesize[$i]\t%0.0f\n", $2, $$completearr[$i] * 100, $gc[$i];
			#printf SUMOUT "$out_f.$1.fasta\t";
			#printf SUMOUT "%0.2f\t", $2;
			#printf SUMOUT "%0.1f%%\t", $$completearr[$i] * 100;
			#printf SUMOUT "%d\t", $genomesize[$i];
			#printf SUMOUT "%0.0f\n", $gc[$i];
			$i++;
		}
	}
	close(FILE);
	close(SUMOUT);
	rename("$out_f.summary.tmp", "$out_f.summary");

	if ($RPLOT == 1)
	{
		$cmd = "$RSCRIPT $HEATMAP_R $out_f.marker $out_f.marker.pdf";
		print "$cmd\n";
		system($cmd);
	}

	if ($preserve == 0)
	{
		writeLOG("Deleting intermediate files.\n");
		unlink("$out_f.sam");
		$cmd = "rm -rf $out_f.tmp";
		system($cmd);
		unlink("$contig_f.prob_dist");
		#unlink("$contig_f.seed.m8");
		unlink("$contig_f.seed");
		if ($abund_f ne "" && $abund_f ne "$contig_f.abund")
		{
			unlink("$contig_f.abund");
		}
		else
		{
			rename("$contig_f.abund", "$out_f.abund");
		}
		unlink("$contig_f.hmmout");
		if (-e "$contig_f.frag.faa")
		{
			unlink("$contig_f.frag.faa");
		}
		if ($REASSEMBLY == 1)
		{
			unlink("$out_f.all.fasta");
			unlink("$out_f.all.fasta.frag.faa");
		}
	}
	unlink("$old_contig_name.1.bt2");
	unlink("$old_contig_name.2.bt2");
	unlink("$old_contig_name.3.bt2");
	unlink("$old_contig_name.4.bt2");
	unlink("$old_contig_name.rev.1.bt2");
	unlink("$old_contig_name.rev.2.bt2");
	if ($remodel_contig == 1)
	{
		unlink($contig_f);
	}

	# Write post instruction to users
	writeLOG("\n\n========== Job finished ==========\nYielded $maxnum bins for contig (scaffold) file $contig_name\n\n");
	writeLOG("Here are the output files for this run.\nPlease refer to the README file for further details.\n\n");
	$i = $maxnum;
	if ($i < 10)
	{
		$j = "00" . $i;
	}
	elsif ($i < 100)
	{
		$j = "0" . $i;
	}
	else
	{
		$j = $i;
	}
	writeLOG("Summary file: $out_f.summary\nMarker file: $out_f.marker\nBin files: $out_f.001.fasta - $out_f.$j.fasta\nUnbinned sequences: $out_f.noclass\n");
	if ($RPLOT == 1 && -e "$out_f.marker.pdf")
	{
		writeLOG("Marker plot: $out_f.marker.pdf\n");
	}
	writeLOG("\n\n========== Elapsed Time ==========\n");
	$endtime = time();
	$line = getElapsedTime($endtime - $starttime);
	writeLOG("$line\n");

	closeLOG("$out_f.log");
}

sub runBowtie2
{
	my $contig_f = $_[0];
	my $reads_f = $_[1];
	my $out_f = $_[2];
	my $cmd;
	my $isfastq;

	# Run bowtie2 to find abundance information
	$cmd = "$BOWTIE2BUILD $contig_f $contig_f 1>$out_f.bowtie2build.out 2>$out_f.bowtie2build.err";
	system($cmd);
	$isfastq = checkFastq($reads_f);
	if ($isfastq == 1)
	{
		$cmd = "$BOWTIE2 -p $THREADNUM -x $contig_f -U $reads_f -S $out_f 1>$out_f.bowtie2.out 2>$out_f.bowtie2.err";
	}
	else
	{
		$cmd = "$BOWTIE2 -f -p $THREADNUM -x $contig_f -U $reads_f -S $out_f 1>$out_f.bowtie2.out 2>$out_f.bowtie2.err";
	}
	system($cmd);
	if (checkResult("$out_f") == -1)
	{
		print "Error running Bowtie2. Bowtie2 log recorded in $out_f.bowtie2build.out, $out_f.bowtie2build.err, $out_f.bowtie2.out, and $out_f.bowtie2.err\nProgram stop.\n";
		exit(-1);
	}
	else
	{
		unlink("$out_f.bowtie2build.out");
		unlink("$out_f.bowtie2build.err");
		unlink("$out_f.bowtie2.out");
		unlink("$out_f.bowtie2.err");
	}
}

sub getHMM
{
	my $contig_f = $_[0];
	my $out_f = $_[1];
	my $geneAA = $_[2];
	my $cutmethod = $_[3];
	my $cmd;
	if (!-e $geneAA){ die "Can't file AA file $geneAA\n"; }
	#if (!(-e "$contig_f.frag.faa"))
	#{
#		print "Running FragGeneScan....\n";
#		$cmd = "$RUNFRAG -genome=$contig_f -out=$contig_f.frag -complete=0 -train=complete -thread $THREADNUM 1>$contig_f.frag.out 2>$contig_f.frag.err";
#		system($cmd);
#	}
	if (checkResult("$geneAA") == -1)
	{
		print "Error running FragGeneScan. Output recorded in $contig_f.frag.out and $contig_f.frag.err.\n";
		print "=== Please make sure that you are running FragGeneScan 1.18 or above. ===\n";
		print "=== There are known bugs in v1.17 and before that will crash FragGeneScan program. ===\n";
		exit(-1);
	}
	else
	{
		unlink("$contig_f.frag.out");
		unlink("$contig_f.frag.err");
	}

=cut # Prodigal gene prediction
	print "Running Prodigal....\n";
	$cmd = "$PRODIGAL -a $contig_f.faa -i $contig_f -m -o $contig_f.prodigal.tmp -p meta -q 1>$contig_f.prodigal.out 2>$contig_f.prodigal.err";
	system($cmd);
	if (checkResult("$contig_f.faa") == -1)
	{
		print "Error running Prodigal. Output recorded in $contig_f.prodigal.out and $contig_f.prodigal.err.\nProgram Stop.\n";
		exit(-1);
	}
	else
	{
		unlink("$contig_f.prodigal.out");
		unlink("$contig_f.prodigal.err");
		unlink("$contig_f.prodigal.tmp");
	}
=cut

	print "Running HMMER hmmsearch....\n";
	if ($MARKERNUM == 107)
	{
		#$cmd = "$HMMSEARCH --tblout $out_f --cut_tc --cpu $THREADNUM $MARKERHMM $contig_f.frag.faa 1>$out_f.out 2>$out_f.err";
		#if (defined($cutmethod) && $cutmethod eq "cuttc")
		#{
			$cmd = "$HMMSEARCH --domtblout $out_f --cut_tc --cpu $THREADNUM $MARKERHMM $geneAA 1>$out_f.out 2>$out_f.err";
		#}
		#else
		#{
		#	$cmd = "$HMMSEARCH --domtblout $out_f -E 1e-5 --cpu $THREADNUM $MARKERHMM $contig_f.frag.faa 1>$out_f.out 2>$out_f.err";
		#}
	}
	else
	{
		$cmd = "$HMMSEARCH --domtblout $out_f -E 1e-3 --cpu $THREADNUM $MARKERHMM $geneAA 1>$out_f.out 2>$out_f.err";
	}
	system($cmd);
	if (checkResult("$out_f") == -1)
	{
		print "Error running Hmmer3. Output recorded in $out_f.out and $out_f.err.\nProgram Stop.\n";
		exit(-1);
	}
	else
	{
		unlink("$out_f.out");
		unlink("$out_f.err");
	}
	#unlink("$contig_f.frag.faa");
	unlink("$contig_f.frag.ffn");
	unlink("$contig_f.frag.gff");
	unlink("$contig_f.frag.out");
}

sub checkResult
{
	my $s = -s $_[0];
	if (-e $_[0] && $s > 0)
	{
		return 0;
	}
	else
	{
		return -1;
	}
}

sub checkContig
{
	my $contig_f = $_[0];
	my $min_length = $_[1];
	my $below_f = $_[2];
	my $i;
	my $len;
	my $faline;
	my $header;
	my $seq;
	my $FA_LIMIT = 128;
	my $FASTA_LEN = 70;
	my $faout = "";
	open(FILE, "<$contig_f") || die "Cannot open file $contig_f\n";
	$i = 0;

	# Remodel input contig file
	seek(FILE, 0, 0);
	$header = "";
	$faout = $contig_f . ".tmp";
	open(FAOUT, ">$faout");
	open(FABELOW, ">$below_f");
	while(defined($faline = <FILE>))
	{
		chomp($faline);
		if ($faline =~ /^>/)
		{
			if ($header ne "" && ($len = length($seq)) >= $min_length)
			{
				print FAOUT "$header\n";
				$i = 0;
				while ($i < $len)
				{
					if ($i + $FASTA_LEN > $len)
					{
						print FAOUT substr($seq, $i);
						print FAOUT "\n";
						$i = $len;
					}
					else
					{
						print FAOUT substr($seq, $i, $FASTA_LEN);
						print FAOUT "\n";
						$i += $FASTA_LEN;
					}
				}
			}
			elsif ($header ne "" && ($len = length($seq)) < $min_length)
			{
				print FABELOW "$header\n";
				$i = 0;
				while ($i < $len)
				{
					if ($i + $FASTA_LEN > $len)
					{
						print FABELOW substr($seq, $i);
						print FABELOW "\n";
						$i = $len;
					}
					else
					{
						print FABELOW substr($seq, $i, $FASTA_LEN);
						print FABELOW "\n";
						$i += $FASTA_LEN;
					}
				}
			}
			$header = $faline;
			$seq = "";
		}
		else
		{
			$seq = $seq . $faline;
		}
	}
	$len = length($seq);
	if ($len >= $min_length)
	{
		print FAOUT "$header\n";
		$i = 0;
		while ($i < $len)
		{
			if ($i + $FASTA_LEN > $len)
			{
				print FAOUT substr($seq, $i);
				print FAOUT "\n";
				$i = $len;
			}
			else
			{
				print FAOUT substr($seq, $i, $FASTA_LEN);
				print FAOUT "\n";
				$i += $FASTA_LEN;
			}
		}
	}
	else
	{
		print FABELOW "$header\n";
		$i = 0;
		while ($i < $len)
		{
			if ($i + $FASTA_LEN > $len)
			{
				print FABELOW substr($seq, $i);
				print FABELOW "\n";
				$i = $len;
			}
			else
			{
				print FABELOW substr($seq, $i, $FASTA_LEN);
				print FABELOW "\n";
				$i += $FASTA_LEN;
			}
		}
	}
	close(FILE);
	close(FAOUT);
	close(FABELOW);

	return $faout;
}

# (Genome size, GC content) = getBinInfo(fasta file)
sub getBinInfo
{
	my $all = 0;
	my $gc = 0;
	my $i;
	my $line;
	open(FAFILE, "<$_[0]");
	while(defined($line = <FAFILE>))
	{
		if ($line !~ /^>/)
		{
			chomp($line);
			$i = ($line =~ tr/[cgCG]//);
			$gc += $i;
			$i = ($line =~ tr/[atcgATCG]//);
			$all += $i;
		}
	}
	close(FAFILE);

	if ($all == 0)
	{
		$gc = 0;
	}
	else
	{
		$gc = $gc / $all * 100;
	}
	return ($all, $gc);
}

# (Read and Check program directory)
sub checkProgram
{
	my $line;
	my $tmpstr;

	open(FILE, "<$Bin\/$SETTING_FILE");
	while(defined($line = <FILE>))
	{
		chomp($line);
		if ($line =~ /\[([A-Za-z0-9]+)\] ([A-Za-z0-9._\(\)\[\]\{\}\|\$\!\=\-\+\\\/]+)/)
		{
			if ($1 eq "FragGeneScan")
			{
				if (-d $2 && -e "$2\/$RUNFRAG")
				{
					$RUNFRAG = $2 . "\/" . $RUNFRAG;
				}
			}
			elsif ($1 eq "Bowtie2")
			{
				if (-d $2 && -e "$2\/$BOWTIE2" && -e "$2\/$BOWTIE2BUILD")
				{
					$BOWTIE2 = $2 . "\/" . $BOWTIE2;
					$BOWTIE2BUILD = $2 . "\/" . $BOWTIE2BUILD;
				}
			}
			elsif ($1 eq "HMMER3")
			{
				if (-d $2 && -e "$2\/$HMMSEARCH")
				{
					$HMMSEARCH = $2 . "\/" . $HMMSEARCH;
				}
			}
			elsif ($1 eq "Velvet")
			{
				if (-d $2 && -e "$2\/$VELVETH" && -e "$2\/$VELVETG")
				{
					$VELVETG = $2 . "\/" . $VELVETG;
					$VELVETH = $2 . "\/" . $VELVETH;
				}
			}
		}
	}
	close(FILE);
#I don't need these programs..
	# Check program
	# FragGeneScan
#	if (-d "tmp")	{		print "FragGeneScan cannot run when there is a subdirectory named \"tmp\". Please remove or rename this directory and try again.\n";		exit;	}
	$line = "$RUNFRAG 1>tmp 2>/dev/null";
	system($line);
	$tmpstr = "";
	open(FILE, "<tmp");
	while(<FILE>)
	{
		$tmpstr .= $_;
	}
#	if ($tmpstr !~ /FragGeneScan/)	{		print "Cannot run FragGeneScan. Please indicate the file directory in \'setting\' file.\n";		exit;	}

	# Bowtie2
	$line = "$BOWTIE2 1>/dev/null 2>tmp";
	system($line);
	$tmpstr = "";
#	open(FILE, "<tmp");	while(<FILE>)	{		$tmpstr .= $_;	}
#	if ($tmpstr !~ /bowtie/)	{		print "Cannot run Bowtie2. Please indicate the file directory in \'setting\' file.\n";		exit;	}

	# HMMER3
	$line = "$HMMSEARCH 1>tmp 2>/dev/null";
	system($line);
	$tmpstr = "";
	open(FILE, "<tmp");
	while(<FILE>)
	{
		$tmpstr .= $_;
	}
	if ($tmpstr !~ /hmmsearch/)
	{
		print "Cannot run HMMER3. Please indicate the file directory in \'setting\' file.\n";
		exit;
	}

	# velvet
	$line = "$VELVETH 1>tmp 2>/dev/null";
	system($line);
	$tmpstr = "";
	open(FILE, "<tmp");
	while(<FILE>)
	{
		$tmpstr .= $_;
	}
	if ($tmpstr !~ /velveth/)
	{
		print "Cannot run Velvet. Please indicate the file directory in \'setting\' file.\n";
		exit;
	}

	unlink("tmp");
}

sub checkFastq
{
	my $input_f = $_[0];
	my $fastq_line;
	my $fastq = -1;
	open(CHECK, "<$input_f") || die "Reads file [$input_f] not found. Check path again.\n";
	while(defined($fastq_line = <CHECK>))
	{
		chomp($fastq_line);
		if ($fastq_line =~ /^>/)
		{
			$fastq = 0;
		}
		elsif ($fastq_line =~ /^\@/)
		{
			$fastq = 1;
		}
		elsif ($fastq_line ne "")
		{
			print "Error header found in reads file.\n======\n$fastq_line\n======\nMaxBin stop.\n";
			exit;
		}
		if ($fastq != -1)
		{
			last;
		}
	}
	close(CHECK);
	return($fastq);
}

sub checkAbundFile
{
	my $input_f = $_[0];
	my $out_f = $_[1];
	my $line;

	if ($input_f eq $out_f)
	{
		return;
	}

	#Check header
	open(ABUNDIN, "<$input_f") || die "Cannot open abund file $input_f\n";
	open(ABUNDOUT, ">$out_f") || die "Cannot open tmp abund file $out_f\n";
	while(defined($line = <ABUNDIN>))
	{
		if ($line =~ /^>/)
		{
			print ABUNDOUT substr($line, 1);
		}
		else
		{
			print ABUNDOUT $line;
		}
	}
	close(ABUNDOUT);
	close(ABUNDIN);
}

sub getElapsedTime
{
	my $input = $_[0];
	my $hour;
	my $min;
	my $sec;
	my $str;

	$sec = $input % 60;
	$input = int($input / 60);
	$min = $input % 60;
	$input = int($input / 60);
	$hour = $input;

	$str = $hour . " hours " . $min . " minutes and " . $sec . " seconds.\n";
	return $str;
}

sub openLOG
{
	my $log_f = $_[0];
	open($LOGOUT, ">$log_f.tmp");
}

sub writeLOG
{
	my $msg = $_[0];
	print $msg;
	print $LOGOUT $msg;
}

sub closeLOG
{
	my $log_f = $_[0];
	close($LOGOUT);
	rename("$log_f.tmp", $log_f);
}

