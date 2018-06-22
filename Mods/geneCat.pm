package Mods::geneCat;
use warnings;
use strict;
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(  systemW );

use Exporter qw(import);
our @EXPORT_OK = qw(readSam rewriteFastaHdIdx readGeneIdx);



sub readGeneIdx($){
	my ($in) = @_;
	my %ret; my $gCnt=0;
	open I,"<$in" or die "Can't read gene index file $in\n";
	while(my $line=<I>){
		next if ($line =~ m/#/);
		chomp $line;
		my @spl = split(/\t/,$line);
		$ret{$spl[2]} = $spl[0];
		$gCnt++;
		#print $spl[2] ." ". $spl[0]."\n";
	}
	close I;
	print "Gene Index read\n";
	return (\%ret,$gCnt);
}


