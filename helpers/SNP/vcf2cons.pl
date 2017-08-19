#!/usr/bin/perl -w
use strict;
use warnings;

sub v2q_post_process;
my $indelWin=5; 

my $line1tmp = <>;
my $statfile="";
if ($line1tmp =~ m/#depthStat (\S+)/){
	$statfile = $1;
}

my %het = (AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y',
             GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K');
my ($last_chr, $seq, $qual, $last_pos, @gaps, @spos);
$last_pos=0;$last_chr="";
my $_Q=20; my $_d=1; my $_D=30000; 
my $bcnts=0; my $ccnts =0;
my $lbcnts=0; my $lccnts =0;
my $lcnt=0;
my $chromL=0;
#depth related stats
my %depthStat;
while (<>) {
   $lcnt++;
   next if (/^#/);
    my @t = split;
    if ($last_chr eq ""){$last_chr = $t[0];$last_chr =~ m/L=(\d+)=/;  $chromL = $1; }#die "$chromL\n";}
	die @t if ($last_pos == $t[1] );
	
 	if ( $last_chr ne $t[0]) {
		#print "$last_pos > ".length ($seq)."\n";
		$last_pos = $chromL+1 if ($chromL > $last_pos);
		#die "$last_pos $chromL \n";
		  if ($last_pos > length ($seq)){
			my $ext = $last_pos - length($seq) - 1;
			$seq .= 'n' x ($ext);
		  }
		  &v2q_post_process($last_chr, \$seq, \$qual, \@gaps, $indelWin,$lbcnts,$lccnts,\@spos) if ($last_chr);
		  ($last_chr, $last_pos) = ($t[0], 0);
		  $seq = $qual = '';@spos=();
		  @gaps = ();$lbcnts=0; $lccnts =0;
		  $last_chr =~ m/L=(\d+)=/;  $chromL = $1;
    }
	if ($t[1] - $last_pos > 1) {
      $seq .= 'n' x ($t[1] - $last_pos - 1);
      #$qual .= '!' x ($t[1] - $last_pos - 1);
    }
    die("[vcf2cons] unsorted input\n$t[1] - $last_pos\non line $lcnt\n") if ($t[1] - $last_pos < 0);
	if (length($t[3]) == 1 && $t[7] !~ /INDEL/ && $t[4] =~ /^([A-Za-z.])(,[A-Za-z])*$/) { # a SNP or reference
	  my ($ref, $alt) = ($t[3], $1);
	  my ($b, $q);
	  $q=0;
	  $b = $ref;
      #die "@t\n$alt  $ref\n";
	#  if ($t[7] =~ /FQ=(-?[\d\.]+)/){
	#	$q = $1;
	#  } 
    #  if ($q < 0) {
    #    $_ = ($t[7] =~ /AF1=([\d\.]+)/)? $1 : 0;
    #    $b = ($_ < .5 || $alt eq '.')? $ref : $alt;
    #    $q = -$q;
    #  } else {
    #    $b = $het{"$ref$alt"};
    #    $b ||= 'N';
    #  }
	  my $dep=0;
	  if ($t[7] =~ /DP=(\d+)/){$dep= $1;}
	  
	  if ($dep<=1){$b = "N";}
	  if ($alt =~ m/,/ || length($alt)>1){
		die "@t\n";
	  }
	  
	  if ($t[7] =~ /AF=(-?\d[\.\d+]?)/){
		$q = $1;
		if ($q > 0.501 && $alt ne '.' ){
			if ( $dep>1 ) { 
				#this is an alternate allele
				$b = $alt;$ccnts++; $lccnts++;push(@spos,$t[1]);
				$depthStat{$dep}{alt}++;$depthStat{$dep}{norm}--;
			}
		}
	  }
	  $depthStat{$dep}{norm}++;
	  #die "$q $b\n";
      $b = lc($b);
      $b = uc($b) if (($t[7] =~ /QA=(\d+)/ && $1 >= $_Q && $1 >= $_d && $1 <= $_D) );
     $seq .= $b;
	  #die "yes";
      #$q = int($q + 33 + .499);
      #$q = chr($q <= 126? $q : 126);
      #$qual .= $q;
	  $bcnts++;$lbcnts++; 
    } elsif ($t[4] ne '.') { # an INDEL
	  die "@t\n";
      push(@gaps, [$t[1], length($t[3])]);
    }
	#print "$last_pos\n ";
    $last_pos = $t[1];
  }


$last_pos = $chromL+1 if ($chromL > $last_pos);
if ($last_pos > length ($seq)){
	my $ext = $last_pos - length($seq) - 1;
	$seq .= 'n' x ($ext);
}
&v2q_post_process($last_chr, \$seq, \$qual, \@gaps, $indelWin,$lbcnts,$lccnts,\@spos) if ($lcnt>0);
#print STDERR "$bcnts $ccnts\n";
#print depth stat instead
my $dsStr=""; my @cumSum;
my @deps = sort {$a <=> $b} keys %depthStat;
for (my $i=0; $i<$#deps;$i++){
	if (exists($depthStat{$i}){
		push (@cumSum,$#cumSum+$depthStat{$i}{norm}+$depthStat{$i}{alt});
	} else {
		push(@cumSum,$#cumSum);
	}
}
die "@cumSum\n";
foreach my $dep (@deps){
	$dsStr .= $dep."\t";
	if (exists($depthStat{$dep}{alt})){
		$dsStr.=$depthStat{$dep}{alt}."\t";
	} else{
		$dsStr .= "0\t";
	}
	if (exists($depthStat{$dep}{norm})){
		$dsStr.=$depthStat{$dep}{norm}."\n";
	} else{
		$dsStr .= "0\n";
	}
	
}
print STDERR $dsStr;


exit;
  
sub v2q_post_process {
  my ($chr, $seq, $qual, $gaps, $l,$reports,$replaces, $ARpos) = @_;
 # print $chr." ".length($$seq)." ".@gaps."\n";
  my @pos = @{$ARpos};
  for my $g (@$gaps) {
	#print "@{$g} ".length($$seq)."\n";
    my $beg = $g->[0] > $l? $g->[0] - $l : 0;
    my $end = $g->[0] + $g->[1] + $l;
    $end = length($$seq) if ($end > length($$seq));
    substr($$seq, $beg, $end - $beg) = lc(substr($$seq, $beg, $end - $beg));
  }
  
  print ">$chr COV=$reports REPL=$replaces POS=".join(",",@pos)."\n$$seq\n"; #&v2q_print_str($seq);
  #print "+\n"; &v2q_print_str($qual);
}