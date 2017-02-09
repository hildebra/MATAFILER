#!/usr/bin/perl -w
use strict;

foreach my $file (@ARGV) {
open (FILE, $file);

my $output_file = "$file.decision.txt";
open (OUTPUT, ">$output_file");
		
while (<FILE>) {
chomp;
	$_ =~ s/^\s+//;
	$_ =~ s/\s+$//;

	my ($seqname, $euk, $bac, $blast) = split (/\t/);

	my ($kingdom,$phylum,$class,$order,$family,$genus,$species) = split(/;/, $blast);
	my ($assem, $tax, $taxdesc, $call, $calldesc);

#nothing matched
	if ($euk eq "nocall" & $bac eq "nocall" & $blast eq "nocall") {
	$assem = $seqname;
	$tax =  "nocall"; 
	$taxdesc =  "unambig_nocall";
	$call =  "nocall";
	$calldesc =  "unambig_nocall";
	}

#no gene calls, blast exists	
	elsif ($euk eq "nocall" & $bac eq "nocall" & $blast ne "nocall") {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc = "unambig_any";
	$call =  "nocall";
	$calldesc =  "unambig_nocall";
	}


#unambig bact call, blast = bact
	elsif ($euk eq "nocall" & $bac ne "nocall" & $blast =~ m/Bacteria/) {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc =  "unambig_bac";
	$call = $bac;
	$calldesc =  "unambig_bac";
	}

#unambig bact call, blast not bact
	elsif ($euk eq "nocall" & $bac ne "nocall" & $blast !~ m/Bacteria/) {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc =  "ambig_notbac";
	$call =  $bac;
	$calldesc =  "unambig_bac";
	}

#double call, blast = bact
	elsif ($euk ne "nocall" & $bac ne "nocall" & $blast =~ m/Bacteria/) {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc =  "ambig_bac";
	$call =  $bac;
	$calldesc =  "ambig_2call";
	}

#double call, blast = euk
	elsif ($euk ne "nocall" & $bac ne "nocall" & $blast =~ m/Eukaryota/) {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc =  "ambig_euk";
	$call =  $euk;
	$calldesc =  "ambig_2call";
	}

#double call, blast not bact or euk
	elsif ($euk ne "nocall" & $bac ne "nocall" & $blast !~ m/Eukaryota/ & $blast !~ m/Bacteria/) {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc =  "ambig_any";
	$call =  $bac;
	$calldesc =  "ambig_2call";
	}

#unambig euk call, blast is eukaryota
	elsif ($euk ne "nocall" & $bac eq "nocall" & $blast =~ m/Eukaryota/) {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc =  "unambig_euk";
	$call =  $euk;
	$calldesc =  "unambig_euk";
	}

#unambig euk call, blast = bact
	elsif ($euk ne "nocall" & $bac eq "nocall" & $blast =~ m/Bacteria/) {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc =  "ambig_bac";
	$call =  $euk;
	$calldesc =  "unambig_euk";
	}

#unambig euk call, blast neither bact nor eukaryota
	elsif ($euk eq "nocall" & $bac eq "nocall" & $blast !~ m/Eukaryota/ & $blast !~ m/Bacteria/) {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc =  "ambig_any";
	$call =  $euk;
	$calldesc =  "unambig_euk";
	}

#unambig euk, no blast call
	elsif ($euk ne "nocall" & $bac eq "nocall" & $blast eq "nocall") {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc =  "ambig_any";
	$call =  $euk;
	$calldesc =  "unambig_euk";
	}

#unambig euk, no blast call is not eukaryota or bacteria
	elsif ($euk ne "nocall" & $bac eq "nocall"  & $blast !~ m/Eukaryota/ & $blast !~ m/Bacteria/) {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc =  "ambig_any";
	$call =  $euk;
	$calldesc =  "unambig_euk";
	}
	
#everything else (hopefully nothing)
	else {
	$assem = $seqname;
	$tax =  $kingdom; 
	$taxdesc = "unk";
	$call =  join(":", $bac, $euk);
	$calldesc =  "unk";
	}

print OUTPUT "$assem\t$tax\t$taxdesc\t$call\t$calldesc\n";
}
}
exit; 
