#!/usr/bin/env perl
use warnings;
use strict;

# This program checks convergence between the pair of chains produced by phylobayes for each significant gene 



my $file = $ARGV[0];#open dSSLS file to find the significant genes
if (!defined $file) {
  die "Usage: $0 dSSLSfile\n";
}
my $count = 0;
print "Checking convergence between paired chains...\n";
my $ext1 = ".phy_1";
my $ext2 = ".phy_2";
my $extout = ".check";
open (INFIL,$file);

while( my $line = <INFIL>){#read in dSSLS file
		my @cols = split /\s+/, $line;
		my $fileid = $cols[0];
		my $pvalue = $cols[1];
		
		if ($pvalue < 0.05){#for each file id with dvalue>0 submit bpcomp and trace comp
			$count++; 
			my $chainid1 = $fileid.$ext1;
			my $chainid2 = $fileid.$ext2;
			my $outfile = $fileid.$extout;
			print "Checking convergence of $fileid \n";
			my $bpcmd = "bpcomp -x 5000 1 -o $outfile  /xdisk/rlapoint/SimsV3/$chainid1 /xdisk/rlapoint/SimsV3/$chainid2";
			my $tracecmd = "tracecomp -x 5000  1 -o $outfile  /xdisk/rlapoint/SimsV3/$chainid1 /xdisk/rlapoint/SimsV3/$chainid2";
			system($bpcmd);
			system($tracecmd);
		}
}
close INFIL;
print "Finished cehcking $count chain pairs\n\n";

