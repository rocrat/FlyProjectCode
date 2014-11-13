#!/usr/bin/env perl
use warnings;
use strict;

my @rhmmrs = `ls /xdisk/rlapoint/ | grep '\.Rhmmr1'`;
#find each Rhmmr1 file
foreach my $file (@rhmmrs){
	chomp $file;
	print "$file\n";
	#read tab delimited blast output to find best hits for each query gene
	
	my $outfile = substr($file,0,-7);
	my @cols;
	my %hash = ();
	open (INFIL,$file) or die "Cannot open $file\n";

 
	while (my $line = <INFIL>){#get unique query sequence IDs
		if ($line !~ /^[\#]/){
			@cols = split /\s+/, $line;
			my $score = $cols[5];
			my $id = $cols[0]; #the query name
			if (!exists $hash{$id}){#if a new ID then assign the score to the hash key
				$hash{$id}[0] = $cols[0];
				$hash{$id}[1] = $cols[2];
				$hash{$id}[2] = $cols[4];
				$hash{$id}[3] = $cols[5];
			
			}elsif ($hash{$id}[3] < $score){#if an existing ID then check to see if score is higher
				$hash{$id}[0] = $cols[0];
				$hash{$id}[1] = $cols[2];
				$hash{$id}[2] = $cols[4];
				$hash{$id}[3] = $cols[5];
			}
		}
	}
	close INFIL;
	open (OUTFIL, ">$outfile.Rbesthits") or die;
	my @loci = keys %hash;
	foreach my $id (@loci){
		print OUTFIL "$hash{$id}[0];$hash{$id}[1];$hash{$id}[2];$hash{$id}[3]\n";
	}	
	close OUTFIL;
}


