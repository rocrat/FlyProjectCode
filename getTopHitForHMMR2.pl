#!/usr/bin/env perl
use warnings;
use strict;

#ensure the directory only contains the hmmr1 files you want to find the best hits in
#find each hmmr1 file
my @hmmrs = `ls /xdisk/rlapoint/ | grep '\.hmmr1'`;

foreach my $file (@hmmrs){
	if($file !~ /.+Rhmmr1$/){#in case there are some Rhmmrs lying about!
	chomp $file;
	print $file;
#read tab delimited blast output to find best hits for each query gene
	my $outfile = substr($file,0,-7);
	my @cols;
	my %hash = ();
	open (INFIL,$file) or die "Cannot open $file\n";

 
	while (my $line = <INFIL>){#get unique query sequence IDs
		if ($line !~ /^[\#]/){
			@cols = split /\s+/, $line;
			my $score = $cols[5];
			my $id = $cols[2]; #the query name
			if (!exists $hash{$id}){#if a new ID then assign the score to the hash key
				$hash{$id}[0] = $cols[2];
				$hash{$id}[1] = $cols[0];
				$hash{$id}[2] = $cols[4];
				$hash{$id}[3] = $cols[5];
			
			}elsif ($hash{$id}[3] < $score){#if an existing ID then check to see if score is higher
				$hash{$id}[0] = $cols[2];
				$hash{$id}[1] = $cols[0];
				$hash{$id}[2] = $cols[4];
				$hash{$id}[3] = $cols[5];
			}
		}
	}
	close INFIL;
	open (OUTFIL, ">$outfile.besthits") or die;
	my @loci = keys %hash;
	foreach my $id (@loci){
		print OUTFIL "$hash{$id}[0];$hash{$id}[1];$hash{$id}[2];$hash{$id}[3]\n";
	}	
	close OUTFIL;
}
}

