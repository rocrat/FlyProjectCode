#!usr/bin/env perl

#This script takes one HMMER file as a base and then appends each additional HMMER file onto the base file

use strict;
use warnings;
my $filedir = $ARGV[0];#get file directory of HMMER Files from command line
my @hfiles = `ls $filedir | grep '.\.hmmr1'`;#get the names of all the files in the directory

#initiate OUTFILput and make header
open(OUTFIL,">CombinedHMMER.csv");
print OUTFIL "id1;id2;evalue\n";

foreach my $file (@hfiles){
	open(INFIL,$file);
	my $qspp = substr($file,0,4);
	my $tspp = substr($file,5,4);
	while (my $line = <INFIL>){
		if ($line !~ /^[\#]/){
			my @cols = split /\s+/, $line;
			if($cols[4]<0.00001){#limit to strong hits, I would like to leave this out in the future and let the weighted clustering algorithm do the work
				print OUTFIL "$tspp.$cols[0];$qspp.$cols[2];$cols[4]\n";
			}
		}
	}
	close(INFIL);
}
close(OUTFIL);
		
