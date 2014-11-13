#!/usr/bin/env perl

#This script is designed to combine besthits files for the same species pair to get reciprocal best blasthits

my @allfiles = @ARGV;
my $file1 = $ARGV[0]; #1st besthits  file
my $file2 = $ARGV[1]; #second besthits file
my $outname;

$outname = substr($file1,0,-9); #get file output name

sub makehash {#make a hash from the first file with the first gene id as the key and the second as the value
	my %hash = ();
	my $file = shift;
	open (INFIL,"$file");
		while (my $line = <INFIL>){
			@cols = split /;/, $line;
			my $id = $cols[0];
			my $otherid = $cols[1];
			if (!exists $hash{$id}){#if a new ID then assign the score to the hash key
				$hash{$id} = $otherid;#pair the ids in the hash table
				
			}
		}
	close INFIL;
	#print "$file\n";
	return %hash;
}
my %hash = &makehash($file1);

open (INFIL,"$file2");
open (OUTFIL,">$outname.Recipbest");
open (ERRFIL,">$outname.Recipdiff");

while (my $line = <INFIL>){
	my @cols = split /;/, $line;
	my $lineid = $cols[1];#the species orders are identical in the second besthits file so this is the id of the species ot in the hash
	my $hashid = $cols[0];#the id of the hash key

	if ($hash{$hashid} eq $lineid){
		print OUTFIL $line;
	}else{
		print ERRFIL $line;
	}
}

close INFIL;
close OUTFIL;
close ERRFIL;
