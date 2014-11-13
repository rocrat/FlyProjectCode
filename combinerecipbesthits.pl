#!/usr/bin/env perl
use warnings;
use strict;
#This script is designed to combine Rbesthits files for multiple species into a single file which lists all the related gene IDs on each line

my @allfiles = `ls | grep '\.Recipbest'`; #reads in all Recipbest files
chomp @allfiles;
my $n = $#allfiles; #get the number of files
#my $outname;
print $n+1,"\n";
#print "Please enter an output file name:  \n";
#chomp($outname = <STDIN>);
my $outname = "RecipBesthitsV3";

sub makehash {#make a hash from an individual file using the original ref spp which should be listed first
	my %hash = ();
	my $file = shift;
	open (INFIL,"$file") or die;
	my @cols;
		while (my $line = <INFIL>){
			@cols = split /;/, $line;
			my $id = $cols[0];
			
			
			if (!exists $hash{$id}){#if a new ID then assign the score to the hash key
				$hash{$id} = $cols[1];#the matched id is assigned as the value and the refspp is assigned as the key
				
			}
		}
	#print $#cols,"\n";
	close INFIL;
	print "$file\n";
	return %hash;
}

my $count=0;

my %hash1 = &makehash($allfiles[0]);

my %hash2 = &makehash($allfiles[1]);
my %hash3 = &makehash($allfiles[2]);
my %hash4 = &makehash($allfiles[3]);
my %hash5 = &makehash($allfiles[4]);
my %hash6 = &makehash($allfiles[5]);
my %hash7 = &makehash($allfiles[6]);
my %hash8 = &makehash($allfiles[7]);

my @hasharray = ({%hash1},{%hash2},{%hash3},{%hash4},{%hash5},{%hash6},{%hash7},{%hash8});
my $numfind = $n-1;

my %idhash;
foreach my $hashref (@hasharray){#get all the id's that are unique to all of the sequence files
	my @ids = keys %$hashref;
	foreach my $id (@ids){
		if (!exists $idhash{$id}){#make a unique hash entry for each ref id in each hash
			$idhash{id} = 1;#make one so if in all files then id will have a 7
		}
		$idhash{$id}++;#add one for each time the id is encountered in a file (should only be max 1 time per file)
		#print $idhash{$id};
	}
}



my @allIncluded = grep { $idhash{$_} == $n+1 } keys %idhash; #grab every key with a value of $n indicating it was found $n times (once in each of seven Rbesthits files)

open (OUTFIL,">RecipBesthitsV3.list");
foreach my $every (@allIncluded){#loop through each id that was found in all  species and print the matched genes
	print OUTFIL "$every;$hasharray[0]{$every};$hasharray[1]{$every};$hasharray[2]{$every};$hasharray[3]{$every};$hasharray[4]{$every};$hasharray[5]{$every};$hasharray[6]{$every};$hasharray[7]{$every}\n";
}
close OUTFIL;


