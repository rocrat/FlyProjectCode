#!/usr/bin/env perl
use warnings;
use strict;
#This script is designed to combine Rbesthits files for multiple species into a single file which lists all the related gene IDs on each line

my @allfiles = `ls | grep '\.Recipbest'`; #reads in all Recipbest files
chomp @allfiles;
my $n = $#allfiles; #get the number of files
#my $outname;
print $n+1," .Recipbest files were found.\n";
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
my @hasharray;
foreach my $file (@allfiles){
	
	if($file !~ /Snoth/){
		my %hash = &makehash($file);
		
		push (@hasharray, {%hash});
	}
}



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

my $compnum = $#hasharray+1;
print $compnum," species were added to each line of the list file. \n";
my @allIncluded;
foreach my $id (keys %idhash){ #grab every key with a value of $n indicating it was found $n times (once in each of seven Rbesthits files)
	if ($idhash{$id} == $compnum){
		push (@allIncluded, $id);
		
	}
}

open (OUTFIL,">RecipBesthitsV3NoSnot.list");
foreach my $every (@allIncluded){#loop through each id that was found in all  species and print the matched genes
	
	print OUTFIL "$every;";

	foreach my $href (@hasharray){#loop through each hash in the array and print the 
		print OUTFIL $href->{$every},";";
	}

	print OUTFIL "\n";
}
close OUTFIL;


