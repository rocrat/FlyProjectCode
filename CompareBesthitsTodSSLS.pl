#!/usr/bin/env perl
use warnings;
use strict;

my $infile1 = $ARGV[0];#RBesthits file input
my $infile2 = $ARGV[1];#dSSLS file input
my $outfile = substr($infile1,0,10);#file to place results
open (INFIL1,$infile1);
my $numgenes = 0;
my $numgenes2 = 0;
my %Rhash;
my $duplicates=0;

while(my $line = <INFIL1>){#make hash of duck gene id's in the Rbesthits file
	my @cols = split /,/, $line;
	if ($cols[0] =~ /(ENSAPLG)(\d+)\|(ENSAPLT\d+)/){
		$numgenes++;
		my $id = "$1$2_$3";
		
		if (!exists $Rhash{$id}){
			$numgenes2++;
			$Rhash{$id} = 1;
		}elsif(exists $Rhash{$id}){
			$duplicates++;
		}
	}
	if ($cols[1] =~ /(ENSAPLG)(\d+)\|(ENSAPLT\d+)/){
		$numgenes++;
		my $id = "$1$2_$3";
		
		if (!exists $Rhash{$id}){
			$numgenes2++;
			$Rhash{$id} = 1;
		}elsif(exists $Rhash{$id}){
			$duplicates++;
		}
	}
}

close INFIL1;
open (INFIL2,$infile2);
open (ERROR,">$infile1.errfile.csv");
my @errorarray = ();
my %errhash =();
my $count = 0;
my $ecnt=0;
my $matched = 0;
my $pscnt = 0;
while (my $line = <INFIL2>){#Check the hash against each line of the old dSSLS file
	if ($line !~ /^[geneID]/){#skip the header
		$count++;
		my @cols = split /,/, $line;
		my $id = $cols[0];		
		if (exists $Rhash{$id}){
			$matched++;#count up the files that match
		}elsif(!exists $Rhash{$id}){		
			$ecnt++;
			print ERROR "$cols[0],$cols[1]";#add errors to error file with score
			$errhash{$id} = $cols[1];
		if ($cols[1]>0){
			$pscnt++;
		}
		}
	}
}
close INFIL2;
close ERROR;
open (OUTFIL,">$outfile.compareBesthitsResults.txt");
print OUTFIL "The Rbesthits file has $numgenes ($numgenes2) genes.\n";
print OUTFIL "There were $duplicates duplicates in the Rbesthits file\n";
print OUTFIL "There were $count files in the first dSSLS file entered,\n";
print OUTFIL "of those, $matched were also in the gene list.\n";
print OUTFIL "Therefore, $ecnt were not included in the second analysis ($pscnt pos).\n";
print OUTFIL "These genes and their scores are found in $infile1.errfile.csv/n";
close OUTFIL;

print "The Rbesthits file has $numgenes ($numgenes2) genes.\n";
print "There were $duplicates duplicates in the Rbesthits file\n";
print "There were $count files in the first dSSLS file entered,\n";
print "of those, $matched were also in the gene list.\n";
print "Therefore, $ecnt were not included in the second analysis ($pscnt pos).\n";
print "These genes and their scores are found in $infile1.errfile.csv\n";
