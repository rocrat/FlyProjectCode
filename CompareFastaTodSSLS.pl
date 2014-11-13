#!/usr/bin/env perl
use warnings;
use strict;

my $infile1 = $ARGV[0];#Fasta file input
my $infile2 = $ARGV[1];#dSSLS file input
my $outfile = $ARGV[2];#file to place results
open (INFIL1,$infile1);
my $numgenes = 0;
my $numgenes2 = 0;
my $largest = 0;
my $smallest = 9999999;
my %fastahash;
my $duplicates=0;

while(my $line = <INFIL1>){
	if ($line =~ /\>(ENSAPLG)(\d+)\|(ENSAPLT\d+)/){
		$numgenes++;
		my $linenum = $2;
		my $id = "$1$2_$3";
		if ($linenum > $largest){
			$largest = $linenum;
		}
		if ($linenum < $smallest){
			$smallest = $linenum;
		}
		if (!exists $fastahash{$id}){
			$numgenes2++;
			$fastahash{$id} = 1;
		}elsif(exists $fastahash{$id}){
			$duplicates++;
		}
	}
}
close INFIL1;
open (INFIL2,$infile2);

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
		if (exists $fastahash{$id}){
			$matched++;#count up the files that match
		}elsif(!exists $fastahash{$id}){		
			$ecnt++;
			#print ERROR "$cols[0],$cols[1]";#add errors to error file with score
			#$errhash{$id} = $cols[1];
		if ($cols[1]>0){
			$pscnt++;
		}
		}
	}
}
close INFIL2;
open (OUTFIL,">$outfile.txt");
print OUTFIL "The fasta file has $numgenes ($numgenes2) genes.\n";
print OUTFIL "The smallest numbered gene is ENSFALG$smallest\n";
print OUTFIL "The largest numbered gene is ENSFALG$largest\n";
print OUTFIL "There were $duplicates duplicates in the fasta file\n";
print OUTFIL "There were $count files in the first dSSLS file entered,\n";
print OUTFIL "of those, $matched were also in the gene list.\n";
print OUTFIL "Therefore, $ecnt were not included in the second analysis.\n";
close OUTFIL;

print "The fasta file has $numgenes ($numgenes2) genes.\n";
print "The smallest numbered gene is ENSFALG$smallest\n";
print "The largest numbered gene is ENSFALG$largest\n";
print "There were $duplicates duplicates in the fasta file\n";
print "There were $count files in the first dSSLS file entered,\n";
print "of those, $matched were also in the gene list.\n";
print "Therefore, $ecnt were not included in the second analysis.\n";
