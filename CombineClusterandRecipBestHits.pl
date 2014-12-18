#! usr/bin/env perl

use strict;
use warnings;

#This script is designed to combine a clustering result with a reciprocal best hit .list file

my $listfile = $ARGV[0];
my $clustfile = $ARGV[1];

open(LIST,$listfile);
open(CLUST,$clustfile);
open(OUT,">CombinedRecipClust.list");
open(KEY,">CombinedRecipClust.key");
#make hash for cluster file
my %clusthsh;
my $n=0;
while (my $line=<CLUST>){
	chomp($line);
	my $num;
	my $spp;
	#file starts with row numbers which are all unique
	if($line =~ /(\d+),(.+)/){
		$num = $1;
	}
	my $oldline = $line;
	$line =~ s/\d+,//;#remove row number from line
	$line =~ s/"//g;#remove quatation marks
	my @cols = split /;/, $line;
	foreach my $col (@cols){#this will make a key->field pair for each gene id in the cluster so the cluster will be represented as many times as there are genes, this is necessary because we don't know which gene will be in the recip file but we will need to account for this when we make the combined file
		if(!exists $clusthsh{$col}){
			$clusthsh{$col} = $line;
		}else{
			#$clusthsh{$col} = $clusthsh{$col}.";".$line;
			$n++;
		}
	}	
}
#go through each line of the list file and check for a matching cluster 
while (my $line = <LIST>){
	chomp($line);
	my @clusts;#shell for all clusters found
	my @cols = split /;/, $line;
	foreach my $geneid (@cols){#go through each gene id and check for a matching cluster
		if(exists $clusthsh{$geneid}){
			push @clusts, $clusthsh{$geneid};#if there is a cluster with the gene id add it to the clusts array
		}
	}
	#get just the unique clusters in the array since some may have been added more than once
	my @unique;
	my %seen;
	foreach my $value (@clusts) {
		if (! $seen{$value}++ ) {
			push @unique, $value;
		}
	}
	
	#create array of geneIDs to find only unique IDs before printing
	my @outgenes = split /;/, $line;
	
	foreach my $spp (@unique){
		my @cols = split /;/, $spp;
		foreach my $sp (@cols){
			push @outgenes, $sp;#print each unique cluster associated with that line
		}
	}
	#find only the unique genes
	my @uniquegenes;
	my %seenids;
	foreach my $id (@outgenes) {
		if (! $seenids{$id}++ ) {
			push @uniquegenes, $id;
		}
	}
	print scalar @outgenes,"\t",scalar @uniquegenes,"\n";
	print OUT join(";",@uniquegenes);#print the line from the recipbesthits.list file
	print KEY scalar @uniquegenes, "\n";
	print OUT "\n";
}

print $n;

close(LIST);
close(CLUST);
close(OUT);
close(KEY);