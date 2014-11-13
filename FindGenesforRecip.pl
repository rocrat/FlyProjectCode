#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
#use Bio::Index::Fasta;
use Bio::DB::Fasta;
#use Bio::Index::GenBank;
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";
$ENV{BIOPERL_INDEX} = ".";

#my $dmel = Bio::DB::Fasta->new("DmelanogasterPEP.fasta");#create Dmel index

my $Best = $ARGV[0];#the best hits file
my $outfile = substr($Best,0,-9);#same as the best hits file bit without the 
#print $outfile,"\n";#check outfile
sub get_id {#get the ID of each gene
   my $header = shift;
   my $id;
   if ($header =~ /([A-Z].{13})/){
   $id = $1;
   }
   return $id;
}

open (INFIL,"$Best") or die "could not open best hits file\n";
open (OUTFIL,">$outfile.rquery.fasta");

while (my $line = <INFIL>){
	my @cols = split /,/, $line;
	my $sppid = $cols[1];#the reference gene id
	my $refsppseq = $dmel->get_Seq_by_id($sppid);
	print OUTFIL ">$sppid\n",$refsppseq->seq,"\n";
}
close INFIL;
close OUTFIL;
	
	
	

