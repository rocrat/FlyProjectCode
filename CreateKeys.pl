#! usr/bin/env perl

use strict;
use warnings;
#this script reads in fasta files from a folder and creates two key files one with unique spp count and one with total numbers of each spp


my $folder = $ARGV[0];

my @files = `ls $folder | grep '\.fasta'` ;
open(SKEY,">UniqueSPPKeyV4.key");
open(NKEY,">SppNumKey.key");
foreach my $file (@files){
	chomp($file);
	my %spphsh;
	my $genecount;
	open(FIL,"./$folder/$file");
	while (my $line = <FIL>){
		chomp($line);
		if ($line =~ /^\>/){
			$genecount++;
			my @cols = split /,/, $line;
			if(!exists $spphsh{$cols[1]}){
				$spphsh{$cols[1]} = 1;
			}else{
				$spphsh{$cols[1]}++;
			}
		}
	}
	close FIL;
	my @uniquespp = keys %spphsh;
	my $numspp = scalar @uniquespp;
	print SKEY $file, ";",$numspp,";",join(";",@uniquespp),"\n";
	print NKEY $file,";$genecount;";
	foreach my $spp (@uniquespp){
		print NKEY "$spp,$spphsh{$spp};";
	}
	print NKEY "\n";
} 
close SKEY;
close NKEY;
