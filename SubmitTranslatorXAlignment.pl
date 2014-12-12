#!/usr/bin/env perl
use warnings;
use strict;

#This script creates a PBS script to align all the sequences in a particular file folder using TranslatorX

my $folder = $ARGV[0]; #the first argument is the path of the folder where the fasta files live
my @files = `ls $folder | grep '.\.fasta'`;
my $outfolder ="/xdisk/rlapoint/TranslatorClustAlignments";
open(OUT,">submitClustAlignments.pbs");
print OUT "#!/bin/bash\n";
print OUT "#PBS -N TransAlignCl\n";
print OUT "#PBS -W group_list=whiteman\n";
print OUT "#PBS -l jobtype=small_mpi\n";
print OUT "#PBS -l select=1:ncpus=12:mem=24000mb\n";
print OUT "#PBS -l place=pack:shared\n";
print OUT "#PBS -l pvmem=25598mb\n";
print OUT "#PBS -l cput=240:00:0\n";
print OUT "#PBS -l walltime=20:00:0\n";
print OUT "#PBS -q standard\n";
print OUT "#PBS -M dlaroche\@email.arizona.edu\n\n";
print OUT "cd \$PBS_O_WORKDIR # this cd's into the directory you submitted the script from\n";
print OUT "module load muscle\n";
	
for my $file (@files){
	chomp($file);
	my $outfile = substr($file,0,12);
	print OUT "perl TranslatorX.pl -i $folder/$file -o $outfolder/$outfile -p M -t F -g '-b1=5 -p=y -p=s'\n";
}

print OUT "wait";
close(OUT);

	
	