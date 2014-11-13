#!/usr/bin/env perl
use warnings;
use strict;

#use this program to submit the reciprocal hmmers.  
my @allfiles = `ls /xdisk/rlapoint | grep '\.rquery\.fasta'`; #get list of rquery databases
my $refspp = $ARGV[0]; #The orginial query database from the first HMMER

my $count=0;

foreach my $file (@allfiles){
	chomp $file;
	my $query = $file;#these are the sequences of the original target species with the best hits in the first HMMER
	
	
	#my $newdata;
	#foreach my $data (@pepfiles){
	#	chomp($data);
	#	if ($data !~ /Dmelanogaster/){
	#		if ($data =~ /$datapre/){
	#			$newdata = $data;
	#		}
#
#		}
#	}
#	print $newdata,"\n";
		
		my $hmmrid = substr($file,0,-13);
		my $subid = substr($hmmrid,5);
		open (NEWSub,">/xdisk/rlapoint/RHMMR2/$hmmrid.RHMMR1") or die;
		
		print "Making new submit for $hmmrid\n";
		print NEWSub "#!/bin/bash\n";
		print NEWSub "#PBS -N $subid\n";
		print NEWSub "#PBS -W group_list=whiteman\n";
		print NEWSub "#PBS -l jobtype=small_mpi\n";
		print NEWSub "#PBS -l select=1:ncpus=12:mem=24000mb\n";
		print NEWSub "#PBS -l place=pack:shared\n";
		print NEWSub "#PBS -l pvmem=25598mb\n";
		print NEWSub "#PBS -l cput=240:00:0\n";
		print NEWSub "#PBS -l walltime=20:00:0\n";
		print NEWSub "#PBS -q standard\n";
		print NEWSub "#PBS -M dlaroche\@email.arizona.edu\n\n";
		print NEWSub "cd \$PBS_O_WORKDIR # this cd's into the directory you submitted the script from\n";
		print NEWSub "module load hmmer\n";
		print NEWSub "time phmmer --tblout /xdisk/rlapoint/$hmmrid.Rhmmr1 -o /dev/null --noali $query $refspp \n";

		close NEWSub;
		my $cmd = "qsub ./RHMMR2/$hmmrid.RHMMR1";#submit pbs script just made
		system($cmd);#actual submission of command
		print "Submitted $hmmrid.RHMMR1 to the HPC\n";
		$count++;
	}
		

	
print "Submitted $count RHMMR jobs\n";
