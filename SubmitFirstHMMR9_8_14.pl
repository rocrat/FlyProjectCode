#!/usr/bin/env perl
use warnings;
use strict;

#find CDS fasta files from a list of files and submit an initial HMMR using one spp as a hub
#use perl SubmitFirstHMMR.pl to run from the directory where the fasta files are located
#Then enter the pep file used as a reference as the first argument to the call

#chomp(my $refspp = <STDIN>);#get the reference spp
my @allfiles = `ls /xdisk/rlapoint | grep 'PEP\.fasta'`;
my $reffile = $ARGV[0];
my $refspp = "Dmel";

my $count=0;

foreach my $file (@allfiles){
	chomp $file;
	print "$file\n";
	if($file !~ /$refspp/){
		my $sppid = substr($file,0,-9);
		my $hmmrid = $refspp."v".$sppid;
		open (NEWSub,">/xdisk/rlapoint/HMMR1/$hmmrid.HMMR1") or die;
		
		print "Making new submit for $hmmrid\n";
		print NEWSub "#!/bin/bash\n";
		print NEWSub "#PBS -N $sppid\n";
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
		print NEWSub "time phmmer --tblout /xdisk/rlapoint/$hmmrid.hmmr1 -o /dev/null --noali $reffile $file \n";

		close NEWSub;
		my $cmd = "qsub ./HMMR1/$hmmrid.HMMR1";#submit pbs script just made
		system($cmd);#actual submission of command
		print "Submitted $hmmrid.HMMR1 to the HPC\n";
		$count++;
	}
		
}
	
print "Submitted $count HMMR jobs\n";
