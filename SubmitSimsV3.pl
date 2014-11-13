#!/usr/bin/env perl
use warnings;
use strict;
#use File::copy#shouldn't be necessary

# This program finds files with p-value <0.05 from the dSSLS .csv file and then converts them from .fasta alignment files
#to phylip alignment files for use with phylobayes 



my $file = $ARGV[0];#open dSSLS file
if (!defined $file) {
  die "Usage: $0 dSSLSfile\n";
}
my $count = 0;
my $fasext = "V.nt_cleanali.fasta";#the file extension of the cleaned alignments
print "Writing files and submitting simulations...\n";
open (INFIL,$file);
open (PBSpred,">ppredcmds.pbs");
while( my $line = <INFIL>){#read in dSSLS file
		my @cols = split /\s+/, $line;
		my $fileid = $cols[0];
		my $pvalue = $cols[1];
	
		if ($pvalue < 0.05){#for each file id with dvalue>0 make pbs script and submit it
			$count++;
			#First convert 
			my $fasid = $fileid.$fasext;
			open (INALI,"/xdisk/rlapoint/CleanAlignmentsV3/$fasid") or die;#Open alignment file from clean alignment folder
			open (OUTALI,">$fileid.phy") or die;#create phylip alignment file
			my @lines;
			while (my $line = <INALI>){#make array of the lines in the alignment file
				push @lines, $line;
			}
			my $nlines = @lines;
			my $nspp = $nlines/2; #get the number of species in the fasta file
			my $length = length($lines[1])-1;#find the length of the sequence minus the newline character
			print OUTALI "$nspp $length\n";#print the header of the alignment file
			foreach my $lin (@lines){
				#if ($lin !~ /\d+/){#skip header
					if ($lin =~ /\>/){#find species names and print first four letters
						my @spcols = split /,/, $lin;
						print OUTALI substr($spcols[1],0,4),"\t";
					}
					if ($lin !~ /\>/){#find sequences and print them after name
						print OUTALI $lin;
					}
				#}
			}
			close INALI;
			close OUTALI;
			#end of file conversion portion.
			#begin making pbs scripts to submit simulations 
			print "Converted $fileid.$fasext to phylip format\n";
			
			
			my $sim = "Sim";
			my $id1 = $fileid;
			
			open(PBSscript,">$fileid.simsubmit");
			print PBSscript "#!/bin/bash\n";
			print PBSscript "#PBS -N $id1\n";
			print PBSscript "#PBS -W group_list=whiteman\n";
			print PBSscript "#PBS -l jobtype=small_mpi\n";
			print PBSscript "#PBS -l select=1:ncpus=2:mem=3600mb\n";
			print PBSscript "#PBS -l place=pack:shared\n";
			print PBSscript "#PBS -l pvmem=25598mb\n";
			print PBSscript "#PBS -l cput=300:00:0\n";
			print PBSscript "#PBS -l walltime=150:00:0\n";
			print PBSscript "#PBS -q standard\n";
			print PBSscript "#PBS -M dlaroche\@email.arizona.edu\n\n";
			print PBSscript "cd \$PBS_O_WORKDIR # this cd's into the directory you submitted the script from\n";
			print PBSscript "time pb -d $fileid.phy -T SppTreeForSims.tre -cat -gtr -x 10 10000 -s /xdisk/rlapoint/SimsV3/$fileid.phy_1 &\n";#?
			print PBSscript "time pb -d $fileid.phy -T SppTreeForSims.tre -cat -gtr -x 10 10000 -s /xdisk/rlapoint/SimsV3/$fileid.phy_2 &\n";#?
			print PBSscript "wait";
			close PBSscript;
			#end making pbs script
			print PBSpred "time ppred -x 1000 10 /xdisk/rlapoint/SimsV3/$fileid.phy_1\n";
			print PBSpred "time ppred -x 1000 10 /xdisk/rlapoint/SimsV3/$fileid.phy_2\n";
			
			
			#end making ppred script
			my $cmd = "qsub $fileid.simsubmit";#submit pbs script just made
			system($cmd);#actual submission of command
			print "Submitted $fileid.simsubmit to the HPC\n";
		}
}
close INFIL;
print "Finished submiting $count simulations\n\n";

close PBSpred;

