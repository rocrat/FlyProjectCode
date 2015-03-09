#!usr/env/bin perl

#This script is designed to make new garli conf files for testing a non-sensical alternative hypothesis

my @filelist = `ls /xdisk/rlapoint/garliConfNoSnothaV3/ `;
my %tophsh;
open(TOP,"FinalnonsenseGeneList.csv");
while (my $gene = <TOP>){
	chomp($gene);
	if($gene =~ /FBgn/){
		$tophsh{$gene} = $gene;
		if(exists $tophsh{$gene}){
		#print "$tophsh{$gene}\n";
		}
	}
}
my @glist = keys %tophsh;
foreach my $gen (@glist){
	#print "$gen\n";
}
#create 5 pbs scripts for submitting many many garli commands
for (my $i = 1; $i <= 1; $i++){#just changed to 1 iteration
	my $sub = *SUB.$i;
	my $name = "nonsenseTree".$i;
	open($sub,">SubmitNonsenseGarli".$i.".pbs");#create file
	print $sub "#!/bin/bash\n";#create the header in each file
	print $sub "#PBS -N $name\n";
	print $sub "#PBS -W group_list=whiteman\n";
	print $sub "#PBS -l jobtype=small_mpi\n";
	print $sub "#PBS -l select=1:ncpus=12:mem=24000mb\n";
	print $sub "#PBS -l place=pack:shared\n";
	print $sub "#PBS -l pvmem=25598mb\n";
	print $sub "#PBS -l cput=600:00:0\n";
	print $sub "#PBS -l walltime=50:00:0\n";
	print $sub "#PBS -q standard\n";
	print $sub "#PBS -M dlaroche\@email.arizona.edu\n\n";
	print $sub "cd \$PBS_O_WORKDIR # this cd's into the directory you submitted the script from\n";
}
#create conf files and populate 5 submission scripts
my $count = 0;
my $numfls = $#filelist;
print "$numfls\n";
foreach my $oldfile (@filelist){
	chomp($oldfile);
	my $geneid;
	if($oldfile =~ /(FBgn\d{7})/){
		$geneid = $1;
		#print "$geneid\n";
	}
	#print "$tophsh{$geneid}\n";
	if(exists $tophsh{$geneid}){
	print "$geneid found \n";
		open(INFIL,"/xdisk/rlapoint/garliConfNoSnothaV3/".$oldfile);
		open(OUTFIL,">/xdisk/rlapoint/simTreeConfNonsenseV1/".$oldfile);
		while (my $line = <INFIL>){
			if ($line =~ /garliConfFiles/){
				print OUTFIL "streefname = /home/u30/dlaroche/NonsenseSubFiles/simTreesNonsenseV1.nexus\n";
			}elsif($line =~ /ofprefix/){
				my $fl = substr($oldfile,0,12);
				print OUTFIL "ofprefix = /xdisk/rlapoint/simTreeNonsenseV1/$fl.simTree\n";
			}else{
				print OUTFIL $line;
			}
		}
		print SUB1 "Garli-2.01 /xdisk/rlapoint/simTreeConfNonsenseV1/$oldfile &\n";
	
		$count++;
		close (INFIL);
		close (OUTFIL);
	}
}
#close each pbs script
for (my $i = 1; $i <= 1; $i++){
	my $sub = *SUB.$i;
	print $sub "wait\n";
	close ($sub);
}
