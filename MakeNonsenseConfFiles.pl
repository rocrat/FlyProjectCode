#!usr/env/bin perl

#This script is designed to make new garli conf files for testing a non-sensical alternative hypothesis

my @filelist = `ls /xdisk/rlapoint/garliConfNoSnothaV3/ `;

#create 5 pbs scripts for submitting many many garli commands
for (my $i = 1; $i <= 5; $i++){
	my $sub = *SUB.$i;
	my $name = "nonsense".$i;
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
	open(INFIL,"/xdisk/rlapoint/garliConfNoSnothaV3/".$oldfile);
	open(OUTFIL,">/xdisk/rlapoint/garliConfNonsenseV1/".$oldfile);
	while (my $line = <INFIL>){
		if ($line =~ /garliConfFiles/){
			print OUTFIL "streefname = /home/u30/dlaroche/NonsenseSubFiles/DsuzvDgriV1.tre\n";
		}elsif($line =~ /ofprefix/){
			my $fl = substr($oldfile,0,12);
			print OUTFIL "ofprefix = /xdisk/rlapoint/garliNonsenseV1/$fl\n";
		}else{
			print OUTFIL $line;
		}
	}
	if ($count < $numfls/5){#populate commands into 5 qsub files
		print SUB1 "Garli-2.01 /xdisk/rlapoint/garliConfNonsenseV1/$oldfile &\n";
	}
	if ($count > $numfls/5 & $count < 2*$numfls/5){#populate commands into 5 qsub files
		print SUB2 "Garli-2.01 /xdisk/rlapoint/garliConfNonsenseV1/$oldfile &\n";
	}
	if ($count > 2*$numfls/5 & $count < 3*$numfls/5){#populate commands into 5 qsub files
		print SUB3 "Garli-2.01 /xdisk/rlapoint/garliConfNonsenseV1/$oldfile &\n";
	}
	if ($count > 3*$numfls/5 & $count < 4*$numfls/5){#populate commands into 5 qsub files
		print SUB4 "Garli-2.01 /xdisk/rlapoint/garliConfNonsenseV1/$oldfile &\n";
	}
	if ($count > 4*$numfls/5){#populate commands into 5 qsub files
		print SUB5 "Garli-2.01 /xdisk/rlapoint/garliConfNonsenseV1/$oldfile &\n";
	}
	$count++;
	close (INFIL);
	close (OUTFIL);
}
#close each pbs script
for (my $i = 1; $i <= 5; $i++){
	my $sub = *SUB.$i;
	print $sub "wait\n";
	close ($sub);
}
