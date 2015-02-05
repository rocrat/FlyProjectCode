#! usr/bin/env perl


use strict;
use warnings;

my ($from, $to) = @ARGV;

open(LIST,"finalGeneList.csv");

while (my $id = <LIST>){
	if ($id =~ /(FBgn\d+)/){
		my $ntalign = $1."V.nt_cleanali.fasta";
		my $aaalign = $1."V.aa_cleanali.fasta";
		my $cmmd1 = "cp $from/$ntalign $to";
		my $cmmd2 = "cp $from/$aaalign $to";
		system($cmmd1);
		system ($cmmd2);
	}
}