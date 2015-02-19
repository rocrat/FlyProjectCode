#! usr/bin/env perl


use strict;
use warnings;

my ($from, $to) = @ARGV;

open(LIST,"finalGeneList.csv");

while (my $id = <LIST>){
	if ($id =~ /(FBgn\d+)/){
		my $ntalign = $1."V3.fasta";
		my $cmmd1 = "cp $from/$ntalign $to";
		
		system($cmmd1);
		
	}
}