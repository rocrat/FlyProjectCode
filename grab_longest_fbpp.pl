#!/usr/bin/perl;
use strict;
use warnings;
use Bio::SeqIO;

######################################################################################
#
# author: Andy Gloss
# date: 3/5/14
# purpose: takes FlyBase fbpp file (e.g. dmoj-all-translation-r1.3.fasta),and outputs 
#	    a new fasta file containing on the longest fbpp for
#          each parent fbgn, with seqs named by fbgn# only
#
######################################################################################


my $infile = $ARGV[0]; #fasta file of total database
my $outfile = $ARGV[1];

open (OUTFIL, ">$outfile") or die "Cannot open $outfile\n";

#open the total database using seqIO
my $inseqIO = Bio::SeqIO->new(-file => "$infile", -format => 'fasta');

my %fbgn_maxlen;
my %fbgn_maxseq;
my %fbgn_maxfbpp;
my $num_read;

while ( my $seqIO = $inseqIO->next_seq() ) {

	my $name = ();
	my $fbgn_name = ();
	my $fbpp_name = ();
	my $length = ();
	my $seq = ();
	
	$name = $seqIO->display_id() . "\t" . $seqIO->desc();

	if ( $name =~ /(FBpp\d+).+parent\=(FBgn\d+)\,/ ) {

		$fbpp_name = $1;

		$fbgn_name = $2;

	} else {

		print "no pattern match: $name\n\n";

	}

	$length = $seqIO->length;

	$seq = $seqIO->seq();

#	if ( $fbgn_maxlen{$fbgn_name} ) {
#		print "multiple entries for: $fbgn_name \n";
#	}

	if ( (! $fbgn_maxlen{$fbgn_name}) || ($length > $fbgn_maxlen{$fbgn_name}) ) {

		$fbgn_maxlen{$fbgn_name} = $length;

		$fbgn_maxseq{$fbgn_name} = $seq;

		$fbgn_maxfbpp{$fbgn_name} = $fbpp_name;

	}

	$num_read++;

}

my @kept = keys %fbgn_maxfbpp;
my $num_kept = @kept;

print "num read: $num_read\n";
print "num unique: $num_kept\n";

for my $key (keys %fbgn_maxseq) {

	print OUTFIL ">$key\n$fbgn_maxseq{$key}\n";

}





