#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Index::Fasta;
use Bio::DB::Fasta;
use Bio::Index::GenBank;
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";
$ENV{BIOPERL_INDEX} = ".";

sub makehash {#make a hash from an individual file using the original ref spp which should be listed first
	my %hash = ();
	my $file = shift;
	my $col = shift;
	open (INFIL,$file) or die;
		while (my $line = <INFIL>){
			my @cols;
			my $pepid;
			my $geneid;
			if ($line =~ /^[\>]/){
				@cols = split /\s+/, $line;
				foreach $col (@cols){
					if ($col =~ /(gene)\:(FBgn\d+)/){
					$pepid = substr($cols[0],1);
					$geneid = $2;
			
						if (!exists $hash{$pepid}){#if a new ID then assign the score to the hash key
							$hash{$pepid} = $geneid;#the matched id is assigned as the value and the refspp is assigned as the key
						}
					}		
				}	
			}
		}
	
	close INFIL;
	print "$file\n";
	return %hash;
}
sub makehashc {#make a hash from a cds file to match a peptide id or gene id to a coding sequence
	my %hash = ();
	my $file = shift;
	my $col = shift;
	open (INFIL,$file) or die;
	while (my $line = <INFIL>){
		if ($line =~ /^[\>]/){
			
			my @cols = split /\>/, $line;
			my $id;
			#first use peptide id if available 
			if ($cols[1] =~ /FlyBase:(FBpp\d+)/){#match peptide id
				$id = $1;
				
			}
			#if there is no peptide id then try to match gene id
			if (!defined $id){	
				if ($cols[1] =~ /(FBgn\d+)/){#match gene id
					$id = $1;
					
				}
			}
			if (defined $id){
				if (!exists $hash{$id}){
					$hash{$id} = $cols[2];
				}
			}
		}
	}
	close INFIL;
	print "$file\n";
	return %hash;
}
my %dmelH = &makehash("DmelanogasterPEP.fasta");
my %dmelHC = &makehashc("DmelforHashCDSStopsRemoved.txt");
my %dbiaH = &makehash("DbiarmipesPEP.fasta");
my %dbiaHC = &makehashc("DbiaforHashCDSStopsRemoved.txt");
my %dgriH = &makehash("DgrimshawiPEP.fasta");
my %dgriHC = &makehashc("DgriforHashCDSStopsRemoved.txt");
my %dmojH = &makehash("DmojavensisPEP.fasta");
my %dmojHC = &makehashc("DmojforHashCDSStopsRemoved.txt");
my %dpseH = &makehash("DpseudoobscuraPEP.fasta");
my %dpseHC = &makehashc("DpseforHashCDSStopsRemoved.txt");
my %dsuzH = &makehash("DsuzukiiPEP.fasta");
my %dsuzHC;#
my %dyakH = &makehash("DyakubaPEP.fasta");
my %dyakHC = &makehashc("DyakforHashCDSStopsRemoved.txt");
my %sflaHC;
my %snotHC;
open (SFLA,"SflavaforHashCDSStopsRemoved.txt") or die;
while (my $line = <SFLA>){
	if ($line =~ /^[\>]/){
		my @cols = split /,/, $line;
		my $id = substr($cols[0],1);
		if (!exists $sflaHC{$id}){
			$sflaHC{$id} = $cols[1];
		}
	}
}
open (SNOT,"SnothaforHashCDSStopsRemoved.txt") or die;
while (my $line = <SNOT>){
	if ($line =~ /^[\>]/){
		my @cols = split /,/, $line;
		my @idcol = split /\s+/, $cols[0];
		my $id = substr($idcol[0],1);
		if (!exists $snotHC{$id}){
			$snotHC{$id} = $cols[1];
		}
	}
}
open (DSUZ,"DsuzforHashNEWCDSStopsRemoved.txt") or die;
while (my $line = <DSUZ>){
	if ($line =~ /^[\>]/){
		my @cols = split />/, $line;
		my @idcol = split /\s+/, $cols[1];
		my $id = substr($idcol[0],0,13);
		
		if (!exists $dsuzHC{$id}){
			$dsuzHC{$id} = $cols[2];
		}
	}
}



my $dbi = Bio::DB::Fasta->new("DbiarmipesCDSStopsRemoved.fasta");


my $loc_factory = Bio::Factory::FTLocationFactory->new;

sub getCodingSeq {#subroutine to get the coding sequence from the protein gi
	my $pdb = shift;
	my $id = shift;
	my $ndb = shift;
	my $prot_obj = $pdb->fetch($id);
	my $nuc_obj;
	my $nucid;
	if (defined $prot_obj){
		foreach my $feat ( $prot_obj->top_SeqFeatures ) {
		if ( $feat->primary_tag eq 'CDS' ) {
			# example: 'coded_by="U05729.1:1..122"'
			my @coded_by = $feat->each_tag_value('coded_by');
			my ($nuc_acc,$loc_str) = split /\:/, $coded_by[0];
			$nucid = $nuc_acc;
			#print $nuc_acc,"\n";
			$nuc_obj = $ndb->fetch($nuc_acc);
			
		}
	}
	}
	return ($nucid,$nuc_obj);
}

open (INFIL,"CombinedRecipClust.list") or die "Cannot open the match file";
open (CNTFIL,">FastaFilesKeyV4NoSnot.txt");
open (CNTFILSPP,">FastaFilesKeySppV4NoSnot.txt");
open (SPPNUMKEY,">SppNumberKey.csv");
open (FASTKEY,">FastaUniqueSPPV4");
my $n = 1;#counter for Sfla sequence id since the long name won't play nice with Gblocks

while (my $line = <INFIL>){
	my %spphsh;
	my %fastahsh;
	my $filename;
	if($line=~/Dmel\.(FBgn\d+),FBpp\d+/){  
		$filename = $1."V4";#names the outgoing fasta file with Dmel anchor geneID
	}
	#first need to parse out the row name
	#$line =~ s/\d+,//;
	my @cols = split /;/, $line;
	my $numspp = 0; #counter for the number of species in the fasta file (to be printed in the key)
	my $spplist = "";#initiate string to be printed in species key
	open (NFIL,">$filename.fasta");

	#use species ids at start of id to identify appropriate columns and keep track of the number of species included on a line
	foreach my $geneid (@cols){#loop through each entry, identify the gene id and print the sequence in the fasta file
	if($geneid !~ /^\d+/){
		
		my $dmelgene;
		my $dmelpep;# peptide id;
		if($geneid =~ /(Dmel)\./){
			if($geneid =~ /(FBgn\d+),(FBpp\d+)/){
				$dmelgene = $1;
				$dmelpep = $2;
			}
			if(!exists $spphsh{$1}){
				$spphsh{$1} = 1;
			}else{
				$spphsh{$1}++;
			}
			#get sequences
			my $dmelseq;
			if(exists $dmelHC{$dmelpep}){
				$dmelseq = $dmelHC{$dmelpep}; 
			}else{
				$dmelseq = $dmelHC{$dmelgene};
			}
		if (defined $dmelseq){
			print NFIL ">$dmelgene,Dmel\n";
			print NFIL $dmelseq,;
			$numspp++;
			$spplist = "Dmel";
			if(!exists $fastahsh{"Dmel"}){
				$fastahsh{"Dmel"} = 1;
			}else{
				$fastahsh{"Dmel"} = 1;
			}
		}
		}	
		my $dbiid;
		if($geneid =~ /(Dbia)\./){
			$dbiid = substr($geneid,5,-3)."-RA";
			if(!exists $spphsh{$1}){
				$spphsh{$1} = 1;
			}else{
				$spphsh{$1}++;
			}
		my $dbiseq = $dbi->get_Seq_by_id($dbiid);# Dbia is a little difference than the other spp's
		if (defined $dbiseq){
			print NFIL ">$dbiid,Dbia\n";
			print NFIL $dbiseq->seq,"\n";
			$numspp++;
			$spplist = $spplist.",Dbia";
			if(!exists $fastahsh{"Dbia"}){
				$fastahsh{"Dbia"} = 1;
			}else{
				$fastahsh{"Dbia"} = 1;
			}
		}
		}
		my $dgrigene;
		my $dgriid;
		if($geneid =~ /(Dgri)\./){
			$dgriid = substr($geneid,5);
			print "$dgriid \n";
			$dgrigene = $dgriH{$dgriid};#get the gene id from the peptide id in the list file
			if(!exists $spphsh{$1}){
				$spphsh{$1} = 1;
			}else{
				$spphsh{$1}++;
			}
			my $dgriseq; 
			#try to get the sequences
			if(exists $dgriHC{$dgriid}){
				$dgriseq = $dgriHC{$dgriid}; 
			}else{
				$dgriseq = $dgriHC{$dgrigene};
			}
			if (defined $dgriseq){
				print NFIL ">$dgriid,Dgri\n";
				print NFIL $dgriseq;
				$numspp++;
				$spplist = $spplist.",Dgri";
				if(!exists $fastahsh{"Dgri"}){
					$fastahsh{"Dgri"} = 1;
				}else{
					$fastahsh{"Dgri"} = 1;
				}
			}
		}
		my $dmogene;
		my $dmoid;
		if($geneid =~ /(Dmoj)\./){
			$dmoid = substr($geneid,5);
			$dmogene = $dmojH{$dmoid};#get the gene id from the peptide id in the list file
			if(!exists $spphsh{$1}){
				$spphsh{$1} = 1;
			}else{
				$spphsh{$1}++;
			}
			my $dmoseq; 
			if(exists $dmojHC{$dmoid}){
				$dmoseq = $dmojHC{$dmoid}; 
			}else{
				$dmoseq = $dmojHC{$dmogene};
			}
		if (defined $dmoseq){
			print NFIL ">$dmoid,Dmoj\n";
			print NFIL $dmoseq;
			$numspp++;
			$spplist = $spplist.",Dmoj";
			if(!exists $fastahsh{"Dmoj"}){
				$fastahsh{"Dmoj"} = 1;
			}else{
				$fastahsh{"Dmoj"} = 1;
			}
		}
		}
		my $dpsegene;
		my $dpseid;
		if($geneid =~ /(Dpse)\./){
			$dpseid = substr($geneid,5);
			$dpsegene = $dpseH{$dpseid};#get the gene id from the peptide id in the list file
			if(!exists $spphsh{$1}){
				$spphsh{$1} = 1;
			}else{
				$spphsh{$1}++;
			}
			my $dpseseq;
			if(exists $dpseHC{$dpseid}){
				$dpseseq = $dpseHC{$dpseid}; 
			}else{
				$dpseseq = $dpseHC{$dpsegene};
			}
		if (defined $dpseseq){
			print NFIL ">$dpseid,Dpse\n";
			print NFIL $dpseseq;
			$numspp++;
			$spplist = $spplist.",Dpse";
			if(!exists $fastahsh{"Dpse"}){
				$fastahsh{"Dpse"} = 1;
			}else{
				$fastahsh{"Dpse"} = 1;
			}
		}
		}
		my $dsuzid;
		if($geneid =~ /(Dsuz)\./){
			$dsuzid = substr($geneid,5);
			if(!exists $spphsh{$1}){
				$spphsh{$1} = 1;
			}else{
				$spphsh{$1}++;
			}
			my $dsuseq = $dsuzHC{$dsuzid}; 
			if (defined $dsuseq){
				print NFIL ">$dsuzid,Dsuz\n";
				print NFIL $dsuseq;
				$numspp++;
				$spplist = $spplist.",Dsuz";
				if(!exists $fastahsh{"Dsuz"}){
					$fastahsh{"Dsuz"} = 1;
				}else{
					$fastahsh{"Dsuz"} = 1;
				}
			}
		}
		my $dyakgene; 
		my $dyakid;
		if($geneid =~ /(Dyak)\./){
			$dyakid = substr($geneid,5);
			$dyakgene = $dyakH{$dyakid};#get the gene id from the peptide id in the list file
			if(!exists $spphsh{$1}){
				$spphsh{$1} = 1;
			}else{
				$spphsh{$1}++;
			}
			my $dyakseq; 
			if(exists $dyakHC{$dyakid}){
				$dyakseq = $dyakHC{$dyakid}; 
			}else{
				$dyakseq = $dyakHC{$dyakgene};
			}
			if (defined $dyakseq){
			print NFIL ">$dyakid,Dyak\n";
			print NFIL $dyakseq;
			$numspp++;
			$spplist = $spplist.",Dyak";
				if(!exists $fastahsh{"Dyak"}){
					$fastahsh{"Dyak"} = 1;
				}else{
					$fastahsh{"Dyak"} = 1;
				}
			}
		}
		my $sflaid;
		if($geneid =~ /(Sfla)\./){
			$sflaid = substr($geneid,5);
			if(!exists $spphsh{$1}){
				$spphsh{$1} = 1;
			}else{
				$spphsh{$1}++;
			}
			my $sflaseq = $sflaHC{$sflaid};
			if (defined $sflaseq){
				print NFIL ">$n,Sfla\n";
				print NFIL $sflaseq;
				$numspp++;
				$spplist = $spplist.",Sfla";
				if(!exists $fastahsh{"Sfla"}){
					$fastahsh{"Sfla"} = 1;
				}else{
					$fastahsh{"Sfla"} = 1;
				}
			}
		}
		$n = $n+1;
	}
	}#close loop through each column
	
	print CNTFIL "$filename,$numspp\n";
	print CNTFILSPP "$filename,$spplist\n";
	#print spp ids and number of genes from each spp
	my @sppids = keys %spphsh;
	my @sorted = sort@sppids;
	print SPPNUMKEY scalar @sppids,",";#print the numbher of unique spp
	close NFIL;
	foreach my $id (@sorted){
		print SPPNUMKEY "$id,$spphsh{$id},";
	}
	print SPPNUMKEY "\n";
	my @uniquespp = keys %fastahsh;
	print FASTKEY scalar @uniquespp,";",join(";",@uniquespp),"\n";
}
close CNTFIL;
close CNTFILSPP;
close INFIL;