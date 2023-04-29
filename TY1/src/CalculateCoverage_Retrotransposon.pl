#!/usr/bin/perl
use strict;
use warnings;


#### 
### This script is to calculate each nucleotide coverage


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","b:s","a:s","o:s","g:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{a}|| !defined $opts{b}  || !defined $opts{o} || !defined $opts{g}  ) {
	die "************************************************************************
	Usage: $0.pl	 -i Insertion Retro bed file
				-a Retrotransposon annotation bed file
				-b Genome Retrotransposon Sequence
				-g Genotype
				-o Position with nucleotide and annotation information
************************************************************************\n";
}



my $genome=$opts{b};
my $genotype=$opts{g};

## index of genome postion
my %str; my $i;
open GE, "$genome" or die $!;
while (<GE>){
	chomp;
	if(/>(\S+)/){
		$i=$1;
	}else{
		$str{$i}.=$_;
	}
}
close GE;

#### rDNA annotation information  ####

my $rdna=$opts{a};
my %annot;
open RD,"$rdna" or die $!;
while (<RD>){
	chomp;
	s/\r//g;
	my ($chr,$st1,$end1,$an2)=(split/\t/,$_)[0,1,2,6];
	foreach my $q ($st1..$end1){
		$annot{$q}.=$an2.";";
	}
}
close RD;




### read the insertion sequence
my $in=$opts{i}; 
my %num; 
open IN,"$in" or die $!;
while (<IN>){
	chomp;
	my ($start,$end,$an)=(split/\t/,$_)[1,2,12];
	
	foreach my $n ($start..$end){
		$num{$n}++;
		#$annot{$n}=$an;
		#$nuc{$i}=substr $str{$chr},$postition,1;
	}	
}
close ;

my $output=$opts{o};
open OUT,">$output" or die $!;
my $chr="YPLWTy1-1"; my $n=0;

print OUT "Chromosome\tPosition\tNuclearSequence\tRankNumber\t$genotype.InsertionCounts\tRetrotransposonAnnotation\n";
foreach my $m (1..5924){
	my $cov=(exists $num{$m})?$num{$m}:0;
	my $postion=$m-1;
	my $nuc=substr $str{$chr},$postion,1;
	my $annotation=(exists $annot{$m})?$annot{$m}:0;
	$n++;
	print OUT "YPLWTy1-1\t$m\t$nuc\t$n\t$cov\t$annotation\n";
	
}

close OUT;


