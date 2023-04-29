#!/usr/bin/perl
use strict;
use warnings;

my $version="1.0 version";
# Author: Xin Wang                                                   
# Email: xin.wang@childrens.harvard.edu                            
# Copyright (c) 2022 Dr. Kaifu lab                                   
# PI: Kaifu Chen    
                                               
## Description: 

## This is to identify the best blast result for each Ty1 insertion

## Requirement: 1. If they have the mutiple blast reads, we retained the maximum matches result, negeleting the gap size, mismatches
## If they could align to the Chromsome 3, they'd better locate from 294400 - 294500.


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","g:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} || !defined $opts{g} ) {
	die "************************************************************************
	Usage: $0.pl
		-i Blast results
		-g Genotype name
		-o Best Blast Hit of Insertion
************************************************************************\n";
}

my $output=$opts{o};
my $in=$opts{i};
my $genotype=$opts{g};


my %hash; my %inf; my %length;


open I,"$in" or die $!;

while (<I>) {
	chomp;
    # print "$_" ;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	my $information=$_;

	## check the alignments
	my $cov=($match/$qlength);
	my $start0=($pstart<$pend)?$pstart:$pend;
	my $end0=($pstart<$pend)?$pend:$pstart;
	my $length0=$qlength-$start0;
	next if ($end0<340);

	if (exists $inf{$qid} && ($match - $length{$qid})>=6){
		$inf{$qid}=$_;
		next;
	}elsif (exists $inf{$qid} && ($match - $length{$qid}) <6){
		next;
	}

	$length{$qid}=$match;
	$inf{$qid}=$_;


}
close I;

open O,">$output.blast" or die $!;
open BED,">$output.bed" or die $!;
my $fn=0;
foreach my $i (sort keys %inf){
	print O "$inf{$i}\n";
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$inf{$i};
	$fn++;
	my $start1=($pstart<$pend)?$pstart:$pend;
	my $end1=($pstart<$pend)?$pend:$pstart;
	my $strand=($pstart<$pend)?"+":"-";
	print BED "$pid\t$start1\t$end1\t$genotype.$fn\t0\t$strand\t$qid\n";

}

close O;
