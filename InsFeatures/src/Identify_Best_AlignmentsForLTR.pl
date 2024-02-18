#!/usr/bin/perl
use strict;
use warnings;
### This is to identify the best blast result for each reads:
### requirement: 1. If they have the mutiple blast reads, we retained the maximum matches result, negeleting the gap size, mismatches 
#### If they could align to the Chromsome 3, they'd better locate from 294400 - 294500.
#use Bio::Seq;
#use Bio::SeqIO;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o}  ) {
	die "************************************************************************
	Usage: $0.pl 
		-i Blast result
		-o Best Blast Hit of Insertion
************************************************************************\n";
}

my $output=$opts{o};
my $in=$opts{i};



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

open O,">$output" or die $!;

foreach my $i (sort keys %inf){
	print O "$inf{$i}\n";
}

close O;


