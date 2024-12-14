#!/usr/bin/perl
#author:wangxin
### function: From the insertion events to 

use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","p:s", "a:s","b:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}|| !defined $opts{o} || !defined $opts{p} || !defined $opts{a} || !defined $opts{b} ) {
       	die "************************************************************************
       	Usage: $0.pl
       		-i: Insertion events with gene annotation
			-a: Column line for gene annoation
			-b: Column line for distance (We only considered the inserted DNA overlapped with genes) 
			-p: Gene list of Essential gene defined by Giaever
			-o: output List information for all the inserted genes annotated with or without Essential genes
************************************************************************\n";
}




my %hugo;

my $symbol=$opts{p};

open SYM,"$symbol" or die $!;

while (<SYM>){
	chomp;
	my ($geneid)=(split/\s+/,$_)[0];
	#my $ensembl= (split/\./,$ensembl)[0];
	$hugo{$geneid}=$_;
}

close SYM;

my $input=$opts{i};
my $output=$opts{o};
my $n=0;
my $Cannote=$opts{a};
my $Cdis=$opts{b};
open IN,"$input" or die $!;

my %hash; my %str;

while (<IN>){
	chomp;
	s/\r//;
	my ($id,$annotation,$distance)=(split/\t/,$_)[0,$Cannote,$Cdis];
	next unless ($annotation=~/^gene/ && $distance ==0);
	my $Igene=(split/\|/,$annotation)[1];
	$hash{$Igene}++;
	$str{$Igene}.=$id.";";

}
close IN;

open OUT,">$output" or die $!;

foreach my $i (keys %hash){
	my $inf=(exists $hugo{$i})?"Y":"N";
	
	print OUT "$i\t$inf\t$str{$i}\n";

}

close OUT;