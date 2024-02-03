#!/usr/bin/perl
use strict;
use warnings;


#### To extract the expression of request gene list


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","c:s","o:s","m:s","n:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{c}   || !defined $opts{o} || !defined $opts{m}   || !defined $opts{n}  ) {
        die "************************************************************************
        Usage: $0.pl	
						-i request Gene list
						-c index of Expression of each id
						-m column number of the index ID
						-n column number of the request ID
						-o output of overlapping GO enrichment
************************************************************************\n";
}


my $exp=$opts{c};
my $out=$opts{o};

my $req=$opts{i};

my $mnum=($opts{m})?$opts{m}:0;
my $nnum=($opts{n})?$opts{n}:0;

my %hash;
my %str;


open C,"$exp" or die $!;

while (<C>){
        chomp;
        s/\r//;
        my $geneA=(split/\t/,$_)[$mnum];
		## store multiple line
        $str{$geneA}.= $_."\n";

}

close C;


open I,"$req" or die $!;
open O,">$out" or die $!;
while (<I>){
        chomp;
        s/\r//;
        my $name=(split/\t/,$_)[$nnum];	
	next if (exists $hash{$name});
	$hash{$name}++;
        print O "$str{$name}";
}

close I;
close O;
