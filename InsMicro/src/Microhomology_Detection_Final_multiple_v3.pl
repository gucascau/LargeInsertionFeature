#!/usr/bin/perl
use strict;
use warnings;
#use Data::Dump qw(dump);
#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

#use Bio::Seq;
#use Bio::SeqIO;


#### the script is to compare all the insertion/deletion events across samples

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","r:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o}  || !defined $opts{r} || defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -i the folder that contain the indels -r total reads each sample  -o output of proportion of each indel events 
	
	Request Parameters:
	-i the folder that contain the indels
	-r total reads each sample
	-o output of proportion of each indel events 
	
	-h Help
************************************************************************\n";
}

### Read the total high quality reads counts of each sample

my %totalcounts; 

my $RCount=$opts{r};
open RCOU,"$RCount" or die $!;

while (<RCOU>) {
	chomp;
    my($sample,$TRcount)=split/\s+/,$_;
	$totalcounts{$sample}=$TRcount;
	
}
close RCOU;

### Read the folder of indels. we put the files with name_indel.deletion.stat or name_indel.insertion.stat
my $indels=$opts{i};
my $output=$opts{o};
my @file;

opendir MT,"$indels" or die $!;
 @file=readdir MT;
closedir MT;
my $n=0; my %hash; my %num;
foreach my $i (@file){
	next unless ($i=~/.stat/);
	
	print "$i\n";

	$n++;
	my $type=join "_",(split/\_/,$i)[0,1];

	#### read clinical file in each folder: and then combine the mutaiton, clinical restults. 	
	open FILE,"$indels/$i" or die $!;	
	while (<FILE>){
		chomp;
		s/\r//;
		s/^\s+//;
		my ($event,$cov)=(split/\t/,$_)[0,1];
		#next if ($gene eq "ID");
		$hash{$event}->{$type}=$cov/$totalcounts{$type};
		$num{$event}->{$type}=$cov;
	}
	close FILE;

}


open OUT, ">$output" or die $!;
print OUT "Events\t";

foreach my $j (sort keys %totalcounts){
	print OUT "$j\t";
}
print OUT "\n";

foreach my $i (sort keys %hash){
	print OUT "$i\t";
	foreach my $n (sort keys %totalcounts){
		my $final=(exists $hash{$i}->{$n})?$hash{$i}->{$n}:0;
		print OUT "$final\t";
	}
	print OUT "\n";
	
}
close OUT;

