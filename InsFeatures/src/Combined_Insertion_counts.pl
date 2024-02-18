
#!/usr/bin/perl
#author:wangxin
### function: generating a segements.txt, stucture.lable.txt, a lists of measure files for circos across different cancer types
### 
use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i} || !defined $opts{o}) {
       	die "************************************************************************
       	Usage: $0.pl
       		-i: query folder (mutator folder)
			-o: output name
************************************************************************\n";
}


### read the mutators files in the folder
my $mutators=$opts{i};
my $output=$opts{o};
my @file;

opendir MT,"$mutators" or die $!;
 @file=readdir MT;
closedir MT;
my $n=0; my %hash; my %str;
foreach my $i (@file){
	next unless ($i=~/.LTR.cov.txt$/);
	
	print "$i\n";
	my $genotype=$i;
	$genotype =~ s/\.LTR\.cov\.txt//;
	$str{$genotype}++;
	$n++;
	#### read clinical file in each folder: and then combine the mutaiton, clinical restults. 	
	open FILE,"$mutators/$i" or die $!;
	
	while (<FILE>){
		chomp;
		s/\r//;
		s/^\s+//;
		my ($num,$cov)=(split/\t/,$_)[1,4];
		next if ($num eq "Position");
		$hash{$num}->{$genotype}=$cov;
	}
	

	close FILE;

}

#my @array=("WT-01","WT-02","Q128HTT16-01","Q128HTT16-02");


open OUT,">$output" or die $!;
print OUT "Position";
foreach my $q (sort keys %str){
	print OUT "\t$q";
}
print OUT "\n";

foreach my $i (sort {$a<=>$b} keys %hash){
	
	print OUT "$i";
	foreach my $g (sort keys %str){
		print OUT "\t$hash{$i}->{$g}";
	}
	print OUT "\n";
}

close OUT;
