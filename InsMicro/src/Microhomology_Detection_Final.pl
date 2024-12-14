#!/usr/bin/perl
#author:wangxin
# Date: 10-2-2020
### #### This script is to detect the microhomology and generate similar result as yang
#
use strict;
use warnings;
#use Algorithm::NeedlemanWunsch;
#use Statistics::Test::WilcoxonRankSum;
#use List::Util qw(sum);
#use List::Util 'shuffle';

my $version="1.0 version";
use Getopt::Long;
use List::Util qw( min max );
my %opts;
GetOptions(\%opts,"i:s","g:s","b:s","t:s","s:s","o:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}|| !defined $opts{g}|| !defined $opts{b} || !defined $opts{t}  ||!defined $opts{s}||!defined $opts{o}  ) {
       	die "************************************************************************
       	Usage: $0.pl
			-i: Insertion fasta sequence
			-g: Yeast genome sequence
			-b: Donor Blast results
			-t: Two side of MAT insert site
			-s: The insertion result with annotation
			-o: Output of insertion with microhomology
************************************************************************\n";
}

### insertion index
my $Ifasta=$opts{i};

## genome index
my $fasta=$opts{g};

# two side of MAT site index

my $sideblast=$opts{t};

# Insertion donor blast index
my $insertion=$opts{b};

# Insertion with annotation information

my $Insert=$opts{s};

# Output with both annotation and microhomology information

my $output=$opts{o};


#### Here is the genome sequence of each genome

open FASTA, "$fasta" or die "cannot open file $fasta";
my %str; my $id1;
while (<FASTA>){
	chomp;	
	if ($_=~/>(\S+)/){
		$id1=$1;
	}else{
		$str{$id1}.=$_;
	}    
}
close FASTA;


#### index of insertion sequence
open IFASTA, "$Ifasta" or die "cannot open file $Ifasta";
my %str2; my $id2;
while (<IFASTA>){
	chomp;	
	if ($_=~/>(\S+)/){
		$id2=$1;
	}else{
		$str2{$id2}.=$_;
	}    
}
close IFASTA;


### index of insertion information

open INSERT, "$Insert" or die "cannot open file $Insert";
my %insertion;

while (<INSERT>){	
	chomp; 
	my $string=(split/\t/,$_)[1];
	#my $string=join ":",$array[0],@array[5..16];
	$insertion{$string}=$_;
    
}
close INSERT;


#### Here is the index of the two side blast results
my %LeftB; my %RightB;

open TBL, "$sideblast" or die "cannot open file $sideblast";
while (<TBL>){
	chomp;
	my ($id,$start,$end,$length0,$refS,$refE)=(split/\t/,$_)[0,6,7,8,9,10];
	
	next unless($start <10 || ($length0-$end)<10);
	#next if (exists $LeftB{$id});
	#next if (exists $RightB{$id});
	
	if($start<10 && $refS>294340 && $refS<294360){
		next if (exists $LeftB{$id});
		$LeftB{$id}=$end;
	}elsif(($length0-$end)<10 && $refE >294420 && $refE < 294435){
		next if (exists $RightB{$id});
		$RightB{$id}=$start;
	}	
}

close TBL;


#### Here is the index of insertion blast results
my %hash;
my %ins;

open INS,"$insertion" or die "cannot open file $insertion";
while (<INS>){
	chomp;
	my ($id,$chr,$start,$end,$length,$refS,$refE)=(split/\t/,$_)[0,1,6,7,8,9,10];
	next if (exists $hash{$id});
	#next if($start <10 || ($length-$end)<10);
	
	## Here we choose the best blast to represent the insertion events
	next unless($start <60 && $start > 20 && ($length-$end)> 20 && ($length-$end) <70);
	
	
	$hash{$id}++;
	$ins{$id}->{Chr}=$chr;
	$ins{$id}->{Start}=$start;
	$ins{$id}->{End}=$end;
	
	my $type0=($refE>$refS)?"+":"-";
	my $Rstart=($refE>$refS)?$refS:$refE;
	my $Rend=($refE>$refS)?$refE:$refS;
	
	$ins{$id}->{Strand}=$type0;
	$ins{$id}->{Rstart}=$Rstart;
	$ins{$id}->{Rend}=$Rend;
	$ins{$id}->{Length}=$length;
	
}

close INS;


#### index of the MAT-R and MAT-L
my $stringS="GCATAGTCGGGTTTTTCTTTTAGTTTCAGCTTTCCGCAACA";
my $stringE="AACAGTATAATTTTATAAACCCTGGTTTTGGTTTTGTAGAGTGGTTG";


### scoring function:
sub score_sub {
  if (!@_) {
    return -20;            # gap penalty
  }
  ## mismatch scores -1, match +1
  return ($_[0] eq $_[1]) ? 1 : -1;
}




open OUT,">$output" or die $!;

my $Titlename="RepresentativeRead";

print OUT "$insertion{$Titlename}\tMicro_5_length\tMicro_5_seq\tMicro_5_confidence\tChange_5_length\tChange_5_seq\tJunction_5_seq\tDonor_5_seq\tMicro_3_length\tMicro_3_seq\tMicro_5_confidence\tChange_3_length\tChange_3_seq\tJunction_3_seq\tDonor_3_seq\tInsertion_length\tInserted_seq\n";

#print OUT "ID\tMicro_5_length\tMicro_5_seq\tChange_5_length\tChange_5_seq\tJunction_5_seq\tDonor_5_seq\tMicro_3_length\tMicro_3_seq\tChange_3_length\tChange_3_seq\tJunction_3_seq\tDonor_3_seq\tInsertion_length\tInserted_seq\n";

foreach my $i (keys %hash){
	
	#### generate the basic information of string
	
	my $start1=$LeftB{$i};
	my $end1=$RightB{$i};
	
	my $start2=$ins{$i}->{Start};
	my $end2=$ins{$i}->{End};
	my $length=$ins{$i}->{Length};
	
	my $chr=$ins{$i}->{Chr};
	my $start3=$ins{$i}->{Rstart};
	my $end3=$ins{$i}->{Rend};
	my $type=$ins{$i}->{Strand};
	
	my $DonorSeq=$str{$chr};
	
	my $RawSeq=$str2{$i};
	
	my $lengthUp=15;
	###### Here is to detect the microhomology in the 5' junction ######
	
	
	my $JunctionF=0; my $LenmicroF=0; my $SeqLnmicroF=0; my $DelmicroF=0; my $SeqDelmicroF=0; my $SeqDonorF=0; 
	my $Mismatch_statusF="TRUE";
	
	### generate the final start point of the insertion
	my $finalStart=$start1;
	if ($start1<=41){
		
	#### Here is generate the 5' end Junction
		my $JunctionFA=	substr $RawSeq,$start1-15,15;
		my $JunctionFB= substr $RawSeq,$start1,15;	
		$JunctionF= join '',$JunctionFA,$JunctionFB;
		
		
	#### the length of 5' deletion/insertion and the sequence;	 
		
		if($start1>=($start2-1)){
			$DelmicroF= 41-$start1;
			my $SeqDelmicroFR = substr $stringS,(41-$DelmicroF),$DelmicroF;
			$SeqDelmicroF ="-".$SeqDelmicroFR;
		}else{
			my $DelmicroF1=41-$start1;
			my $InmicroF2=$start2-$start1-1;
			$DelmicroF=$DelmicroF1."\,\+$InmicroF2";
			
			my $SeqDelmicroFR1 = substr $stringS,(41-$DelmicroF1),$DelmicroF1;
			my $SeqDelmicroFR2 = substr $RawSeq,$start1,$InmicroF2;
			
			$SeqDelmicroF="\-$SeqDelmicroFR1"."\,\+$SeqDelmicroFR2";
		}
	
		
		
	#### the length of 5' microhomology and the sequence;	### here is for the potential microhomology if there are no mismatches at the junction:
		my $LenmicroFP= ($start1>=$start2)? ($start1-$start2+1):0;		
		my $SeqLnmicroFP = substr $stringS, ($start1-$LenmicroFP), $LenmicroFP;
		
		
		#### further confirmation
		### some donor might have one or two mismatches at the microhomology, we would like consider them as shorter microhomology
		
		
		#### the donor microhomology sequence 
		my $SeqLnmicroFV1 = substr $DonorSeq, ($start3-1), $LenmicroFP;
	
		my $SeqLnmicroFV2 = substr $DonorSeq, ($end3-$LenmicroFP), $LenmicroFP;
		
		my $SeqLnmicroFV2r = reverse $SeqLnmicroFV2;
		
		$SeqLnmicroFV2r =~ tr/ATGCatgc/TACGtacg/;
		
		my $SeqLnmicroFV= ($type eq "+")?$SeqLnmicroFV1:$SeqLnmicroFV2r;
		
	
		##### check whether donor microhomology is consistent with mat micrhomoloyg
		### print the $LenmicroF and $SeqLnmicroF as the final microhomology
		
		if ($SeqLnmicroFP eq $SeqLnmicroFV || (!$SeqLnmicroFP && !$SeqLnmicroFV)){
			
			$LenmicroF=$LenmicroFP;
			$SeqLnmicroF=$SeqLnmicroFP;
			
			$Mismatch_statusF="Confidence";
			#print "$RawSeq\t$SeqLnmicroFP\t$SeqLnmicroFV\n";
		
		}else{
			
			my @sourceP =split//,$SeqLnmicroFP;
			my @baseP=split//,$SeqLnmicroFV;
			

			# print "$RawSeq\t$SeqLnmicroFP\t$SeqLnmicroFV\t";
			#
			# my ($a_align, $b_align);#


			## callbacks that print something useful
			## prints an 'alignment string' in the order of the  
			## recursion of the dynamic programming algorithm
			## print "-" only on match
			#sub on_align  {  
			#$a_align .= $sourceP[$_[0]];
			#$b_align .= $baseP[$_[1]];
			#};

			#sub on_shift_a { $a_align .= $sourceP[$_[0]]; $b_align .= '-'};
			#sub on_shift_b { $b_align .= $baseP[$_[0]]; $a_align .= '-'};

			### Dumb select, need to return one of the keys for alternative
			### alignments with equal score. Here, we always take the first option, but don't print it.

			#sub on_select_align { print "(select_align)\n"; return (keys (%{$_[0]})) [0]};
			## one gets the same behaviour with not assigning on_select_align, I am too lazy to implement this callback correctly ...

			# my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);
	# 		my $score = $matcher->align(
	# 		                \@sourceP,
	# 		                \@baseP,
	# 		                {   align     => \&on_align,
	# 		                shift_a => \&on_shift_a,
	# 		                shift_b => \&on_shift_b,
	# 		         	    select_align => \&on_select_align
	# 		                });
	##
			 # print "$a_align\t$b_align\t";
			 # 	 			print "$score\n";
			#
			# my @source =split//,$a_align;
			# my @base=split//,$b_align;
			#
			my @MismatchF;
			for my $k (0.. $#sourceP){
				if ($sourceP[$k] ne $baseP[$k]){	
					push @MismatchF,$k;
				}
				
			}
			
			my $maxF= max @MismatchF;

			
	
			$Mismatch_statusF=($#MismatchF>1 )?"Potential":"TRUE";
			
			#print "$#MismatchF\t$Mismatch_statusF\t$a_align\t$b_align\t$score\n";
			
			$start3= $start3+$maxF+1 if ($type eq "+");
 			$end3=$end3-$maxF-1 if ($type eq "-");
			
			
			#print "$insertion{$i}\n" if($#MismatchF >2);
			
			if ($maxF == ($LenmicroFP-1)){
				$LenmicroF=0;
				$SeqLnmicroF="";
				
			}else{
					
				$LenmicroF=$LenmicroFP-$maxF-1;
				$SeqLnmicroF=substr $SeqLnmicroFP,($maxF+1),$LenmicroF;
			}
			
			
		}
		
		
	#### Here is to generate the 5 donor junction
	
		my $DonorFA1 = substr $DonorSeq,($start3+$LenmicroF-16),30;
		
		
		my $DonorFA2 = substr $DonorSeq, ($end3-$LenmicroF-15), 30;	
		
		my $DonorFA2r= reverse $DonorFA2;
		
		$DonorFA2r =~ tr/ATGCatgc/TACGtacg/;
			
		my $FDonorFA= ($type eq "+")?$DonorFA1:$DonorFA2r;
		
		
		$SeqDonorF=$FDonorFA;
		
		
		
	}else{
	### generate the final start point of the insertion	
		$finalStart=41;
	#### Here is generate the 5' end Junction
		my $JunctionFA=	substr $RawSeq,27,15;
		my $JunctionFB= substr $RawSeq,42,15;
		$JunctionF= join '',$JunctionFA,$JunctionFB;
		
	##### The length of 5' deletion/insertion and the sequence;
	
		if($start2<=42){
			$DelmicroF=0;
			$SeqDelmicroF=0;
		}else{
			$DelmicroF=$start2-42;
			my $SeqDelmicroFR = substr $RawSeq, 42, $DelmicroF;
			$SeqDelmicroF="+".$SeqDelmicroFR;
			
		}
	
	######################################################
	#### the length of 5' microhomology and the sequence;\
	#####################################################	
	
		### potential microhomology sequence based on MAT sequence	 
		my $LenmicroFP= (41>$start2)? (41-$start2+1):0;
		my $SeqLnmicroFP = substr $stringS, (41 -$LenmicroFP),$LenmicroFP;
		

		#### potential  microhomology sequence based on donor sequence
		my $SeqLnmicroFV1 = substr $DonorSeq, ($start3-1), $LenmicroFP;
		
		my $SeqLnmicroFV2 = substr $DonorSeq, ($end3-$LenmicroFP),$LenmicroFP;
		
		my $SeqLnmicroFV2r = reverse $SeqLnmicroFV2;
		
		$SeqLnmicroFV2r =~ tr/ATGCatgc/TACGtacg/;
		
		my $SeqLnmicroFV= ($type eq "+")?$SeqLnmicroFV1:$SeqLnmicroFV2r;
		
		
		
		#### confirmation if they are conisistent ####
		
		##### check whether donor microhomology is consistent with mat micrhomoloyg
		### print the $LenmicroF and $SeqLnmicroF as the final microhomology
		
		
		
		if ($SeqLnmicroFP eq $SeqLnmicroFV || (!$SeqLnmicroFP && !$SeqLnmicroFV) ){
			
			$LenmicroF=$LenmicroFP;
			$SeqLnmicroF=$SeqLnmicroFP;
			$Mismatch_statusF="Confidence";
			
		}else{
			
			my @source =split//,$SeqLnmicroFP;
			my @base=split//,$SeqLnmicroFV;
			
			my @MismatchF;
			for my $k (0.. length $SeqLnmicroFP){
				if ($source[$k] ne $base[$k]){	
					push @MismatchF,$k;
				}
				
			}
			
			my $maxF= max @MismatchF;
			
			$start3= $start3+$maxF+1 if ($type eq "+");
 			$end3=$end3-$maxF-1 if ($type eq "-");
			
			
			
			
		#
	# 		my ($a_align, $b_align);#
	#
	#
	# 		## callbacks that print something useful
	# 		## prints an 'alignment string' in the order of the
	# 		## recursion of the dynamic programming algorithm
	# 		## print "-" only on match
	# 		sub on_align  {
	# 		$a_align .= $source[$_[0]];
	# 		$b_align .= $base[$_[1]];
	# 		};
	#
	# 		sub on_shift_a { $a_align .= $source[$_[0]];};
	# 		sub on_shift_b { $b_align .= $base[$_[0]];};
	#
	# 		### Dumb select, need to return one of the keys for alternative
	# 		### alignments with equal score. Here, we always take the first option, but don't print it.
	#
	# 		sub on_select_align { print "(select_align)\n"; return (keys (%{$_[0]})) [0]};
	# 		## one gets the same behaviour with not assigning on_select_align, I am too lazy to implement this callback correctly ...
	#
	# 		my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);
	# 		my $score = $matcher->align(
	# 		                \@source,
	# 		                \@base,
	# 		                {   align     => \&on_align,
	# 		                shift_a => \&on_shift_a,
	# 		                shift_b => \&on_shift_b,
	# 		         #select_align => \&on_select_align
	# 		                });
	# #
	# # 		print "$RawSeq\t$a_align\t$b_align\t";
	# # 		print "$score\n";
	# # #
	#
			
			$Mismatch_statusF=($#MismatchF>1)?"Potential":"TRUE";
			#print "$insertion{$i}\n" if($#MismatchF >2);
			
			if ($maxF == ($LenmicroFP-1)){
				$LenmicroF=0;
				$SeqLnmicroF="";
				
			}else{
					
				$LenmicroF=$LenmicroFP-$maxF-1;
				$SeqLnmicroF=substr $SeqLnmicroFP,($maxF+1),$LenmicroF;
			}
			
			
		}
		
	
	###############################################
	#### Here is to generate the 5 donor junction
	###############################################	
		my $DonorFA1 = substr $DonorSeq,($start3+$LenmicroF-16),30;
		my $DonorFA2 = substr $DonorSeq, ($end3-$LenmicroF-15), 30;
		my $DonorFA2r= reverse $DonorFA2;
		
		$DonorFA2r =~ tr/ATGCatgc/TACGtacg/;
	
		my $FDonorFA= ($type eq "+")?$DonorFA1:$DonorFA2r;
	
		$SeqDonorF=$FDonorFA;
		 
	}
	
	#print OUT "$i\t$start1\t$end1\t$start2\t$end2\t$length\t$chr\t$start3\t$end3\t$type\t$LenmicroF\t$SeqLnmicroF\t$DelmicroF\t$SeqDelmicroF\t$JunctionF\t$SeqDonorF\t";
	
	print OUT "$insertion{$i}\t$LenmicroF\t$SeqLnmicroF\t$Mismatch_statusF\t$DelmicroF\t$SeqDelmicroF\t$JunctionF\t$SeqDonorF\t";
	
	
	
	###### Here is to detect the microhomology in the 3' junction 
	my $rightMAT=$length-$end1+1; my $JunctionR=0; my $DelmicroR=0; my $SeqDelmicroR=0; my $SeqLnmicroR=0;	
	my $LenmicroR=0; 
	my $SeqDonorR=0;
	
	### generate the final end point of the insertion
	my $finalEnd=$end1;
	
	my $Mismatch_statusR="TRUE";
			
	
	if($rightMAT<=47){
		
		#### Here is generate the 3' end Junction
		
		my $JunctionRA=	substr $RawSeq,$end1-16,15;
		my $JunctionRB= substr $RawSeq,$end1-1,15;	
		$JunctionR= join '',$JunctionRA,$JunctionRB;	
		
	
		#### the length of 5' deletion/insertion and the sequence;	 
		
		if($end2>=$end1){ ### here do have microhomology at least 1bp
			$DelmicroR= 47-$rightMAT;
			my $SeqDelmicroRR = substr $stringE, 0 , $DelmicroR;
			$SeqDelmicroR="-".$SeqDelmicroRR;
	
		}else{
				
			my $DelmicroR1=47-$rightMAT;
			
			
			my $InmicroR2=$end1-$end2-1;
			$DelmicroR=$DelmicroR1."\,\+$InmicroR2";
			
			my $SeqDelmicroRR1 = substr $stringE,0,$DelmicroR1;
			my $SeqDelmicroRR2 = substr $RawSeq,$end2,$InmicroR2;
			
			$SeqDelmicroR="\-$SeqDelmicroRR1"."\,\+$SeqDelmicroRR2";
		}
	
		######################################################
		#### the length of 3' microhomology and the sequence;\
		#####################################################	
		#### potenital	3' mcirohomology
		my $LenmicroRP= ($end2>=$end1)? ($end2-$end1+1):0;
		my $SeqLnmicroRP = substr $RawSeq, ($end1-1), $LenmicroRP;
	
		
		
		#### potential  3'microhomology sequence based on donor sequence
		my $SeqLnmicroRV1 = substr $DonorSeq, ($end3-$LenmicroRP),$LenmicroRP;
		
		my $SeqLnmicroRV2 = substr $DonorSeq, ($start3-1), $LenmicroRP;
		
		my $SeqLnmicroRV2r = reverse $SeqLnmicroRV2;
		
		$SeqLnmicroRV2r =~ tr/ATGCatgc/TACGtacg/;
		
		my $SeqLnmicroRV= ($type eq "+")?$SeqLnmicroRV1:$SeqLnmicroRV2r;
		
		
		#### confirmation if they are conisistent ####
		
		##### check whether donor microhomology is consistent with mat micrhomoloyg
		### print the $LenmicroF and $SeqLnmicroF as the final microhomology
		
		if ($SeqLnmicroRP eq $SeqLnmicroRV || (!$SeqLnmicroRP  && !$SeqLnmicroRV)){
			
			$LenmicroR=$LenmicroRP;
			$SeqLnmicroR=$SeqLnmicroRP;
			$Mismatch_statusR="Confidence";
			
		}else{
			
			my @sourceRM =split//,$SeqLnmicroRP;
			my @baseRM=split//,$SeqLnmicroRV;
			
			my @MismatchR;
			for my $j (0.. $#sourceRM){
				if ($sourceRM[$j] ne $baseRM[$j]){	
					push @MismatchR,$j;
				}
				
			}
			
			my $minR= min @MismatchR;
			#print "$insertion{$i}\n" if($#MismatchR >2);
			
		
			# my ($a_align, $b_align);#
			#
			#
			# ## callbacks that print something useful
			# ## prints an 'alignment string' in the order of the
			# ## recursion of the dynamic programming algorithm
			# ## print "-" only on match
			# sub on_align  {
			# $a_align .= $sourceRM[$_[0]];
			# $b_align .= $baseRM[$_[1]];
			# };
			#
			# sub on_shift_a { $a_align .= $sourceRM[$_[0]];};
			# sub on_shift_b { $b_align .= $baseRM[$_[0]];};
			#
			# ### Dumb select, need to return one of the keys for alternative
			# ### alignments with equal score. Here, we always take the first option, but don't print it.
			#
			# sub on_select_align { print "(select_align)\n"; return (keys (%{$_[0]})) [0]};
			# ## one gets the same behaviour with not assigning on_select_align, I am too lazy to implement this callback correctly ...
			#
			# my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);
			# my $score = $matcher->align(
			#                 \@sourceRM,
			#                 \@baseRM,
			#                 {   align     => \&on_align,
			#                 shift_a => \&on_shift_a,
			#                 shift_b => \&on_shift_b,
			#          #select_align => \&on_select_align
			#                 });
			#
			# #
			# #
			# # print "$RawSeq\t$SeqLnmicroRP\t$SeqLnmicroRV\t$a_align\t$b_align\t";
			# # print "$score\n";
			#
			#
			
			$Mismatch_statusR=($#MismatchR>1)?"Potential":"TRUE";
		
			$start3= $start3+$LenmicroRP-$minR if ($type eq "-");
 			$end3=$end3-$LenmicroRP+$minR-1  if ($type eq "+");
			
			if ($minR == 0){
				$LenmicroR=0;
				$SeqLnmicroR="";
				
			}else{
					
				$LenmicroR=$minR;
				$SeqLnmicroR=substr $SeqLnmicroRP,0,$LenmicroR;
			}
			
			
		}
		
		
		
		#### Here is to generate the 3 donor junction
	
			my $DonorRA1 = substr $DonorSeq,($end3-$LenmicroR-15),30;
		
		
			my $DonorRA2 = substr $DonorSeq, ($start3+$LenmicroR-16), 30;	
		
			my $DonorRA2r= reverse $DonorRA2;
		
			$DonorRA2r =~ tr/ATGCatgc/TACGtacg/;
			
			my $FDonorRA= ($type eq "+")?$DonorRA1:$DonorRA2r;
		
		
			$SeqDonorR=$FDonorRA;
		
		
			
	}else{
		
		### generate the final end point of the insertion
		$finalEnd=$length-47+1;
		
		#### Generate the 3' end Junction
		
		$JunctionR= substr $RawSeq,$finalEnd-16,30;
		

		
		#### the length of 3' deletion/insertions and the sequence;

	
			if($end2>=$finalEnd){
				$DelmicroR=0;
				$SeqDelmicroR=0;
			}else{
				$DelmicroR=$finalEnd-$end2-1;
				my $SeqDelmicroRR = substr $RawSeq, $end2, $DelmicroR;
				$SeqDelmicroR="+".$SeqDelmicroRR;
			
			}
	
		
		######################################################
		#### the length of 3' microhomology and the sequence;\
		#####################################################	
		
		#### potenital	3' mcirohomology based on the mata sequence		
			 
		my $LenmicroRP= ($end2>=$end1)? ($end2-$finalEnd+1):0;
		my $SeqLnmicroRP = substr $RawSeq, ($finalEnd-1), $LenmicroRP;
		
		
		
		#### potential  3'microhomology sequence based on donor sequence
		my $SeqLnmicroRV1 = substr $DonorSeq, ($end3-$LenmicroRP), $LenmicroRP;
		
		my $SeqLnmicroRV2 = substr $DonorSeq, ($start3-1),$LenmicroRP;
		
		my $SeqLnmicroRV2r = reverse $SeqLnmicroRV2;
		
		$SeqLnmicroRV2r =~ tr/ATGCatgc/TACGtacg/;
		
		my $SeqLnmicroRV= ($type eq "+")?$SeqLnmicroRV1:$SeqLnmicroRV2r;
		
		
		
		# print "$RawSeq\t$LenmicroRP\t$SeqLnmicroRP\t$SeqLnmicroRV\n";
		
		#### confirmation if they are conisistent ####
		
		##### check whether donor microhomology is consistent with mat micrhomoloyg
		### print the $LenmicroF and $SeqLnmicroF as the final microhomology
		
		if ($SeqLnmicroRP eq $SeqLnmicroRV ||(!$SeqLnmicroRP  && !$SeqLnmicroRV)){
			
			$LenmicroR=$LenmicroRP;
			$SeqLnmicroR=$SeqLnmicroRP;
			$Mismatch_statusR="Confidence";
			
			
		}else{
			
			my @sourceRN =split//,$SeqLnmicroRP;
			my @baseRN=split//,$SeqLnmicroRV;
			
			my @MismatchR;
			for my $j (0 .. length $SeqLnmicroRP){
				if ($sourceRN[$j] ne $baseRN[$j]){	
					push @MismatchR,$j;
				}
				
			}
			
			my $minR= min @MismatchR;		
			$start3= $start3+$LenmicroRP-$minR if ($type eq "-");
 			$end3=$end3-$LenmicroRP+$minR-1 if ($type eq "+");
			
			#
		# 	my ($a_align, $b_align);#
		#
		#
		# 	## callbacks that print something useful
		# 	## prints an 'alignment string' in the order of the
		# 	## recursion of the dynamic programming algorithm
		# 	## print "-" only on match
		# 	sub on_align  {
		# 	$a_align .= $sourceRN[$_[0]];
		# 	$b_align .= $baseRN[$_[1]];
		# 	};
		#
		# 	sub on_shift_a { $a_align .= $sourceRN[$_[0]];};
		# 	sub on_shift_b { $b_align .= $baseRN[$_[0]];};
		#
		# 	### Dumb select, need to return one of the keys for alternative
		# 	### alignments with equal score. Here, we always take the first option, but don't print it.
		#
		# 	sub on_select_align { print "(select_align)\n"; return (keys (%{$_[0]})) [0]};
		# 	## one gets the same behaviour with not assigning on_select_align, I am too lazy to implement this callback correctly ...
		#
		# 	my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);
		# 	my $score = $matcher->align(
		# 	                \@sourceRN,
		# 	                \@baseRN,
		# 	                {   align     => \&on_align,
		# 	                shift_a => \&on_shift_a,
		# 	                shift_b => \&on_shift_b,
		# 	         #select_align => \&on_select_align
		# 	                });
		#
		# 	# print "$RawSeq\t$SeqLnmicroRP\t$SeqLnmicroRV\t$a_align\t$b_align\t";
		# # 	print "$score\n";

		
			
			$Mismatch_statusR=($#MismatchR>1)?"Potential":"TRUE";
		
			print "$RawSeq\t$LenmicroRP\t$SeqLnmicroRP\t$SeqLnmicroRV\t$Mismatch_statusR\t$#MismatchR\n";
			#print "$insertion{$i}\n" if($#MismatchR >2);
			
			if ($minR == 0){
				$LenmicroR=0;
				$SeqLnmicroR="";
				
			}else{
					
				$LenmicroR=$minR;
				$SeqLnmicroR=substr $SeqLnmicroRP,0,$LenmicroR;
			}	
			
		}
		
		#### Here is to generate the 3 donor junction		
		my $DonorRA1 = substr $DonorSeq,($end3-$LenmicroR-15),30;
		my $DonorRA2 = substr $DonorSeq, ($start3+$LenmicroR-16), 30;
		my $DonorRA2r= reverse $DonorRA2;
		
		$DonorRA2r =~ tr/ATGCatgc/TACGtacg/;
	
		my $FDonorRA= ($type eq "+")?$DonorRA1:$DonorRA2r;
	
		$SeqDonorR=$FDonorRA;
		
	}

	my $finalLength=$finalEnd-$finalStart-1;
	my $final_insertion=substr $RawSeq, $finalStart,$finalLength;
	
	print OUT "$LenmicroR\t$SeqLnmicroR\t$Mismatch_statusR\t$DelmicroR\t$SeqDelmicroR\t$JunctionR\t$SeqDonorR\t$finalLength\t$final_insertion\n";
}

close OUT;

