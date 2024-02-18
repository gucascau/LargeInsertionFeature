#!/bin/bash
##################################################################
# Author: Xin Wang                                                   
# Email: xin.wang@childrens.harvard.edu                            
# Copyright (c) 2022 Dr. Kaifu lab                                   
# PI: Kaifu Chen                                                   
# Description: 
#  
#   The pipeline were used for the proximity and overlapping feature analyses of large insertion events. These features include ARS, Telomere, R-loop, Centromere, Tandem repeats, tRNA
#	Method Notes:
# 		1. We ignored the insertion events from MATa, LTR, rDNA and Mitochondrion, mostly concentrated on the genomic region.
#		2. 
#	
################################################################





### Usage function tells users how to run the software
helpFunction()
{
	echo "*********************************** how to use Ty1Nuceotide ***********************************"
	echo "Usage: sh $0 -a Sample ID -b Work Directory  -f Single insertion events -r poly insertion events -o Output  -p Software installed Directory"
	echo "Note: In certain mutants and aging cells, we have identified a substantial number of insertions originating from retrotransposons. However, these retrotransposons were found to be extensively distributed throughout the entire genome, posing a challenge in assessing the mechanism of retrotransposon insertion. To delve into the Ty1 nucleotide insertion frequency within each mutant, we have created this package to comprehensively analyze the distribution of insertions across the Ty1 element. This tool aims to provide a detailed understanding of the Ty1 insertion patterns, enabling a more precise investigation into the mechanisms underlying retrotransposon insertions in various mutants and aging cell contexts."
	echo ""
	echo -e "Request Parameters:"
	echo -e "\t-a Sample Id (Example: yYY398-B_S10)"
	echo -e "\t-b The working directory, where you put the insertion events"
	echo -e "\t-f The insertion event (Note: We only considered the insertion defined as single donor as the insertions with multiple donors might not be able to have a clean Ty1 locus.)"
	echo -e "\t-p Software installed Path (Default:"", it is the path where you store your blast software)" 

	
	echo ""
	echo ""
	echo -e "Optional Parameters:"
	echo "" 
	echo -e "Optional Parameters -- Ty1 Blast setting up:"
	echo -e "\t-ms identity (Default: 60)"
	echo -e "\t-mc Gapsize  (Default t)"
	echo -e "\t-mb E value  (Default t)"
	echo "" 
	echo -e "Optional Parameters of Ty1 reference"
	echo -e "\t-c Ty1 blast index (Note: if you don't have blast index, please run the 'makeblastdb -in Ty1.fasta -dbtype nucl') (Default: ./LargeInsertionFeature/Ty1Nuceotide/Database/YPLWTY1-1.fasta)"
	echo -e "\t-d The Ty1 annotation file (Default: ./LargeInsertionFeature/Ty1Nuceotide/Database/YPLWTY1-1.annotation.bed)"
	echo -e "\t-r Script Stored Path (Default: ./LargeInsertionFeature/Ty1Nuceotide/src)" 
	echo ""

	echo -e "\t-h help"
	
	echo "For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu"
	echo "**************************************************"
   exit 1 # Exit script after printing help
}


# Get absolute path for scripts and check if required scripts exist
in=$PWD

### Set default options
# some general paramaters including mata locus thread, and software path

# default software path 
softwarepath=''

# default parameter parameter for selfblast
### 
DepIden=60
DepGapsize=5
DeEvalue=0.01


while getopts "a:b:c:d:p:f:r:ms:mc:mb" opt
do
   case "$opt" in
      a ) SampleID="$OPTARG" ;;
      b ) WDir="$OPTARG" ;;
	  p ) softpath="$OPTARG" ;;
  	  f ) Insert="$OPTARG" ;;


      ms ) DepIden="$OPTARG" ;;
	  mc ) DepGapsize="$OPTARG" ;;
      mb ) DeEvalue="$OPTARG" ;;
	  c ) TYIndex="$OPTARG" ;;
      d ) TYAnn="$OPTARG" ;;
	  r ) Dscript="$OPTARG" ;;

      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


###single insertion:

## updated the name
#cp wt 3 days.One.txt ${SampleID}.One.txt

# default script path

srcDir=${softpath}/LargeInsertionFeature/Ty1Nuceotide/src

# default Ty1-1 reference sequence
Ty1=${softpath}/LargeInsertionFeature/Ty1Nuceotide/Database/YPLWTY1-1.fasta

# default Ty1-1 reference annotation
Ty1Ann=${softpath}/LargeInsertionFeature/Ty1Nuceotide/Database/YPLWTY1-1.annotation.bed

# Print helpFunction in case parameters are empty
if [ -z "${SampleID}" ] 
then
   echo "*** error: input Sample ID must be provided ***";
   helpFunction
fi

if [ -z "${WDir}" ] 
then
   echo "*** error: input work path must be provided ***";
   helpFunction
fi


if [ -z "${Insert}" ] 
then
   echo "*** error: input Insertion events with Retrotransposons must be defined ***";
   helpFunction
fi


if [ -z "${softpath}" ] 
then
   echo "*** error: input software path must be provided ***";
   helpFunction
fi


# Begin script in case all parameters are correct
echo "${SampleID}"
echo "${WDir}"
echo "${TYIndex}"
echo "${TYAnn}"
echo "${Insert}"
echo "${Dscript}"
echo "${softpath}"
echo "${Ty1}"
echo "${Ty1Ann}"
echo "${DepIden}"
echo "${DepGapsize}"
echo "${DeEvalue}"
# 

echo "All the paramter are sucessfully provided, now let's detect the large insertion event"


#### Here is the last version of large insertion events


echo "The job array is started ..."
date




################################################
### pre-process the insertion events, 
################################################

## remove MAT, rDNA, chrmt, chrmp,MATa hotspots, LTR 
cat  wt*.One.txt wt*MultipleClean.txt| perl -ne '{chomp; my ($represent,$id,$chr,$start,$end,$strand,$inf,$ltr,$distdes)=(split/\t/,$_)[1,2,10,11,12,13,41,42,43]; next if ($inf eq "NO" ); next if ($chr eq "Unknown"); next if ($chr eq "chrXII" && $start >=451418 && $end <= 469316); next if ($chr eq "chrmp"||$chr eq "chrmt" ); my $type=(split/\|/,$inf)[0];next if (($type eq "LTR_retrotransposon" && $distdes eq "ENTIRE") || ($type eq "long_terminal_repeat" && $distdes eq "ENTIRE")); next if (($type eq "silent_mating_type_cassette_array" ) || ($type eq "mating_type_region")); next if ($id eq "SampleID"); $hash{$type}->{$id}++; $str{$id}++; print "$chr\t$start\t$end\t$represent\t0\t$strand\t$inf\n"}'|sort -k 1,1 -k 2,2n  >WT_Combined.noMATaTyrDNAChrmtChrmp.bed

## generate the randomized insertion events, here we excluded the LTR, MAT, rDNA and Mitochondrion for the randomization.
for i in {1..1000}; do ${softpath}/bedtools2/bin/bedtools shuffle -i WT_Combined.noMATaTyrDNAChrmtChrmp.bed -g saccharomyces_cerevisiae_R64-2-1_20150113_modified.genome.txt -chrom  -excl LTR_MAT_rDNA_Mit_sorted.bed |sort -k 1V,1 -k 2n,2 > WTAging_Combined_randome_${i}.bed; done


#################################################################
### check the randomized insertion and the closest features, 
###    including ARS, Telomere, rloop, tandem repeat, Centromere and tRNA
#################################################################

# 1. ARS

# For the observed results
~/Software/bedtools2/bin/bedtools  closest -D b -t first -b ./Database/ARS_distribution_modified2.txt -a WT_Combined.noMATaTyrDNAChrmtChrmp.bed  >WT_Combined.noMATaTyrDNAChrmtChrmp.ARS.bed

# The number that are close to the ARS

awk '{if ($12<=1000 &&$12>=-1000){print}}'   WT_Combined.noMATaTyrDNAChrmtChrmp.ARS.bed |wc (1723)

# For the random results

for i in {1..1000}; do ~/Software/bedtools2/bin/bedtools  closest -D b -t first -b .././Database/ARS_distribution_modified2.txt -a ../../RandomInsertion/WTAging_Combined_randome_${i}.bed >Insertion.${i}.ARSanot.bed; done

for i in {1..1000}; do awk '{if ($12<=1000 &&$12>=-1000){n++}} END {print n}' Insertion.${i}.ARSanot.bed >Insertion.${i}.ARSanot.number.txt; done

# generated the random number that appreciate to ARS
cat Insertion.*.ARSanot.number.txt > WT_Combined.Random.ARSannt.number.txt

### Final combined and measure the P value
perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t1723\n"}else{print "Random\t$_\n"}}' WT_Combined.Random.ARSannt.number.txt > WT_Combined.FinalARSComparison.number.txt
cp WT_Combined.* ../

# 2. Telomere

# observed
${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ./Database/Telomere.bed -a  WT_Combined.noMATaTyrDNAChrmtChrmp.bed >WT_Combined.noMATaTyrDNAChrmtChrmp.Telomore.bed
awk '{if ($15<=1000 &&$15>=-1000){print}}'  WT_Combined.noMATaTyrDNAChrmtChrmp.Telomore.bed|wc (243)

## Random
mkdir Random
cd Random
# annotated with the closest telomere element for WT aging
for i in {1..1000}; do ${softpath}/bedtools2/bin/bedtools closest -D b -t first -b .././Database/Telomere.bed -a ../../RandomInsertion/WTAging_Combined_randome_${i}.bed >Insertion.${i}.Telomer.bed; done
# check the number of insertion events that are close to telmore (1kbp)

for i in {1..1000}; do awk '{if ($15<=1000 &&$15>=-1000){n++}} END {print n}' Insertion.${i}.Telomer.bed >Insertion.${i}.Telomer.number.txt; done
cat Insertion*Telomer.number.txt > WT_Combined.Ranome.Telemore.number.txt

## generate a comparision between observed and random insertions
perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t243\n"}else{print "Random\t$_\n"}}'  WT_Combined.Ranome.Telemore.number.txt > WT_Combined.FinalTelomereComparison.number.txt


# 3. For rloop

#observed
${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ./Database/Rloop_Final.bed -a ../WT_Combined.noMATaTyrDNAChrmtChrmp.bed > WT_Combined.noMATaTyrDNAChrmtChrmp.Rloop.bed

awk '{if ($11>=-200 && $11<=200){n++}} END {print n}'  WT_Combined.noMATaTyrDNAChrmtChrmp.Rloop.bed (788)
#Random
for i in {1..1000}; do ${softpath}/bedtools2/bin/bedtools closest -D b -t first -b .././Database/Rloop_Final.bed -a ../../RandomInsertion/WTAging_Combined_randome_${i}.bed >Insertion.${i}.Rloop.bed; done

# 200bp

for i in {1..1000}; do awk '{if ($11>=-200 && $11<=200){n++}} END {print n}' Insertion.${i}.Rloop.bed >Insertion.${i}.Rloop200Close.number.txt; done
cat Insertion*.Rloop200Close.number.txt > WT_Combined.Random.Rloop200Close.number.txt

# compare between Random and observed
perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t788\n"}else{print "Random\t$_\n"}}' WT_Combined.Random.Rloop200Close.number.txt  > WT_Combined.Rloop200Close.numberCom.txt



# 4. For Centromere
# observed
${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ./Database/Centromere.bed -a ../WT_Combined.noMATaTyrDNAChrmtChrmp.bed> WT_Combined.noMATaTyrDNAChrmtChrmp.Centromere.bed
awk '{if ($15<=10000 &&$15>=-10000){print}}'  WT_Combined.noMATaTyrDNAChrmtChrmp.Centromere.bed

# random
for i in {1..1000}; do ${softpath}/bedtools2/bin/bedtools closest -D b -t first -b .././Database/Centromere.bed  -a ../../RandomInsertion/WTAging_Combined_randome_${i}.bed >Insertion.${i}.Centro.bed; done
for i in {1..1000}; do awk '{if ($15<=10000 &&$15>=-10000){n++}} END {print n}' Insertion.${i}.Centro.bed>Insertion.${i}.Centro10k.number.txt; done
cat Insertion*.Centro10k.number.txt > WT_Combined.Random.Centro10k.number.txt

perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t265\n"}else{print "Random\t$_\n"}}'  WT_Combined.Random.Centro10k.number.txt > WT_Combined.Centro10k. numberCom.txt

## 
## 5. For Tandem repeats
# observed
${softpath}/bedtools2/bin/bedtools closest -D b -t first -b  ./Database/TRF.bed -a ../WT_Combined.noMATaTyrDNAChrmtChrmp.bed> WT_Combined.noMATaTyrDNAChrmtChrmp.TRF.bed

awk '{if ($24==0){n++}} END {print n}' WT_Combined.noMATaTyrDNAChrmtChrmp.TRF.bed (226)

## extract insertions with Tandem repeats


# random
for i in {1..1000}; do ${softpath}/bedtools2/bin/bedtools closest -D b -t first -b  .././Database/TRF.bed -a ../../RandomInsertion/WTAging_Combined_randome_${i}.bed >Insertion.${i}.TRF.bed; done

for i in {1..1000}; do awk '{if ($24==0){n++}} END {print n}' Insertion.${i}.TRF.bed >Insertion.${i}.TRF.number.txt; done

cat Insertion.*TRF.number.txt |perl -ne '{chomp;my $i=$_; if ($i eq ""){print "0\n"}else{print "$i\n"}}' > WT_Combined.Random.TRF.number.txt

perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t226\n"}else{print "Random\t$_\n"}}'   WT_Combined.Random.TRF.number.txt > WT_Combined.TRF.numberComp.txt



## 6.  tRNA 

# observed
${softpath}/bedtools2/bin/bedtools closest -D b -t first -b  ./Database/S288C_tRNA.nonchrM.bed -a ../WT_Combined.noMATaTyrDNAChrmtChrmp.bed> WT_Combined.noMATaTyrDNAChrmtChrmp.tRNA.bed
awk '{if ($15<=200 &&$15>=-200){n++}} END {print n}' WT_Combined.noMATaTyrDNAChrmtChrmp.tRNA.bed (202)


## random
for i in {1..1000}; do ${softpath}/bedtools2/bin/bedtools closest -D b -t first -b  .././Database/S288C_tRNA.nonchrM.bed -a ../../RandomInsertion/WTAging_Combined_randome_${i}.bed >Insertion.${i}.tRNA.bed; done

for i in {1..1000}; do awk '{if ($15<=200 &&$15>=-200){n++}} END {print n}' Insertion.${i}.tRNA.bed > Insertion.${i}.tRNA.number.txt ; done;

cat *.number.txt > WT_Combined.Random.tRNA.txt
# perl -ne '{chomp; my $num=($_>0)?$_:0; print "$num\n"}' WT_Combined.Random.tRNA.txt >Rad27_Random_updated.tRNA.txt
perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t202\n"}else{print "Random\t$_\n"}}'   WT_Combined.Random.tRNA.txt > WT_Combined.tRNA.numberComp.txt





### Set up path file:

echo "Change to the Working Path, where you store your Large events"

echo "Congratulation! Insertion detection and deduplication job array is finished !"