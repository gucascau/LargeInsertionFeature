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
#		2. You need to update the observed number for the comparision to randamized number.
#		3. You can also set up the requested upstream/downstream approximity distance and random number accordingly.
#		4. You can add other sorted features for comparision. i.e., LTR, Proteincoding, UTR et.al.
################################################################





### Usage function tells users how to run the software
helpFunction()
{
	echo "*********************************** how to use InsFeatures ***********************************"
	echo "Usage: sh $0 -a Sample ID -b Work Directory  -f Single insertion events -r poly insertion events -o Output  -p Software installed Directory"
	echo "Note: In certain mutants and aging cells, we have identified a substantial number of insertions originating from retrotransposons. However, these retrotransposons were found to be extensively distributed throughout the entire genome, posing a challenge in assessing the mechanism of retrotransposon insertion. To delve into the Ty1 nucleotide insertion frequency within each mutant, we have created this package to comprehensively analyze the distribution of insertions across the Ty1 element. This tool aims to provide a detailed understanding of the Ty1 insertion patterns, enabling a more precise investigation into the mechanisms underlying retrotransposon insertions in various mutants and aging cell contexts."
	echo ""
	echo -e "Request Parameters:"
	echo -e "\t-a Genotype name(Example: WT)"
	echo -e "\t-b The working directory, where you put the insertion features"
	echo -e "\t-f The single donor insertion event generated by iDSBins"
	echo -e "\t-r The poly-donor insertion event generated by iDSBins"
	echo -e "\t-p Software installed Path (Default:"", it is the path where you store your bedtools software)" 

	
	echo ""
	echo ""
	echo -e "Optional Parameters:"
	
	echo ""
	echo -e "\t-d The yeast feature folder (Default: ./LargeInsertionFeature/InsFeatures/Database/)"
	echo -e "\t-r Script Stored Path (Default: ./LargeInsertionFeature/InsFeatures/src)" 
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



while getopts "a:b:p:f:r:d:r" opt
do
   case "$opt" in
      a ) SampleID="$OPTARG" ;;
      b ) WDir="$OPTARG" ;;
	  p ) softpath="$OPTARG" ;;
  	  f ) SInsert="$OPTARG" ;;
  	  r) MIsert="$OPTARG" ;;



	  d ) DataIndex="$OPTARG" ;;
	  r ) Dscript="$OPTARG" ;;

      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


###single insertion:

## updated the name
#cp wt 3 days.One.txt ${GenotypeID}.One.txt

# default script path

srcDir=${softpath}/LargeInsertionFeature/InsFeatures/src

# default R-loop datasets
Rloop=${softpath}/LargeInsertionFeature/InsFeatures/Database/Rloop_Final.bed

# default ARS datasets
ARS=${softpath}/LargeInsertionFeature/InsFeatures/Database/ARS_distribution_modified2.bed

# default Centromere datasets
Centromere=${softpath}/LargeInsertionFeature/InsFeatures/Database/Centromere.bed

# default Telomere datasets
Telomere=${softpath}/LargeInsertionFeature/InsFeatures/Database/Telomere.bed
# default tRNA datasets

tRNA=${softpath}/LargeInsertionFeature/InsFeatures/Database/TRF.bed

# default whole chromosome 
Chromosome =${softpath}/LargeInsertionFeature/InsFeatures/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.genome.txt

# default 
Redudant=${softpath}/LargeInsertionFeature/InsFeatures/Database/LTR_MAT_rDNA_Mit_sorted.bed


# Print helpFunction in case parameters are empty
if [ -z "${GenotypeID}" ] 
then
   echo "*** error: input Sample ID must be provided ***";
   helpFunction
fi

if [ -z "${WDir}" ] 
then
   echo "*** error: input work path must be provided ***";
   helpFunction
fi


if [ -z "${SInsert}" ] 
then
   echo "*** error: input single Insertion events with Retrotransposons must be defined ***";
   helpFunction
fi

if [ -z "${MInsert}" ] 
then
   echo "*** error: input poly Insertion events with Retrotransposons must be defined ***";
   helpFunction
fi


if [ -z "${softpath}" ] 
then
   echo "*** error: input software path must be provided ***";
   helpFunction
fi


# Begin script in case all parameters are correct
echo "${GenotypeID}"
echo "${WDir}"
echo "${SInsert}"
echo "${MInsert}"
echo "${softpath}"

echo "${DataIndex}"

# 

echo "All the paramter are sucessfully provided, now let's detect the large insertion event"


#### Here is the last version of large insertion events


echo "The job array is started ..."
date


### Set up path file:

echo "Change to the Working Path, where you store your Large events"


################################################
### pre-process the insertion events, 
################################################

## remove MAT, rDNA, chrmt, chrmp,MATa hotspots, LTR 
cd ${WDir}

cat  ${Sinsert} ${MIsert}| perl -ne '{chomp; my ($represent,$id,$chr,$start,$end,$strand,$inf,$ltr,$distdes)=(split/\t/,$_)[1,2,10,11,12,13,41,42,43]; next if ($inf eq "NO" ); next if ($chr eq "Unknown"); next if ($chr eq "chrXII" && $start >=451418 && $end <= 469316); next if ($chr eq "chrmp"||$chr eq "chrmt" ); my $type=(split/\|/,$inf)[0];next if (($type eq "LTR_retrotransposon" && $distdes eq "ENTIRE") || ($type eq "long_terminal_repeat" && $distdes eq "ENTIRE")); next if (($type eq "silent_mating_type_cassette_array" ) || ($type eq "mating_type_region")); next if ($id eq "SampleID"); $hash{$type}->{$id}++; $str{$id}++; print "$chr\t$start\t$end\t$represent\t0\t$strand\t$inf\n"}'|sort -k 1,1 -k 2,2n  >${WDir}/${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.bed

mkdir -p ${WDir}/RandomInsertion
cd ${WDir}/RandomInsertion
## generate the randomized insertion events, here we excluded the LTR, MAT, rDNA and Mitochondrion for the randomization.
for i in {1..10000}; do ${softpath}/bedtools2/bin/bedtools shuffle -i ../${WDir}/${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.bed -g ${Chromosome} -chrom  -excl ${Redudant} |sort -k 1V,1 -k 2n,2 > ${GenotypeID}_Combined_randome_${i}.bed; done


#################################################################
### check the randomized insertion and the closest features, 
###    including ARS, Telomere, rloop, tandem repeat, Centromere and tRNA
#################################################################

# 1. ARS

cd ${WDir}
mkdir -p ${WDir}/ARS

cd ${WDir}/ARS

# For the observed results
${softpath}/bedtools2/bin/bedtools  closest -D b -t first -b ${ARS} -a ${WDir}/${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.bed  >${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.ARS.bed

# The number that are close to the ARS

awk '{if ($12<=1000 &&$12>=-1000){print}}'   ${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.ARS.bed  #(1723)

# For the random results
mkdir -p ${WDir}/ARS/Random
cd ${WDir}/ARS/Random

for i in {1..10000}; do ${softpath}/bedtools2/bin/bedtools  closest -D b -t first -b ${ARS} -a ${WDir}/RandomInsertion/${GenotypeID}_Combined_randome_${i}.bed >Insertion.${i}.ARSanot.bed; done

for i in {1..10000}; do awk '{if ($12<=1000 &&$12>=-1000){n++}} END {print n}' Insertion.${i}.ARSanot.bed >Insertion.${i}.ARSanot.number.txt; done

# generated the random number that appreciate to ARS

cat Insertion.*.ARSanot.number.txt > ${GenotypeID}_Combined.Random.ARSannt.number.txt

### Final combined and measure the P value

# attention: please write the correct one for the observed number: 1723
perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t1723\n"}else{print "Random\t$_\n"}}' ${GenotypeID}_Combined.Random.ARSannt.number.txt > ${GenotypeID}_Combined.FinalARSComparison.number.txt
cp ${GenotypeID}_Combined.* ../

# 2. Telomere

# observed

cd ${WDir}
mkdir -p ${WDir}/Telomere

cd ${WDir}/Telomere

${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ./Database/Telomere.bed -a  ${WDir}/${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.bed >${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.Telomore.bed
awk '{if ($15<=1000 &&$15>=-1000){print}}'  ${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.Telomore.bed  #(243)

## Random
mkdir -p ${WDir}/Telomere/Random
cd ${WDir}/Telomere/Random
# annotated with the closest telomere element for ${GenotypeID} aging
for i in {1..10000}; do ${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ${Telomere} -a ${WDir}/RandomInsertion/${GenotypeID}_Combined_randome_${i}.bed >Insertion.${i}.Telomer.bed; done
# check the number of insertion events that are close to telmore (1kbp)

for i in {1..10000}; do awk '{if ($15<=1000 &&$15>=-1000){n++}} END {print n}' Insertion.${i}.Telomer.bed >Insertion.${i}.Telomer.number.txt; done
cat Insertion*Telomer.number.txt > ${GenotypeID}_Combined.Ranome.Telemore.number.txt

## generate a comparision between observed and random insertions

# attention: please write the correct one for the observed number: 243
perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t243\n"}else{print "Random\t$_\n"}}'  ${GenotypeID}_Combined.Ranome.Telemore.number.txt > ${GenotypeID}_Combined.FinalTelomereComparison.number.txt
cp ${GenotypeID}_Combined.* ../

# 3. For rloop
mkdir -p ${WDir}/Rloop
cd ${WDir}/Rloop

#observed
${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ./Database/Rloop_Final.bed -a ${WDir}/${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.bed > ${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.Rloop.bed

awk '{if ($11>=-200 && $11<=200){n++}} END {print n}'  ${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.Rloop.bed   #(788)

## Random
mkdir -p ${WDir}/Rloop/Random
cd ${WDir}/Rloop/Random
for i in {1..10000}; do ${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ${Rloop} -a ${WDir}/RandomInsertion/${GenotypeID}_Combined_randome_${i}.bed >Insertion.${i}.Rloop.bed; done

# 200bp
for i in {1..10000}; do awk '{if ($11>=-200 && $11<=200){n++}} END {print n}' Insertion.${i}.Rloop.bed >Insertion.${i}.Rloop200Close.number.txt; done
cat Insertion*.Rloop200Close.number.txt > ${GenotypeID}_Combined.Random.Rloop200Close.number.txt

# compare between Random and observed

# attention: please write the correct one for the observed number: 788
perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t788\n"}else{print "Random\t$_\n"}}' ${GenotypeID}_Combined.Random.Rloop200Close.number.txt  > ${GenotypeID}_Combined.Rloop200Close.numberCom.txt



# 4. For Centromere
# observed
mkdir -p ${WDir}/Centromere
cd ${WDir}/Centromere

${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ${Centromere} -a ${WDir}/${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.bed> ${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.Centromere.bed
awk '{if ($15<=10000 &&$15>=-10000){print}}'  ${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.Centromere.bed  #(265)

mkdir -p ${WDir}/Centromere/Random
cd ${WDir}/Centromere/Random
# random
for i in {1..10000}; do ${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ${Centromere}  -a ${WDir}/RandomInsertion/${GenotypeID}_Combined_randome_${i}.bed >Insertion.${i}.Centro.bed; done
for i in {1..10000}; do awk '{if ($15<=10000 &&$15>=-10000){n++}} END {print n}' Insertion.${i}.Centro.bed>Insertion.${i}.Centro10k.number.txt; done
cat Insertion*.Centro10k.number.txt > ${GenotypeID}_Combined.Random.Centro10k.number.txt


# attention: please write the correct one for the observed number: 265

perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t265\n"}else{print "Random\t$_\n"}}'  ${GenotypeID}_Combined.Random.Centro10k.number.txt > ${GenotypeID}_Combined.Centro10k. numberCom.txt

## 
## 5. For Tandem repeats
# observed
mkdir -p ${WDir}/TRF
cd ${WDir}/TRF
${softpath}/bedtools2/bin/bedtools closest -D b -t first -b  ${Trf} -a ${WDir}/${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.bed> ${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.TRF.bed

awk '{if ($24==0){n++}} END {print n}' ${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.TRF.bed #(226)

## extract insertions with Tandem repeats
mkdir -p ${WDir}/TRF/Random
cd ${WDir}/TRF/Random

# random
for i in {1..10000}; do ${softpath}/bedtools2/bin/bedtools closest -D b -t first -b  ${Trf} -a ${WDir}/RandomInsertion/${GenotypeID}_Combined_randome_${i}.bed >Insertion.${i}.TRF.bed; done

for i in {1..10000}; do awk '{if ($24==0){n++}} END {print n}' Insertion.${i}.TRF.bed >Insertion.${i}.TRF.number.txt; done

cat Insertion.*TRF.number.txt |perl -ne '{chomp;my $i=$_; if ($i eq ""){print "0\n"}else{print "$i\n"}}' > ${GenotypeID}_Combined.Random.TRF.number.txt

# attention: please write the correct one for the observed number: 226

perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t226\n"}else{print "Random\t$_\n"}}'   ${GenotypeID}_Combined.Random.TRF.number.txt > ${GenotypeID}_Combined.TRF.numberComp.txt



## 6.  tRNA 

# observed
mkdir -p ${WDir}/tRNA
cd ${WDir}/tRNA

${softpath}/bedtools2/bin/bedtools closest -D b -t first -b  ${tRNA} -a ${WDir}/${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.bed> ${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.tRNA.bed
awk '{if ($15<=200 &&$15>=-200){n++}} END {print n}' ${GenotypeID}_Combined.noMATaTyrDNAChrmtChrmp.tRNA.bed  #(202)

## random
mkdir -p ${WDir}/tRNA/Random
cd ${WDir}/tRNA/Random
for i in {1..10000}; do ${softpath}/bedtools2/bin/bedtools closest -D b -t first -b  ${tRNA} -a ${WDir}/RandomInsertion/${GenotypeID}_Combined_randome_${i}.bed >Insertion.${i}.tRNA.bed; done

for i in {1..10000}; do awk '{if ($15<=200 &&$15>=-200){n++}} END {print n}' Insertion.${i}.tRNA.bed > Insertion.${i}.tRNA.number.txt ; done;

cat *.number.txt > ${GenotypeID}_Combined.Random.tRNA.txt
# perl -ne '{chomp; my $num=($_>0)?$_:0; print "$num\n"}' ${GenotypeID}_Combined.Random.tRNA.txt >Rad27_Random_updated.tRNA.txt

# attention: please write the correct one for the observed number: 202
perl -ne '{chomp; my $num=$_; $n++; if ($n==1){print "Type\tTotalNumber\n"}elsif($n==2){print "Observed\t202\n"}else{print "Random\t$_\n"}}'   ${GenotypeID}_Combined.Random.tRNA.txt > ${GenotypeID}_Combined.tRNA.numberComp.txt


echo "Congratulation! Insertion feature job array has finished !"