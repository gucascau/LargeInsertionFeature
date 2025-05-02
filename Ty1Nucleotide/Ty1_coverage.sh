#!/bin/bash
##################################################################
# Author: Xin Wang                                                   
# Email: xin.wang@childrens.harvard.edu                            
# Copyright (c) 2022 Dr. Kaifu lab                                   
# PI: Kaifu Chen                                                   
# Description: 
#  Retrotransposon -- Ty1 nucleotide insertion map analysis
#    In certain mutants and aging cells, we have identified a substantial number of insertions originating from retrotransposons. However, these retrotransposons were found to be extensively distributed throughout the entire genome, posing a challenge in assessing the mechanism of retrotransposon insertion. To delve into the Ty1 nucleotide insertion frequency within each mutant, we have created this package to comprehensively analyze the distribution of insertions across the Ty1 element. This tool aims to provide a detailed understanding of the Ty1 insertion patterns, enabling a more precise investigation into the mechanisms underlying retrotransposon insertions in various mutants and aging cell contexts.

# 	1. We extracted all insertion events that are approximately to or locataed within retrotransposons. 
#	2. We only took one of the longest Ty1-1 from yeast genome in the table, which is the most dominant Ty and easy for the mapping
#	3. We mapped the insertions with Ty retrotransposon annotaton against the Ty1-1 reference
# 		Request: 
#			If there is mismatch between the query sequence and subject sequence during blast, we use the whole subject sequence to do mapping.
# 			If they were mapped to both LTR, use the right LTR because it was reported that cDNA is generated from the right one. However, there are some extra nucleotides (1-3) which are perfectly match the left one. 
#			 However, it could be misleading to be mapped to the left one considering the mechanisms how cDNA is generated from Ty.
#	
################################################################


### Usage function tells users how to run the software
helpFunction()
{
	echo "*********************************** how to use Ty1Nuceotide ***********************************"
	echo "Usage: sh $0 -a Sample ID -b Work Directory  -f Insertion events -r Output  -p Software installed Directory"
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


### Set up path file:

echo "Change to the Working Path, where you store your Large events"

cd ${WDir}


### Workflow for Ty1 analyses ####

# if you hav not done the blast index, please run it before.

#source /programs/biogrids.shrc
#/home/ch220812/software/ncbi-blast-2.8.1+/bin/makeblastdb -in Retrotransposon_Seq_Yang.fasta -dbtype nucl

# ~/Software/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in YPLWTY1-1.fasta


# get the Ty1 insertion fasta files
cat ${WDir}/${Insert}| perl -ne '{chomp; my ($read,$id,$chr,$start,$end,$strand,$seq,$inf,$ltr,$distdes)=(split/\t/,$_)[1,2,10,11,12,13,40,41,42,43];next if ($inf eq "NO");   next if ($chr eq "chrXII" && $start >=451418 && $end <= 469316); my ($type,$details)=(split/\|/,$inf)[0,1];next unless (($type eq "LTR_retrotransposon" && $distdes eq "ENTIRE" && $details =~/Ty1-1/)|($type eq "long terminal repeats")); print ">$read\n$seq\n" }'  > ${SampleID}.LTR.fasta
 
# print the insertion events:

cat ${WDir}/${Insert}| perl -ne '{chomp; my ($read,$id,$chr,$start,$end,$strand,$seq,$inf,$ltr,$distdes)=(split/\t/,$_)[1,2,10,11,12,13,40,41,42,43];next if ($inf eq "NO");   next if ($chr eq "chrXII" && $start >=451418 && $end <= 469316); my ($type,$details)=(split/\|/,$inf)[0,1];next unless (($type eq "LTR_retrotransposon" && $distdes eq "ENTIRE" && $details =~/Ty1-1/)|($type eq "long terminal repeats")); print "$_\n" }'  > ${SampleID}.Ty1.insert.txt
 
# cat nuc1\ 3\ days.One.txt| perl -ne '{chomp; my ($read,$id,$chr,$start,$end,$strand,$seq,$inf,$ltr,$distdes)=(split/\t/,$_)[1,2,10,11,12,13,40,41,42,43];next if ($inf eq "NO");   next if ($chr eq "chrXII" && $start >=451418 && $end <= 469316); my ($type,$details)=(split/\|/,$inf)[0,1];next unless (($type eq "LTR_retrotransposon" && $distdes eq "ENTIRE" && $details =~/Ty1-1/)|($type eq "long terminal repeats")); print ">$read\n$seq\n" }' >${SampleID}.LTR.fasta

#grep "LTR" ${i}_combined_single_final.Microhomology.txt|perl -ne '{chomp; my @array=split/\t/,$_; if ($array[16]==0){print ">$array[0]\n$array[3]\n"}}' >${i}_combined_single_final.LTR.fasta

# blast these Ty1 insertion againt Ty1 reference, Here we used a blastn-short to allow short alignments
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}.LTR.fasta -db ${Ty1} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${SampleID}.LTR.blast -task blastn-short  -word_size 11 -evalue ${DeEvalue} -dust no -soft_masking false -gapopen ${DepGapsize} -penalty -1 -perc_identity ${DepIden}

# checking the overlapping stats
perl ${srcDir}/Identify_Best_AlignmentsForLTR_Bed.pl -i ${SampleID}.LTR.blast -o ${SampleID}.LTR.uniq -g ${SampleID}

### generate the bed file and for the mapping 

sort -k 1,1V -k2,2n -k 3,3n ${SampleID}.LTR.uniq.bed >${SampleID}.LTR.uniq.sorted.bed

### measure the insertion number
perl ${srcDir}/CalculateCoverage_Retrotransposon.pl -i ${SampleID}.LTR.uniq.sorted.bed -a ${Ty1Ann} -b ${Ty1} -g ${SampleID} -o ${SampleID}.LTR.cov.txt

#### After get all the insertion, we combined them
# perl  ${srcDir}/Combined_Insertion_counts.pl  -i $PWD -o Combined_InsertionCounts.txt
#
# ### Seperate into forward and reverse insertion
# awk '{if ($6=="+"){print}}' ${SampleID}_combined.LTR.bed >${SampleID}_combined_forward.LTR.bed
# awk '{if ($6=="-"){print}}' ${SampleID}_combined.LTR.bed >${SampleID}_combined_reverse.LTR.bed
#
# perl  ${srcDir}/CalculateCoverage_Retrotransposon.pl -i ${SampleID}_combined_reverse.LTR.bed -a ${Ty1Ann} -b ${Ty1} -g ${SampleID}_combined_reverse -o ${OutputDir}/${SampleID}_combined_reverse.LTR.cov.txt
# perl ${srcDir}/CalculateCoverage_Retrotransposon.pl -i ${SampleID}_combined_forward.LTR.bed -a ${Ty1Ann} -b ${Ty1} -g ${SampleID}_combined_forward -o ${OutputDir}/${SampleID}_combined_forward.LTR.cov.txt

echo "Congratulation! Insertion detection and deduplication job array is finished !"
