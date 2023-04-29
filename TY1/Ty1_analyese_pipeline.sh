---
title: "Coder Upgrade 2023: Ty1 nucleotide insertion map analysis"
author: "Xin Wang"
date: "`r format(Sys.time(), '%m/%d/%Y')`"
output:
    html_document:
        toc: true
        toc_float: true
        toc_depth: 5
        number_sections: false
        code_folding: show
---

``` {sh}
# Set the start time
message("Starting serial processing for Ty1 nucleotide insertion map:\n")
t <- c(Sys.time())

#!/bin/bash

### Workflow for Ty1 analyses ####

### It is very interesting that the pattern of Ty insertion looks different in dna2, many other mutants and aging. This indicates the mechanism could be different to generate the substrates from Ty. It would ### be good if you generate some figures of Ty like what you did for rDNA. It will not show every events using bars but show the percent time nucleotide among all insertions.  

### The following are some thoughts to do this analysis.

### We only take Ty1-1 from Annotation (Retrotransposon) in the table because it is the most dominant Ty and easier for mapping.
### We will only analyze the Ty from single donor events because it is more accurate for the edge than multiple donor.
### If there is mismatch between the query sequence and subject sequence during blast, we use the whole subject sequence to do mapping.
### If they were mapped to both LTR, use the right LTR because it was reported that cDNA is generated from the right one. However, there are some extra nucleotides (1-3) which are perfectly match the left one. ### However, it could be misleading to be mapped to the left one considering the mechanisms how cDNA is generated from Ty.
### We need chronological aging, replicative aging and all these mutants. pif1-m2 dna2, sgs1 exo1, rtt105, rnh and derivatives, sgs1, rad51, rad52, mec1 tel1 sml1, rad5 and rad27. The criteria to select these ### mutants is the frequency of Ty insertion is increased 5 fold more than WT and the Ty insertion number is more than 15 in single donor events. The exception is rad27. In rad27, the frequency of Ty insertion ### is only 1.4 fold of WT. However, we have 18 events. We may not put in the paper. Let's see if there is something interesting for rad27.
### For chronological aging, we need the samples from wt, nuc1, atg1 and spt3. You could separate the 1d, 2d, 3d, 8d and 16d. Maybe there will be some difference at different time point. You could also combine 
### all time points for each genotype.

### LTR analaysis

#source /programs/biogrids.shrc
#/home/ch220812/software/ncbi-blast-2.8.1+/bin/makeblastdb -in Retrotransposon_Seq_Yang.fasta -dbtype nucl

# ~/Software/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in YPLWTY1-1.fasta

### change name ID
###
### for f in *' '*; do mv "$f" "`echo $f | sed -e 's/ /_/g'`"; done

Fid= Name

###single insertion:

## updated the name
#cp wt 3 days.One.txt ${Fid}.One.txt

Ty1=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Ty1Sequence/YPLWTY1-1.fasta
Ty1Ann=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Ty1Sequence/YPLWTY1-1.annotation.bed

cat ${Fid}.One.txt| perl -ne '{chomp; my ($read,$id,$chr,$start,$end,$strand,$seq,$inf,$ltr,$distdes)=(split/\t/,$_)[1,2,10,11,12,13,40,41,42,43];next if ($inf eq "NO");   next if ($chr eq "chrXII" && $start >=451418 && $end <= 469316); my ($type,$details)=(split/\|/,$inf)[0,1];next unless (($type eq "LTR_retrotransposon" && $distdes eq "ENTIRE" && $details =~/Ty1-1/)); print ">$read\n$seq\n" }'  > ${Fid}.LTR.fasta
 
# print the insertion events:

cat ${Fid}.One.txt| perl -ne '{chomp; my ($read,$id,$chr,$start,$end,$strand,$seq,$inf,$ltr,$distdes)=(split/\t/,$_)[1,2,10,11,12,13,40,41,42,43];next if ($inf eq "NO");   next if ($chr eq "chrXII" && $start >=451418 && $end <= 469316); my ($type,$details)=(split/\|/,$inf)[0,1];next unless (($type eq "LTR_retrotransposon" && $distdes eq "ENTIRE" && $details =~/Ty1-1/)); print "$_\n" }'  > ${Fid}.Ty1.insert.txt
 
# cat nuc1\ 3\ days.One.txt| perl -ne '{chomp; my ($read,$id,$chr,$start,$end,$strand,$seq,$inf,$ltr,$distdes)=(split/\t/,$_)[1,2,10,11,12,13,40,41,42,43];next if ($inf eq "NO");   next if ($chr eq "chrXII" && $start >=451418 && $end <= 469316); my ($type,$details)=(split/\|/,$inf)[0,1];next unless (($type eq "LTR_retrotransposon" && $distdes eq "ENTIRE" && $details =~/Ty1-1/)); print ">$read\n$seq\n" }' >${Fid}.LTR.fasta

#grep "LTR" ${i}_combined_single_final.Microhomology.txt|perl -ne '{chomp; my @array=split/\t/,$_; if ($array[16]==0){print ">$array[0]\n$array[3]\n"}}' >${i}_combined_single_final.LTR.fasta

~/Software/ncbi-blast-2.8.1+/bin/blastn -query ${Fid}.LTR.fasta -db ${Ty1} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${Fid}.LTR.blast -task blastn-short  -word_size 11 -evalue 0.01 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 60

# checking the overlapping stats
perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Scripts/Identify_Best_AlignmentsForLTR_Bed.pl -i ${Fid}.LTR.blast -o ${Fid}.LTR.uniq -g ${Fid}

### generate the bed file and for the mapping 

sort -k 1,1V -k2,2n -k 3,3n ${Fid}.LTR.uniq.bed >${Fid}.LTR.uniq.sorted.bed

### measure the insertion number
perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Scripts/CalculateCoverage_Retrotransposon.pl -i ${Fid}.LTR.uniq.sorted.bed -a ${Ty1Ann} -b ${Ty1} -g ${Fid} -o ${Fid}.LTR.cov.txt

#### After get all the insertion, we combined them
perl  /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Scripts/Combined_Insertion_counts.pl  -i $PWD -o Combined_InsertionCounts.txt

### Seperate into forward and reverse insertion
awk '{if ($6=="+"){print}}' ${Fid}_combined.LTR.bed >${Fid}_combined_forward.LTR.bed
awk '{if ($6=="-"){print}}' ${Fid}_combined.LTR.bed >${Fid}_combined_reverse.LTR.bed

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Scripts/CalculateCoverage_Retrotransposon.pl -i ${Fid}_combined_reverse.LTR.bed -a ${Ty1Ann} -b ${Ty1} -g ${Fid}_combined_reverse -o ${Fid}_combined_reverse.LTR.cov.txt
perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Scripts/CalculateCoverage_Retrotransposon.pl -i ${Fid}_combined_forward.LTR.bed -a ${Ty1Ann} -b ${Ty1} -g ${Fid}_combined_forward -o ${Fid}_combined_forward.LTR.cov.txt

message("Congratulations! Ty1 nucleotide insertion map is completed:\n")
t <- c(Sys.time())

```