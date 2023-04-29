### Ty1 analyses ####
#
# The request 
### It is very interesting that the pattern of Ty insertion looks different in dna2, many other mutants and aging. This indicates the mechanism could be different to generate the substrates from Ty. It would ### be good if you generate some figures of Ty like what you did for rDNA. It will not show every events using bars but show the percent time nucleotide among all insertions.  

 

### The following are some thoughts to do this analysis.

### TWe only take Ty1-1 from Annotation (Retrotransposon) in the table because it is the most dominant Ty and easier for mapping.
###  We will only analyze the Ty from single donor events because it is more accurate for the edge than multiple donor.
### If there is mismatch between the query sequence and subject sequence during blast, we use the whole subject sequence to do mapping.
### If they were mapped to both LTR, use the right LTR because it was reported that cDNA is generated from the right one. However, there are some extra nucleotides (1-3) which are perfectly match the left one. ### However, it could be misleading to be mapped to the left one considering the mechanisms how cDNA is generated from Ty.
### We need chronological aging, replicative aging and all these mutants. pif1-m2 dna2, sgs1 exo1, rtt105, rnh and derivatives, sgs1, rad51, rad52, mec1 tel1 sml1, rad5 and rad27. The criteria to select these ### mutants is the frequency of Ty insertion is increased 5 fold more than WT and the Ty insertion number is more than 15 in single donor events. The exception is rad27. In rad27, the frequency of Ty insertion ### is only 1.4 fold of WT. However, we have 18 events. We may not put in the paper. Let's see if there is something interesting for rad27.
### For chronological aging, we need the samples from wt, nuc1, atg1 and spt3. You could separate the 1d, 2d, 3d, 8d and 16d. Maybe there will be some difference at different time point. You could also combine 
### all time points for each genotype.


# 
### LTR analaysis

#source /programs/biogrids.shrc
#/home/ch220812/software/ncbi-blast-2.8.1+/bin/makeblastdb -in Retrotransposon_Seq_Yang.fasta -dbtype nucl

# ~/Software/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in YPLWTY1-1.fasta

### change name ID
###
###     for f in *' '*; do mv "$f" "`echo $f | sed -e 's/ /_/g'`"; done

Fid=

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

## test for nuc1 

# checking the overlapping stats
perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Scripts/Identify_Best_AlignmentsForLTR_Bed.pl -i ${Fid}.LTR.blast -o ${Fid}.LTR.uniq -g ${Fid}

### generate the bed file and for the mapping 

sort -k 1,1V -k2,2n -k 3,3n ${Fid}.LTR.uniq.bed >${Fid}.LTR.uniq.sorted.bed

### measure the insertion number
perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Scripts/CalculateCoverage_Retrotransposon.pl -i ${Fid}.LTR.uniq.sorted.bed -a ${Ty1Ann} -b ${Ty1} -g ${Fid} -o ${Fid}.LTR.cov.txt




#### After get all the insertion, we combined them
perl  /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Scripts/Combined_Insertion_counts.pl  -i $PWD -o Combined_InsertionCounts.txt



Fid=wt

### Seperate into forward and reverse insertion
awk '{if ($6=="+"){print}}' ${Fid}_combined.LTR.bed >${Fid}_combined_forward.LTR.bed
awk '{if ($6=="-"){print}}' ${Fid}_combined.LTR.bed >${Fid}_combined_reverse.LTR.bed

perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Scripts/CalculateCoverage_Retrotransposon.pl -i ${Fid}_combined_reverse.LTR.bed -a ${Ty1Ann} -b ${Ty1} -g ${Fid}_combined_reverse -o ${Fid}_combined_reverse.LTR.cov.txt
perl /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/Ty1Analyses/Scripts/CalculateCoverage_Retrotransposon.pl -i ${Fid}_combined_forward.LTR.bed -a ${Ty1Ann} -b ${Ty1} -g ${Fid}_combined_forward -o ${Fid}_combined_forward.LTR.cov.txt

###
# perl -ne '{chomp; my ($matches,$ilength,$start,$end,$length)=(split/\t/,$_)[3,8,9,10,11]; my $cov=($matches/$ilength);my $start0=($start<$end)?$start:$end; my $end0=($start<$end)?$end:$start; my $length0=$length-$start0; print "$_\n" if ($cov>0.8 && ($end0<=341 || $length0 <=341))}'


### idenfy the best blast hit choose the top1 and transfer into bedfile
# $perl -ne '{chomp; my @array=split/\t/,$_; my $id =$array[0];  next if(exists $hash{$id}); my ($start,$end,$strand)=($array[9]>$array[10])? ($array[10],$array[9],"-"):($array[9],$array[10],"+"); my $string=join "\t",($array[1],$start,$end,$id,"0",$strand); print "$string\n"; $hash{$id}++;}' ${i}_combined_single_final.LTR.blast |sort -k 1,1 -k 2,2n >${i}_combined_single_final.LTR.besthit.bed


# cat ${i}_combined_single_final.LTR.blast ${i}_combined_mutiple_final.LTR.blast |perl -ne '{chomp; my @array=split/\t/,$_; my $id =$array[0];  next if(exists $hash{$id}); my ($start,$end,$strand)=($array[9]>$array[10])? ($array[10],$array[9],"-"):($array[9],$array[10],"+"); my $string=join "\t",($array[1],$start,$end,$id,"0",$strand); print "$string\n"; $hash{$id}++;}'  |sort -k 1,1 -k 2,2n >${i}_combined_all_final.LTR.besthit.bed


#### Generate the site of end insertion from bed file:
cat ../wt_*_day* |sort -k 1,1V -k 2,2n -k 3,3n >wt_chronologicalaging.bed
cat ../nuc1_3_days.LTR.uniq.bed ../nuc1_8_days.LTR.uniq.bed ../nuc1_16_days.LTR.uniq.bed |sort -k 1,1V -k 2,2n -k 3,3n >nuc1_chronologicalaging.bed
cat ../spt3_*_days*.bed|sort -k 1,1V -k 2,2n -k 3,3n >spt3_chronologicalaging.bed
cp ../atg1_8_days.LTR.uniq.bed ./
mv atg1_8_days.LTR.uniq.bed atg1_chronologicalaging.bed
cp ../pif1-m2_dna2.LTR.uniq.bed ./
cp ../rtt105.LTR.uniq.bed ./

#### add the color information

perl -ne '{chomp;$n++; print "track itemRgb=On\n" if ($n==1); my ($chr,$start,$end,$id,$zero,$strand,$inf)=split/\t/,$_; my $color= ($strand eq "+")?"255,0,0":"0,0,255"; print "$chr\t$start\t$end\t$id\t$zero\t$strand\t$start\t$end\t\"$color\"\n"}' atg1_chronologicalaging.bed|head >atg1_chronologicalaging_color.bed



### Using UCSC 

 Chromosome XVI 56452..62375
 
 chrVII	56452	
 
 
perl -ne '{chomp; $n++; if ($n==1){print "track name=\"WT\" description=\"Item RGB demonstration\" visibility=2 itemRgb=\"On\"\n"; next }; my ($chr,$start,$end,$id,$zero,$strand,$start,$end,$color)=split/\t/,$_; my $fstart=$start+56451;my $fend=$end+56451;$color=~s/\"//g; print "chrXVI\t$fstart\t$fend\t$id\t$zero\t$strand\t$fstart\t$fend\t$color\n"}' wt_chronologicalaging_color.bed >wt_chronologicalaging_color_chromosme.bed
 
 
perl -ne '{chomp; $n++; if ($n==1){print "track name=\"nuc1\" description=\"Item RGB demonstration\" visibility=2 itemRgb=\"On\"\n"; next }; my ($chr,$start,$end,$id,$zero,$strand,$start,$end,$color)=split/\t/,$_; my $fstart=$start+56451;my $fend=$end+56451;$color=~s/\"//g; print "chrXVI\t$fstart\t$fend\t$id\t$zero\t$strand\t$fstart\t$fend\t$color\n"}' nuc1_chronologicalaging_color.bed >nuc1_chronologicalaging_color_chromosome.bed

 
perl -ne '{chomp; $n++; if ($n==1){print "track name=\"atg1\" description=\"Item RGB demonstration\" visibility=2 itemRgb=\"On\"\n"; next }; my ($chr,$start,$end,$id,$zero,$strand,$start,$end,$color)=split/\t/,$_; my $fstart=$start+56451;my $fend=$end+56451;$color=~s/\"//g; print "chrXVI\t$fstart\t$fend\t$id\t$zero\t$strand\t$fstart\t$fend\t$color\n"}' atg1_chronologicalaging_color.bed >atg1_chronologicalaging_color_chromosome.bed

perl -ne '{chomp; $n++; if ($n==1){print "track name=\"pif1-m2 dna2\" description=\"Item RGB demonstration\" visibility=2 itemRgb=\"On\"\n"; next }; my ($chr,$start,$end,$id,$zero,$strand,$start,$end,$color)=split/\t/,$_; my $fstart=$start+56451;my $fend=$end+56451;$color=~s/\"//g; print "chrXVI\t$fstart\t$fend\t$id\t$zero\t$strand\t$fstart\t$fend\t$color\n"}' pif1-m2_dna2.LTR.uniq_color.bed >pif1-m2_dna2.LTR.uniq_color_chromosome.bed


perl -ne '{chomp; $n++; if ($n==1){print "track name=\"rtt105\" description=\"Item RGB demonstration\" visibility=2 itemRgb=\"On\"\n"; next }; my ($chr,$start,$end,$id,$zero,$strand,$start,$end,$color)=split/\t/,$_; my $fstart=$start+56451;my $fend=$end+56451;$color=~s/\"//g; print "chrXVI\t$fstart\t$fend\t$id\t$zero\t$strand\t$fstart\t$fend\t$color\n"}' rtt105.LTR.uniq_color.bed >rtt105.LTR.uniq_color_chromosome.bed



perl -ne '{chomp; $n++; if ($n==1){print "track name=\"spt3\" description=\"Item RGB demonstration\" visibility=2 itemRgb=\"On\"\n"; next }; my ($chr,$start,$end,$id,$zero,$strand,$start,$end,$color)=split/\t/,$_; my $fstart=$start+56451;my $fend=$end+56451;$color=~s/\"//g; print "chrXVI\t$fstart\t$fend\t$id\t$zero\t$strand\t$fstart\t$fend\t$color\n"}' spt3_chronologicalaging_color.bed >spt3_chronologicalaging_color_chromosome.bed



### counting the polarity of the inserts

LTR region: 1-337/5588-5924

perl -ne '{chomp; my ($start,$end,$strand)=(split/\t/,$_)[1,2,5]; next unless ($start>=5588 && $end<=5924); $num++ if ($strand eq "+") }END {print "$num"}' wt_chronologicalaging.bed
perl -ne '{chomp; my ($start,$end,$strand)=(split/\t/,$_)[1,2,5]; next unless ($start>=5588 && $end<=5924); $num++ if ($strand eq "-") }END {print "$num"}' wt_chronologicalaging.bed
perl -ne '{chomp; my ($start,$end,$strand)=(split/\t/,$_)[1,2,5]; next if ($start>=5588 && $end<=5924); $num++ if ($strand eq "-") }END {print "$num"}' wt_chronologicalaging.bed
perl -ne '{chomp; my ($start,$end,$strand)=(split/\t/,$_)[1,2,5]; next if ($start>=5588 && $end<=5924); $num++ if ($strand eq "+") }END {print "$num"}' wt_chronologicalaging.bed

### measure the number of AACA, ACA, CA, AC, AA
# Top strand
perl -ne '{chomp; my $str=$_; my $count= ($str =~ s/CA/CA/g); print "$count\n"}' YPLWTY1_sequence_Top.txt
# bottom strand
perl -ne '{chomp; my $str0=$_; my $str=reverse $str0; $str=~  tr/ATGCatgc/TACGtacg/; my $count= ($str =~ s/AC/AC/g); print "$str\n" }' YPLWTY1_sequence_Top.txt