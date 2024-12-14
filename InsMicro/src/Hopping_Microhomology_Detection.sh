#!/bin/bash
# Sample batchscript to run a simple job on HPC
#SBATCH --partition=mghpcc-compute              # queue to be used
#SBATCH --account=bch-mghpcc                    # account name to be used
#SBATCH --time=12:01:00                         # Running time (in hours-minutes-seconds)
#SBATCH --job-name=WT              # Job name
#SBATCH --output=output_%j.txt                  # Name of the output file
#SBATCH --nodes=1                               # Number of nodes needed
#SBATCH --ntasks=15                             # Number of CPUs needed
#SBATCH --mem=20G                                # Memory needed for the computing session


### Define parameters by hands
source /programs/biogrids.shrc
export JAVA_X=jdk1.8.0_144


SampleID=CAgingV0430

softpath=/scratch/ch220812/Software

## Final combine the samples
#srcDir=

Chr3=${softpath}/iDSBindel/Database/chr3.fasta

Chr3Masked=${softpath}/InsMicro/Database/Chr3.masked.fasta

genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa


cd /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/NC_Kap123_043024/Update_V210924/Hopping_v0519

cp ../*/FinalInsertion/*.txt ./
# generate fasta
cat *.txt|perl -ne '{chomp; my ($name, $id,$sampleID, $string,$identity,$Rcount,$quality)=(split/\t/,$_)[0,1,2,3,7,8,9]; next if (exists $hash{$id} || $id eq "RepresentativeRead" ||$name =~/\./); print ">$id\t$Rcount\t$quality\t$identity\t$sampleID\n$string\n"; $hash{$id}++}' >${SampleID}_Combined.fasta
cat *.txt >${SampleID}_Combined.txt



# Combine the basic information (This part we ignore the rDNA as it would be hard to know whether they are from hopping)


${softpath}/ncbi-blast-2.8.1+/bin/makeblastdb  -in ${SampleID}_Combined.fasta -dbtype nucl
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${SampleID}_Combined.fasta -db ${SampleID}_Combined.fasta -out ${SampleID}_Combined.blast -num_threads 13  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'

#### For the final round deduplication, we briefly adapt self-blast strategy to remove these duplicates. Considering the fact that the longer the reads are, the lower quality the reads will be. MiSeq read length is 300bp, error rate can be as high as 20%. Therefore, we specifically required a loose standard for long insertion deduplicates, the other short still used the same parameters. Regarding the long insertion events, we removed long insertions low-coverage duplicates that are with similar length (difference <=5), high identity (>95%), less mismatches/indels (<= 15% insertion length). 

### 
### please pay attention to this one, even if some reads can be clustered together. 

perl ${softpath}/iDSBins/CombineSamplesrc/SelfBlast_RemoveDuplication_eachsampleIdentitySorted_totalinsertion_finalround_v2_combineSamples_generateRcounts_WT.pl -b ${SampleID}_Combined.blast -i ${SampleID}_Combined.txt -f ${SampleID}_Combined.fasta -o ${SampleID}_Combined


### please pay attention to this one, even if some reads can be clustered together. We only consider the clustering that one with more 2000 another sample less than 5 reads as hopping events.

### here we ranked the read number across samples and check the first >=2000 and second <=5
 perl ${softpath}/iDSBins/CombineSamplesrc/GenerateFinalResultsCombiningSamples_compareReadCovDiff.pl -g ${SampleID}_Combined.cls -i ${SampleID}_Combined.txt -o ${SampleID}_Combined
 

 awk 'NR==1; NR > 1 {print $0 | "sort  -k 11,11V -k 12,12n -k 13,13n"}' ${SampleID}_Combined.One.txt|perl -ne '{chomp; my @array=split/\t/,$_; my ($chr,$start,$end)=@array[10,11,12]; my $string=join "\t",@array[0..20]; my $string2=join "\t",@array[23..24]; my $up; my $down;if (exists $hash{$chr}){$up=$start-$start0; $down=$end-$end0}else{$up=$start; $down=$end} print "$string\t$up\t$down\t$string2\n"; $start0=$start; $end0=$end; $hash{$chr}++;}' >${SampleID}_Combined.One_updated.txt
 awk 'NR==1; NR > 1 {print $0 | "sort  -k 11,11V -k 12,12n -k 13,13n"}' ${SampleID}_Combined.Multiple.txt |perl -ne '{chomp; my @array=split/\t/,$_; my ($chr,$start,$end)=@array[10,11,12];my $string=join "\t",@array[0..20]; my $string2=join "\t",@array[23..24]; if ($start eq "Unknown" || $start  eq "NO"){print "$string\tUnknown\tUnknown\t$string2\n"; next}; my $up; my $down;if (exists $hash{$chr}){$up=$start-$start0; $down=$end-$end0}else{$up=$start; $down=$end} print "$string\t$up\t$down\t$string2\n"; $start0=$start; $end0=$end; $hash{$chr}++;}'  >${SampleID}_Combined.Multiple_updated.txt
 
 #### Then genetate the final result by removin the duplicates that have 0-10 bp shift and low coverage, and add the new read count number
 mkdir -p FinalInsertionNoMicro
 cd FinalInsertionNoMicro
 
 ## remove hopping insertions 
 ### we igonore the multiple donor title of multiple donor:
 
 ## combining the donor for multiple and single
 cat ../${SampleID}_Combined.One_updated.txt  ../${SampleID}_Combined.Multiple_updated.txt | perl -ne '{chomp; my ($chr,$start,$end)=(split/\t/,$_)[10];next if ($chr =~/;/); print "$_\n" }' >${SampleID}_CombinedDonors.txt
 
 
 
 sort -k 11,11V -k 12,12n ${SampleID}_CombinedDonors.txt| perl -ne '{chomp; my ($donor,$start,$end)=(split/\t/,$_)[6,11,12];my $inf=($donor eq "1")?"Single":"Multiple";   my $dif1=$start-$start0; my $dif2=$end-$end0; $start0=$start; $end0=$end;print "$_\t$inf\t$dif1\t$dif2\n"}' >${SampleID}_Checkhopping.txt

perl  ${softpath}/iDSBins/CombineSamplesrc/TestRemoveHopping_norDNALTYMt.pl -i ${SampleID}_Checkhopping.txt -o ${SampleID}_Removehopping >Hoppinginformation.txt


 ### We genrate the final mutliple donor insertion events.

 grep "Single" -v ${SampleID}_Removehopping.Corrected.txt >${SampleID}_final.MultipleIndex.txt
 
 ### We generate the final one insertion event
 grep "Single|CaseID" -E ${SampleID}_Removehopping.Corrected.txt |sort -k 11,11V -k 12,12n |perl -ne '{chomp; my @array=split/\t/,$_; if($array[0] eq "CaseID"){print "$_\n";next;} $n++; my $string=join "\t",("S$n",@array[1..24]); print "$string\n"}' >${SampleID}_final_Removinghopping.One.txt

 ### generate final table with microhomology multiple
 
 perl  ${softpath}/iDSBins/CombineSamplesrc/GnerateMultipleInsertion2Removinghopping.pl -i ${SampleID}_final.MultipleIndex.txt -g ../${SampleID}_Combined.Multiple_updated.txt -o ${SampleID}_final_Removinghopping.Multiple.txt
 
  
 cp ${SampleID}_final_Removinghopping.One.txt ${SampleID}_final.one.txt
 
 cp ${SampleID}_final_Removinghopping.Multiple.txt ${SampleID}_final.multiple.txt
 
 cd ../
 
 mkdir FinalInsertionMicro
 cd FinalInsertionMicro



 ### This combination should be modififed by sample integration as ${SampleID} and ${SampleID} have some insertion hopping events
 # cat ../*.multiple.txt >${SampleID}_all.multiple.txt
 # cat ../*.sinlge.txt >${SampleID}_all.one.txt

 cp ../FinalInsertionNoMicro/${SampleID}_final.one.txt ../FinalInsertionNoMicro/${SampleID}_final.multiple.txt ./

 perl -ne '{chomp; my @array=split/\t/,$_;next if ($array[0] eq "CaseID");my $id=$array[1]; next if ($id eq "");next if (exists $hash{$id}); my $seq= substr ($array[3], 3,-3); print ">$id\n$seq\n"; $hash{$id}++;}'  ${SampleID}_final.one.txt >${SampleID}_final.one.fasta

 # determine the two side locus; DUE TO inaccurate alignment, we performed the blast with -task blastn-short
 # Here we only retain the MAT sequence of Chr3

 #/archive/tmhxxw9/software/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in Chr3.masked.fasta
 ### two round of blast results, some locus cannot sensentively to map against the genome

 ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query   ${SampleID}_final.one.fasta -db ${Chr3Masked} -num_threads 8  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${SampleID}_final.one.twoside.noshortblast.blast

 ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query   ${SampleID}_final.one.fasta -db ${Chr3Masked} -num_threads 8 -task blastn-short  -word_size 5 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${SampleID}_final.one.twoside.short.blast

 cat ${SampleID}_final.one.twoside.noshortblast.blast ${SampleID}_final.one.twoside.short.blast >${SampleID}_final.one.twoside.finalblast.blast

 ### determine the insertion events

 ### two round of blast results

 # We first determine the best blast igonore the short ones
 ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${SampleID}_final.one.fasta  -db  ${genomeseq} -num_threads 15 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${SampleID}_final.one.insertion.noshortblast.blast

 ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${SampleID}_final.one.fasta  -db  ${genomeseq} -num_threads 15 -task blastn-short  -word_size 11 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'  -out ${SampleID}_final.one.insertion.shortblast.blast


 cat ${SampleID}_final.one.insertion.noshortblast.blast ${SampleID}_final.one.insertion.shortblast.blast >${SampleID}_final.one.insertion.finalblast.blast

 ### then detect the microhomoloy and short indels at junction 

 perl ${softpath}/InsMicro/src/Microhomology_Detection_Final_multiple_v2.pl -i  ${SampleID}_final.one.fasta -g ${genomeseq} -b ${SampleID}_final.one.insertion.finalblast.blast -t ${SampleID}_final.one.twoside.finalblast.blast -s ${SampleID}_final.one.txt -o ${SampleID}_final.one.Micro.txt >${SampleID}_single_inf.txt



 ### For multiple donor
 ### two round of blast results, some locus cannot sensentively to map against the genome

 perl -ne '{chomp; my @array=split/\t/,$_;next if ($array[0] eq "CaseID");my $id=$array[1]; next if (exists $hash{$id}); my $seq= substr ($array[3], 3,-3); print ">$id\n$seq\n"; $hash{$id}++;}'  ${SampleID}_final.multiple.txt >${SampleID}_final.multiple.fasta



 ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${SampleID}_final.multiple.fasta -db ${Chr3Masked} -num_threads 8  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${SampleID}_final.multiple.twoside.noshortblast.blast

 ### because it is very short, we reduced the blastn size, this will be very sensitive to the large deletion of MAT 
 ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query   ${SampleID}_final.multiple.fasta -db  ${Chr3Masked} -num_threads 8 -task blastn-short -word_size 5 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out  ${SampleID}_final.multiple.twoside.short.blast

 cat ${SampleID}_final.multiple.twoside.noshortblast.blast ${SampleID}_final.multiple.twoside.short.blast >${SampleID}_final.multiple.twoside.finalblast.blast

 ### determine the insertion events

 ### two round of blast results
 ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_final.multiple.fasta -db  ${genomeseq} -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${SampleID}_final.multiple.insertion.noshortblast.blast

 ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_final.multiple.fasta -db  ${genomeseq} -num_threads 8 -task blastn-short -word_size 11 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${SampleID}_final.multiple.insertion.shortblast.blast

 cat ${SampleID}_final.multiple.insertion.noshortblast.blast ${SampleID}_final.multiple.insertion.shortblast.blast >${SampleID}_final.multiple.insertion.finalblast.blast


 perl ${softpath}/InsMicro/src/Microhomology_Detection_Final_multiple_v2.pl -i ${SampleID}_final.multiple.fasta -g ${genomeseq} -b ${SampleID}_final.multiple.insertion.finalblast.blast -t ${SampleID}_final.multiple.twoside.finalblast.blast -s ${SampleID}_final.multiple.txt -o ${SampleID}_final.multiple.Micro.txt >${SampleID}_multiple_inf.txt


 ## update the insertion

 perl ${softpath}/InsMicro/src/GnerateMultipleInsertion3.pl  -i ${SampleID}_final.multiple.Micro.txt -g ${SampleID}_final.multiple.txt -o ${SampleID}_final_updated.Multiple.Micro.txt


sort -k 11,11V -k 12,12n ${SampleID}_final.one.Micro.txt|perl -ne '{chomp; my @array=split/\t/,$_; if($array[0] eq "CaseID"){print "$_\n"; next} $num++; shift @array; my $string=join "\t",@array; print "S$num\t$string\n"}' > ${SampleID}_final_updated.One.Micro.txt


  mkdir FinalMicro
  cd FinalMicro
  cp ../${SampleID}_final_updated.* ./
 
 
 
 