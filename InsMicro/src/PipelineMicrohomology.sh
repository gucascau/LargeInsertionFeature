#!/bin/bash
# Sample batchscript to run a simple job array to create 5 different files, and modify the files independently with one batchscript
 
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=1:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=yYY521_indel # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=xin.wang@childrens.harvard.edu # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=20 # Number of cpu cores on one node
#SBATCH --mem=20G



# Get absolute path for scripts and check if required scripts exist
source /programs/biogrids.shrc

SampleID=yXC876A_S5

in=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Yang-Yu_miseqyeast_03132020/Alignment_1/20200313_214022/Fastq/Update_V210924


FORWARD=${SampleID}_QC.highquality.R1.fastq
REVERSE=${SampleID}_QC.highquality.R2.fastq
INSERT=${SampleID}_detected_combined.highqual.txt

### Set default options
# some general paramaters including mata locus thread, and software path
# default threads
nproc=15
# default software path and Chromosome 3 index

softpath=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/Updated_V0903/UsedMutants/TestSoftWare/Software

Chr3=${softpath}/iDSBindel/Database/chr3.fasta

Chr3Masked=${softpath}/InsMicro/Database/Chr3.masked.fasta

genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

# default the customed index size
Findexsize=3
Rindexsize=3

# default whole MAT length
Matasize=84
# default Mata chromosome 
Matachr="ChrIII"
# default Mata start site
Matastart=294300
# default Mata end site
Mataend=294500

# default minimum insertion length
Ilength=10

# default parameter for first cuting-edge deduplication
### the size of primer (We allowed 5bp frameshift)
LeftPrimerSize=25
RightPrimerSize=22


### the size of left and right MAT size (we also allowed 5bp frameshift)
LeftMataSize=45
RightMataSize=51

### The read count to define the unique events
ReadCount=5

### The collect upstream and downstream for the unique of deletion or insertions from the raw reads (5bp)
Cutstream=5

# Begin script in case all parameters are correct
echo "${SampleID}"
echo "${in}"
echo "${INSERT}"
echo "${FORWAR}"
echo "${REVERSE}"
echo "${softpath}"

# 

echo "All the paramter are sucessfully provided, now let's detect the large insertion event"


#### Here is the last version of large insertion events


echo "The job array of microhomology detection is started ..."
echo ""
date




### Trim the upstring and downstream 3bp


mkdir Microhomology
cd Microhomology


Type=WT

### This combination should be modififed by sample integration as ${Type} and DNA2 have some insertion hopping events
# cat ../*.multiple.txt >${Type}_all.multiple.txt
# cat ../*.sinlge.txt >${Type}_all.one.txt

cp ../FinalSampleIntegration/FinalInsertion/${Type}_final.* ./

perl -ne '{chomp; my @array=split/\t/,$_;next if ($array[0] eq "CaseID");my $id=$array[1]; next if ($id eq "");next if (exists $hash{$id}); my $seq= substr ($array[3], 3,-3); print ">$id\n$seq\n"; $hash{$id}++;}'  ${Type}_final.one.txt >${Type}_final.one.fasta

# determine the two side locus; DUE TO inaccurate alignment, we performed the blast with -task blastn-short
# Here we only retain the MAT sequence of Chr3

#/archive/tmhxxw9/software/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in Chr3.masked.fasta
### two round of blast results, some locus cannot sensentively to map against the genome

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query   ${Type}_final.one.fasta -db ${Chr3Masked} -num_threads 8  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${Type}_final.one.twoside.noshortblast.blast

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query   ${Type}_final.one.fasta -db ${Chr3Masked} -num_threads 8 -task blastn-short  -word_size 5 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${Type}_final.one.twoside.short.blast

cat ${Type}_final.one.twoside.noshortblast.blast ${Type}_final.one.twoside.short.blast >${Type}_final.one.twoside.finalblast.blast

### determine the insertion events

### two round of blast results

# We first determine the best blast igonore the short ones
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${Type}_final.one.fasta  -db  ${genomeseq} -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${Type}_final.one.insertion.noshortblast.blast

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${Type}_final.one.fasta  -db  ${genomeseq} -num_threads 8 -task blastn-short  -word_size 11 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'  -out ${Type}_final.one.insertion.shortblast.blast


cat ${Type}_final.one.insertion.noshortblast.blast ${Type}_final.one.insertion.shortblast.blast >${Type}_final.one.insertion.finalblast.blast

### then detect the microhomoloy and short indels at junction 

perl ${softpath}/InsMicro/src/Microhomology_Detection_Final_multiple_v2.pl -i  ${Type}_final.one.fasta -g ${genomeseq} -b ${Type}_final.one.insertion.finalblast.blast -t ${Type}_final.one.twoside.finalblast.blast -s ${Type}_final.one.txt -o ${Type}_final.one.Micro.txt >${Type}_single_inf.txt



### For multiple donor
### two round of blast results, some locus cannot sensentively to map against the genome

perl -ne '{chomp; my @array=split/\t/,$_;next if ($array[0] eq "CaseID");my $id=$array[1]; next if (exists $hash{$id}); my $seq= substr ($array[3], 3,-3); print ">$id\n$seq\n"; $hash{$id}++;}'  ${Type}_final.multiple.txt >${Type}_final.multiple.fasta



${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${Type}_final.multiple.fasta -db ${Chr3Masked} -num_threads 8  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${Type}_final.multiple.twoside.noshortblast.blast

### because it is very short, we reduced the blastn size, this will be very sensitive to the large deletion of MAT 
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query   ${Type}_final.multiple.fasta -db  ${Chr3Masked} -num_threads 8 -task blastn-short -word_size 5 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out  ${Type}_final.multiple.twoside.short.blast

cat ${Type}_final.multiple.twoside.noshortblast.blast ${Type}_final.multiple.twoside.short.blast >${Type}_final.multiple.twoside.finalblast.blast

### determine the insertion events

### two round of blast results
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${Type}_final.multiple.fasta -db  ${genomeseq} -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${Type}_final.multiple.insertion.noshortblast.blast

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${Type}_final.multiple.fasta -db  ${genomeseq} -num_threads 8 -task blastn-short -word_size 11 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${Type}_final.multiple.insertion.shortblast.blast

cat ${Type}_final.multiple.insertion.noshortblast.blast ${Type}_final.multiple.insertion.shortblast.blast >${Type}_final.multiple.insertion.finalblast.blast


perl ${softpath}/InsMicro/src/Microhomology_Detection_Final_multiple_v2.pl -i ${Type}_final.multiple.fasta -g ${genomeseq} -b ${Type}_final.multiple.insertion.finalblast.blast -t ${Type}_final.multiple.twoside.finalblast.blast -s ${Type}_final.multiple.txt -o ${Type}_final.multiple.Micro.txt >${Type}_multiple_inf.txt


perl -ne '{chomp; my @array=split/\t/,$_; next if ($array[41] eq "Multiple" && ($array[14] =~/^gene\|/) && (($array[42] <=5 && $array[42]>=0)||($array[43]>=0 && $array[43]<=5))); print "$_\n"}' DNA2_final.Combined.micro_v3.txt >DNA2_final.Combined.micro_v3_nofragement.txt
