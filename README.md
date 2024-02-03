# Large Insertion Feature Analyses

## Feature analyses of large insertion events

Insertions of mobile elements, mitochondrial DNA and fragments of nuclear chromosomes at DNA double-strand breaks sites (DSBs) threaten genome integrity and are common in cancer. Despite extensive efforts, our knowledge of these insertions still remains unknown. These tutorials will focus on the different features of large insertion events.

##  Availability 

  1. Ty1 realignments, Ty nucleotide insertion map
  2. Mitochondrial insertion distribution
  3. Mata insertion detection
  4. Inverted insertion detection
  5. ssDNA transformation
  6. Genomic preference: ARS, Telomere, R-loop, tRNA, tandem repeats, essential genes
  7. Predictions of non-B structures of insertions: 

## Dependencies

Perl is used to run the scripts. The following softwares are also required:

. blast (ncbi-blast-2.8.1+)(https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

. bedtools (bedtools-2.25.0) (https://bedtools.readthedocs.io/en/latest/)

## Install

```
    cd ~
    git clone https://github.com/gucascau/LargeInsertionFeature
```   

## Usage

### Large Insertion Packages:
1. Ty1Nuceotide

1.1 Introduction of Ty1Nuceotide
In some mutants and aging cells, we detected tons of insertions from retrotransposons. However, the retrotransposons were widely spread across the whole genome, which will result in the challenging to evaluate the retrotransposon insertion mechanism. To investigate the Ty1 nucleotide insertion frequency in each mutant, we have developed this package to well address the insertion distribution across Ty1.
 
1.2 Usage: 
```
Usage: sh Ty1_coverage.sh -a Sample ID -b Work Directory  -f Insertion events -r Output  -p Software installed Directory
Note: In some mutants and aging cells, we detected tons of insertions from retrotransposons. However, the retrotransposons were widely spread across the whole genome, which will result in the challenging to evaluate the retrotransposon insertion mechanism. To investigate the Ty1 nucleotide insertion frequency in each mutant, we have developed this package to well address the insertion nucleotide distribution across Ty1.

Request Parameters:
	-a Sample Id (Example: yYY398-B_S10)
	-b The working directory, where you put the insertion events
	-f The insertion event (Note: We only considered the insertion defined as single donor as the insertions with multiple donors might not be able to have a clean Ty1 locus.)
	-p Software installed Path (Default:, it is the path where you store your blast software)


Optional Parameters:

Optional Parameters -- Ty1 Blast setting up:
	-ms identity (Default: 60)
	-mc Gapsize  (Default t)
	-mb E value  (Default t)

Optional Parameters of Ty1 reference
	-c Ty1 blast index (Note: if you don't have blast index, please run the 'makeblastdb -in Ty1.fasta -dbtype nucl') (Default: ./LargeInsertionFeature/Ty1Nuceotide/Database/YPLWTY1-1.fasta)
	-d The Ty1 annotation file (Default: ./LargeInsertionFeature/Ty1Nuceotide/Database/YPLWTY1-1.annotation.bed)
	-r Script Stored Path (Default: ./LargeInsertionFeature/Ty1Nuceotide/src)

	-h help

```
1.3 Output
The output will contain six columns, including the mapped Ty1 position and the nucleotide number of insertions.

| column | explaination |
| ------| ------|
| 1st | Ty1 ID |
| 2nd | Ty1 Position |
| 3rd | Nuclear Sequence |
| 4th | Rank Number |
| 5th | The nucleotide number of Ty1 insertions |
| 6th | Ty1 Annotation |


# Contact


For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu

This project is licensed under the terms of the MIT license.

Copyright (c) 2023 Dr. Kaifu Chen lab

Current version v1.0

