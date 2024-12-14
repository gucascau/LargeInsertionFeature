# Large Insertion Feature Analyses

## Feature analyses of large insertion events

Insertions of mobile elements, mitochondrial DNA and fragments of nuclear chromosomes at DNA double-strand breaks sites (DSBs) threaten genome integrity and are common in cancer. Despite extensive efforts, our knowledge of these insertions still remains unknown. These tutorials will focus on the different features of large insertion events.

##  Availability 

  1. Ty1 realignments, Ty nucleotide insertion map
  2. Genomic preference: ARS, Telomere, R-loop, tRNA, tandem repeats, essential genes

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
####  Ty1Nucleotide

1.1 Introduction of Ty1Nucleotide

In certain mutants and aging cells, we have identified a substantial number of insertions originating from retrotransposons. However, these retrotransposons were found to be extensively distributed throughout the entire genome, posing a challenge in assessing the mechanism of retrotransposon insertion. To delve into the Ty1 nucleotide insertion frequency within each mutant, we have created this package to comprehensively analyze the distribution of insertions across the Ty1 element. This tool aims to provide a detailed understanding of the Ty1 insertion patterns, enabling a more precise investigation into the mechanisms underlying retrotransposon insertions in various mutants and aging cell contexts.
 
1.2 Usage: 

```
Usage: sh Ty1_coverage.sh -a Sample ID -b Work Directory  -f Insertion events -r Output  -p Software installed Directory

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
	-c Ty1 blast index (Note: if you don't have blast index, please run the 'makeblastdb -in Ty1.fasta -dbtype nucl') (Default: ./LargeInsertionFeature/Ty1Nucleotide/Database/YPLWTY1-1.fasta)
	-d The Ty1 annotation file (Default: ./LargeInsertionFeature/Ty1Nucleotide/Database/YPLWTY1-1.annotation.bed)
	-r Script Stored Path (Default: ./LargeInsertionFeature/Ty1Nucleotide/src)

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

####  InsFeatures
2.1 Introduction of InsFeatures

In certain mutants and aging cells, we have identified a substantial number of insertions. Inserted fragments of nuclear genome often originated from fragile regions of the genome such as telomeric, R-loops, or repetitive regions. In DNA2 defeciencty, insertions from nuclear genome tends to be preference for fragile regions such as origins of replication, R-loops, centromeres, telomeres or replication fork barriers (Yang. et.al. Nature. 2018). However, the approximity and ovelapping features of these insertions are still unknown in other mutants and aging cells. Here we introduced the InsFeatures to automatically detect the approximity and overlapping features, including ARS, Telomere, R-loop, Centromere, Tandem repeats, tRNA.

2.2 Usage: 

Please following the instruction of Insertion_Approximity_Overlap_Features.sh and pay attention to the note details.

````
sh Insertion_Approximity_Overlap_Features.sh
````

2.3 Output:

The output will have the feature comparisons between randomed and observed, including ARS, R-loop, centromere, telomere, tRNA, and TRF.

# Citation
_Yu, Y.#, Wang, X.#, Fox, J.#, Yu, R., Thakre, P., Nikoloutsos, N., Yu, Y., Li, Q., Hastings P.J., Dang, W., Chen, K*, Ira, G*_. **Yeast EndoG prevents genome instability by degrading extranuclear DNA species**. _Nat Commun_ 15, 7653 (2024). https://doi.org/10.1038/s41467-024-52147-2

_Yu, Y.#*, Wang, X.#, Fox, J., Li, Q., Yu, Y., Hastings, P.J., Chen, K.* and Ira, G*._ **RPA and Rad27 limit templated and inverted insertions at DNA breaks**. _Nucleic Acids Research_(2024).  https://doi.org/10.1093/nar/gkae1159
# Contact

For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu

This project is licensed under the terms of the MIT license.

Copyright (c) 2023 Dr. Kaifu Chen lab

Current version v1.0

