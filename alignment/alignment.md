# Sequence Alignment, Preprocessing, and Quality Check
Author: Jaehee Kim

## Data
COVID-19 data were downloaded from **GISAID EpiCov** database on March 19, 2020, 2:53am. Human host data were selected. There were 933 entries (**filename: gisaid_cov2020_sequences_all.fasta**), and if filtered for high coverage only, there were 672 entries (**filename: gisaid_cov2020_sequences_hc.fasta**). For our analysis, we used high coverage data. 

The processed metadata are downloaded from [nextstrain/ncov](https://github.com/nextstrain/ncov/blob/master/data/metadata.tsv). The repo also contains some [quality check notes](https://github.com/nextstrain/ncov/blob/master/config/exclude.txt) which will be confirmed with our own quality check. 

Reference genome is downloaded from [**GenBank**](https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/#nucleotide-sequences):
|GenBank | Gene Region | Colleciton Date | Locality|
|--------|-------------|-----------------|---------| 
|[MN908947](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) | complete | Dec-2019| China |


## Preprocessing
Done by running [quality_check.R](https://github.com/JuliaPalacios/Covid19_Analyses/blob/master/alignment/R/quality_check.R). Eventually, all these preprocessing will be a single bash script but for now, it's a sub-optimal hybrid of bash and R. :pensive: 

1. [preprocessing.sh](https://github.com/JuliaPalacios/Covid19_Analyses/blob/master/alignment/data/preprocessing.sh) script:
	* Format the strain name the same as nextstrain
	* Remove short strain (default cut off = 15000 following nextstrain)
	* Remove obvious duplicates

2. The second part of the code does:
	* Remove strains in nextstrain's [exclude.txt](https://bit.ly/33tuKjc) file
	* Remove strains not in nextstrain's [metadata.tsv](https://bit.ly/2QtA9la) file

Note that the nextstrain's metadata.tsv contains low coverage data as well.

The following four sequences are present in the nextstrain's metadata.tsv but based on the discussions in the virological.org [link 1](https://bit.ly/2Ujod6L) and [link 2](https://bit.ly/2xMDqFC), they may have to be removed too:
|       strain        |  gisaid_epi_isl   |
|----------------------|---------------|
|Shenzhen/SZTH-001/2020 | EPI_ISL_406592 |
|Shenzhen/SZTH-004/2020 | EPI_ISL_406595 |
|Wuhan/IPBCAMS-WH-02/2019 | EPI_ISL_403931 |
|Wuhan/IVDC-HB-04/2020 | EPI_ISL_402120 |

>We considered 51 sequences downloaded from the GISAID database **but we excluded four that may have substantial sequencing errors (EPI_ISL_406592, EPI_ISL_406595, EPI_ISL_403931, and EPI_ISL_402120; see discussion here https://bit.ly/2ROCrfU)**.

>(Andrew Rambaut) There are two genomes with substantial sequencing errors and are omitted in my analysis: EPI_ISL_406592, EPI_ISL_406595. One is your extreme outlier. 

> (Louis du Plessis) In the alignment I used (prepared by Kristian Andersen) there are 6 probable sequencing errors on EPI_ISL_403931 that are masked (marked as ambiguous sites).
> The placement of EPI_ISL_402120 is also different on your tree - I also had it on the big polytomy, you have it on a relatively long branch with several unique mutations. On EPI_ISL_402120 Kristian’s alignment had 96 ambiguous sites (probable sequencing errors, on EPI_ISL_402120 many of these were indels that broke ORFs).

## Alignment
To align it reference genome, all lines of sequence except the last line per named sequence are of the same length. The GenBank reference fasta is formatted for 70 characters per line and the EpiCov data is formatted for 80 characters. (Check with the following command) 
```
sed '2q;d' <file name> | awk "{print length}"
```

Reformat the GenBank reference fasta using **NormalizeFasta** from [picard](https://github.com/broadinstitute/picardi). Here, we will reformate the GenBank reference to match the EpiCov fasta:
```
java -jar picard.jar NormalizeFasta \
      INPUT=input_sequence.fasta \
      OUTPUT=normalized_sequence.fasta \
      LINE_LENGTH=80
``` 

Then concatenate the normalized GenBank reference file and EpiCov file.
```
cat file1 file2 > file3 
``` 

The alignment is performed using **MAFFT** with GenBank's **MN908947** as a reference. A few options tried:
```
mafft --thread -1 mafft_in.fasta > mafft_out.fasta
fftnsi --thread -1 mafft_in.fasta > mafft_out.fasta
mafft --thread -1 --retree 2 --maxiterate 1000 mafft_in.fasta > mafft_out.fasta
```

For details, see the [MAFFT manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html).

The alignment can also be done using **MUSCLE**.

Visualize the alignment using [**MEGA**](https://www.megasoftware.net/), [**SEAVIEW**](http://doua.prabi.fr/software/seaview) or [**Base-By-Base**](https://4virology.net/virology-ca-tools/base-by-base/) (recommended by http://virological.org/).





