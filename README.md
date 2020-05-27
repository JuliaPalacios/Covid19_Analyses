# Covid-19 Genetic Epidemiology Analyses

In this repository we report our analyses of publicly available genomic sequences of Human h-Cov19. All the data used here is provided by the laboratories that kindly shared information through GISAID. We cannot share the input datafile, however please refer to Data pre-processing to generate the fasta file: aligned.fasta. Pre-processing relies on https://github.com/nextstrain/ncov . We are thankful for their effort and for making their work publicly available.


This is an ongoing project. Please refer to [todo](https://github.com/JuliaPalacios/Covid19_Analyses/tree/master/todo) if you would like to contribute on a certain topic.

1. Data pre-processing. We need to download the msa and metadata files directly from GISAID. Sequences are already aligned. We need to place those files in the data folder. We then need to run *quality_check.R*. This code will place two new files in your data folder.

2. *sites_ref.R* removes sites that are not in the reference sequence and *subset_filter.R* will filter sites that have more than 20% of missing data.

3. Accessing the sequences is slow. You can subset data with *subset.fasta* function in *subset_data.R*.

4. Indexing a large file and accessing it via the indexed file can be faster and more efficient. For an example, see *subset_filter.R*.

5. Quick estimates are generated with *fast_covid.R*.

6. Estimates of Mutation rate


7. Population Structure





8. Diversity (Effective population size) [Preliminary](https://github.com/JuliaPalacios/Covid19_Analyses/blob/master/phylodynamic/Phylodynamic_Analyses1.pdf)




9. Comparative Analyses



10. Comparison to Surveillance Data


## Contributors:

Julia Palacios (juliapr@stanford.edu)

Jaehee Kim 

James Johndrow

Mackenzie Simper

Vladimir Minin

Leonardo Bonanno

Aaron Behr

Samyak Rajanala

## Other Data sources:
https://midasnetwork.us/covid-19/

https://covidtracking.com/data/ 

https://ourworldindata.org/

https://github.com/nytimes/covid-19-data


