# Covid-19 Genetic Epidemiology Analyses

In this repository we report our code and analyses of publicly available genomic sequences of Human h-Cov19. All the data used here is provided by the laboratories that kindly shared information through GISAID. 

This is an ongoing project. Please refer to [todo](https://github.com/JuliaPalacios/Covid19_Analyses/tree/master/todo) if you would like to contribute on a certain topic.

Data Last update: July 18, 2020.

0. *Data update*. We are processing data on Sherlock. Deposit the msa.fasta and metadata.tsv files from Gisaid into a folder named with the current date in /home/groups/juliapr/covid19/alignment/data/

1. *Data pre-processing*. This code is available on Sherlock: /home/groups/juliapr/covid19/alignment/code/. Edit quality_check.R to change the name of the data folder.
ml R
Rscript quality_check.R > outqual20200718
This code will place two new files in your data folder.

2. *Subset_filter.R* creates a fasta file per population (country and some states)
ml R
Rscript subset_filter.R >outfileter20200718

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

Lorenzo Cappello

## Other Data sources:
https://midasnetwork.us/covid-19/

https://covidtracking.com/data/ 

https://ourworldindata.org/

https://github.com/nytimes/covid-19-data


