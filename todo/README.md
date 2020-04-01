# Covid-19 Genetic Epidemiology Analyses TODO list


- [ ] 1. Access. Everyone should request access to Gisaid.org and have access to viral sequences directly. We will take turns in downloading the data and running [quality_check.R](https://bit.ly/3bIECbE). 

- [ ] 2. Software. Everyone should install 
   1. [BEAST](https://beast.community/) 
   2. [phylodyn](https://github.com/mdkarcher/phylodyn)
   3. [Tracer](https://beast.community/tracer)
   4. TreeAnnotator 
   
- [ ] 3. Fast code for alignment. In quality_check.R, we have used mafft for alignment but as you will see It takes too long. We need to try muscle or clustal to see if it provides a faster alternative.
   - TBD

- [ ] 4. The standard Coalescent model assumes a random sample from the populations. Our samples have population structure (epiclusters). It is better to generate results country-by-country or by epiclusters. One challenge is to find these epiclusters. James started a clustering and there is a movie (with pink background unfortunately). Is there a way to make this "prettier"?
   - TBD

- [ ] 5. Simulation Study. It is not clear how to subsample the whole dataset to have a representative "worldwide" analysis. We want to start with a simulation. We want to simulate two samples of sizes n1 and n2 from a population from a constant trajectory and from a population with exponential growth. What is the optimal assignment of n1 and n2 such that n1+n2=50 (for example) in order to estimate the global population size of the number of infections?
   - Mackenzie

- [ ] 6. Fast Tree. Check out [TreeTime: Maximum-likelihood phylodynamic analysis](https://doi.org/10.1093/ve/vex042), as a first alternative to BEAST runs.
   - TBD

- [ ] 7. Fast initial tree estimation. We want to develop a fast Bayesian EM type of inference. Sequential UPGMA seems to be partially implemented in R. We need someone who can test or implement a version of Sequential UPGMA or Sequential NNI clustering or simply Sequential hierarchical clustering. This will be an initial estimate that can sequentially be optimized according to the posterior distribution. 
   - TBD

- [ ] 8. Estimation of Mutation rate. BEAST *.log files prints joint posterior distribution. We need to have a script that reads all *.log files available in a folder and plots all boxplots of mutation rate per country. 
   - TBD

- [ ] 9. MCCT / Average Phylogenies. We need a script that automatically generates the MCCT with TreeAnnotator from *.trees file output from BEAST.
   - TBD

- [ ] 10. Preferential Sampling. It would be great to accompany all Ne estimations with BNPR_PS.
   - TBD

- [ ] 11. We want to segment genomic regions to see if we can detect some signatures of selection. I would start with USA data only and compare estimations along the genome.
   - Jaehee :nerd_face:
   
- [ ] 12. Comparative Analyses. Can we compare with other Coronavirus data? What is available in Gisaid?
   - Jaehee :nerd_face:
   
- [ ] 13. Surveillance. We want to relate our Ne estimates with actual number of infections. Can someone download a good source of surveillance data? Number of positive and negative tests, number of reported cases, number of hospitalizations? We want a script that compares Ne estimates to surveillance data per country.
   - Jaehee :nerd_face:
   
- [ ] 14. Inference of recombination rate and detection of recombination events. 
   - Jaehee :nerd_face:



