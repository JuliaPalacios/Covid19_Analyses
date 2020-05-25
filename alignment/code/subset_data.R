

subset.fasta <- function(git.dir, ind.vec, fasta.out) {
    ## This function subsets all_seq.fasta (or other fasta files) 
    #  based on the ind.vec given. For now, ind.vec has to be
    #  sorted and has no duplicate entries.
    ## Input:
    #   git.dir: Home dir for Covid19_Analysis.
    #   ind.vec: A vector of indices of sequences to extract.
    #   fasta.out: Name of output fasta file for subsetted data
    ## Output:
    #   Subsetted fasta in the name of fasta.out in data.dir.
    #   The order of the sequences are sort(unique(ind.vec)) order.
    
    stopifnot(is.vector(ind.vec))
    data.dir <- paste(git.dir, 'alignment/data/', sep='')
    script.dir <- paste(git.dir, 'alignment/code/', sep='')
    fasta.in <- paste(data.dir, 'all_seq.fasta', sep='')
    
    cat(paste('data.dir containing all_seq.fasta is', data.dir, '\n'))
    cat(paste('script.dir containing subset_fasta.sh is', script.dir, '\n\n'))
    
    ind.vec <- sort(unique(ind.vec))
    ind.fasta <- c(rbind(2*ind.vec-1, 2*ind.vec))
    extract.ind.f <- paste(data.dir, 'tmp_ind_fasta.txt', sep='')
    write.table(ind.fasta, file=extract.ind.f,
                col.names=FALSE, row.names=FALSE, quote=FALSE)
    
    cat(paste('Subsetting', length(ind.vec), 'number of sequences from all_seq.fasta \n'))
    cat('\nWriting a large file can take a few mins. \n')
    getfasta.str <-  paste('sh ', script.dir, 'subset_fasta.sh ',
                           data.dir, ' ', extract.ind.f, ' ',
                           fasta.in, ' ', fasta.out, sep='')
    system(getfasta.str) # this takes some time (1-2 mins max) to run
    system(paste('rm ', extract.ind.f))
    cat('\n======== Tada! ========= \n\n')
}
# 
# ## ==== simple demonstration ====
# git.dir <- '~/Documents/Covid_Analysis/'
# 
# ind.vec <- c(5,4,2,6)
# fasta.out <- 'test.fasta'
# 
# subset.fasta(git.dir, ind.vec, fasta.out)
# 
# ## ==== subsetting for specific country ====
# meta.f <- paste(git.dir, 'alignment/data/all_meta.tsv', sep='')
# meta <- read.delim(meta.f, as.is=TRUE, sep='\t', header=TRUE)
# 
# country <- 'USA'
# country.ind <- which(meta$country == country)
# fasta.out.2 <- paste(country, '.fasta', sep='')
# subset.fasta(git.dir, country.ind, fasta.out.2)
# 
# ## ===== extracting single sequence also works! =====
# # Let's extract reference genome
# # GenBank: MN908947.3
# # GISAID: Wuhan-Hu-1/2019	EPI_ISL_402125
# ref.ind <- which(meta$gisaid_epi_isl == 'EPI_ISL_402125')
# fasta.out.3 <- 'ref.fasta'
# subset.fasta(git.dir, ref.ind, fasta.out.3)
# 
# 
# 
# 
# 
# 
# 
# 
