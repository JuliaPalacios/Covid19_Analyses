

subset.fasta <- function(script.dir, data.dir, ind.vec, fasta.in, fasta.out) {
    ## This function subsets all_seq.fasta (or other fasta files) based on
    #  the ind given. 
    ## Input:
    #   script.dir: dir where subset_fasta.sh script resides.
    #   data.dir: dir where fasta.in and fasta.out reside.
    #   ind: A vector of indices of sequences to extract.
    #   fasta.in: name of the input fasta file containing all sequences.
    #   fasta.out: name of output fasta file for subsetted data
    ## Output:
    #   Subsetted fasta in the name of fasta.out in data.dir.
    #   The order of the sequences are sorted(ind) order.
    
    stopifnot(is.vector(ind.vec))
    cat(paste('Subsetting', length(ind.vec), 'number of sequences.\n'))
    cat('\nWriting large file can take a few mins. \n')
    ind.vec <- sort(ind.vec)
    ind.fasta <- c(rbind(2*ind.vec-1, 2*ind.vec))
    extract.ind.f <- paste(data.dir, 'tmp_ind_fasta.txt', sep='')
    write.table(ind.fasta, file=extract.ind.f,
                col.names=FALSE, row.names=FALSE, quote=FALSE)
    
    getfasta.str <-  paste('sh ', script.dir, 'subset_fasta.sh ',
                           data.dir, ' ', extract.ind.f, ' ',
                           fasta.in, ' ', fasta.out, sep='')
    system(getfasta.str) # this takes some time (1-2 mins max) to run
    system(paste('rm ', extract.ind.f))
    cat('\nTada! \n\n')
}

## ==== simple demonstration ====
script.dir <- '~/Documents/Covid_Analysis/alignment/data/'
data.dir <- '~/Documents/Covid_Analysis/alignment/data/'
ind.vec <- c(5,4,2,6)
fasta.in <- 'all_seq.fasta'
fasta.out <- 'test.fasta'

subset.fasta(script.dir, data.dir, ind.vec, fasta.in, fasta.out)

## ==== subsetting for specific country ====
meta.f <- paste(data.dir, 'all_meta.tsv', sep='')
meta <- read.delim(meta.f, as.is=TRUE, sep='\t', header=TRUE)

country <- 'USA'
country.ind <- which(meta$country == country)
fasta.out.2 <- paste(country, '.fasta', sep='')
subset.fasta(script.dir, data.dir, country.ind, fasta.in, fasta.out.2)

## ===== extracting single sequence also works! =====
# Let's extract reference genome
# GenBank: MN908947.3
# GISAID: Wuhan-Hu-1/2019	EPI_ISL_402125
ref.ind <- which(meta$gisaid_epi_isl == 'EPI_ISL_402125')
fasta.out.3 <- 'ref.fasta'
subset.fasta(script.dir, data.dir, ref.ind, fasta.in, fasta.out.3)








