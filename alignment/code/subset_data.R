library(Rsamtools)
library(Biostrings)

subset.fasta <- function(git.dir, ind.vec, fasta.out, fasta.in='all_seq.fasta', include.ref=FALSE) {
    ## This function subsets all_seq.fasta (or other fasta files) 
    #  based on the ind.vec given. For now, ind.vec has to be
    #  sorted and has no duplicate entries.
    #  Assumes the original fasta file name is all_seq.fasta
    ## Input:
    #   git.dir: Home dir for Covid19_Analysis.
    #   ind.vec: A vector of indices of sequences to extract.
    #   fasta.out: Name of output fasta file for subsetted data
    #   include.ref: If TRUE, the output file will contain reference sequence
    #                as a first entry in the fasta.out.
    ## Output:
    #   Subsetted fasta in the name of fasta.out in data.dir.
    #   The order of the sequences are sort(unique(ind.vec)) order.
    
    stopifnot(is.vector(ind.vec))
    data.dir <- paste(git.dir, 'alignment/data/', sep='')
    
    cat(paste('data.dir containing all_seq.fasta is', data.dir, '\n'))
    
    setwd(data.dir)
    if (!file.exists(sprintf("%s.fai", fasta.in))) {
        # It needs to be indexed only once.
        cat('\nIndexing all sequences. \n')
        indexFa(fasta.in)
        cat('\nIndexing all sequences completed. \n')
    }
    
    fa <- open(FaFile(fasta.in))
    idx <- scanFaIndex(fa)
    n.seq <- countFa(fa)
    headers.all <- as.character(seqnames(idx))
    
    ref.id <- 'hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125|2019-12-31'
    ref.ind <- which(headers.all == ref.id)
    ref.seq <- scanFa(fa, param=idx[ref.ind])[1][[1]]
    
    tot.seq.len <- length(ref.seq)
    stopifnot(all(width(ranges(idx)) == tot.seq.len)) # make sure all sequences are the same length
    cat(paste('\nTotal number of sequences:', n.seq, '\n'))
    cat(paste('Total length of per sequence:', tot.seq.len, '\n\n'))
    
    ref.seq.len <- 29903 # what it should be based on genbank MN908947
    if (tot.seq.len == ref.seq.len) {
        cat('All sequences have the same length as the MN908947 reference sequence.\n')
        gap.pos <- unlist(gregexpr(pattern='-', ref.seq))
        stopifnot(gap.pos == -1)
    }
    
    cat(paste('Subsetting', length(ind.vec), 'number of sequences from all_seq.fasta \n'))
    if (include.ref) {
        headers <- headers.all[c(ref.ind, ind.vec)] 
    } else {
        headers <- headers.all[ind.vec]
    }
    
    subset.seq <- scanFa(fa, param=GRanges(seqnames=headers, 
                                           IRanges(start=1, width=tot.seq.len)))
    
    cat('Writing file. \n')
    writeXStringSet(subset.seq, file=fasta.out, append=FALSE)
    cat(paste('\nDone! \nOutput is written to', data.dir, fasta.out, '\n', sep=''))
    close(fa)
    cat('\n======== Tada! ========= \n\n')
}
# 
# ## ==== simple demonstration ====
# git.dir <- '~/Documents/Covid_Analysis/'
# 
# ind.vec <- c(5,4,2,6)
# fasta.out <- 'test.fasta'
# fasta.in <- 'all_seq.fasta'
# 
# subset.fasta(git.dir, ind.vec, fasta.out, fasta.in=fasta.in, include.ref=TRUE)
# 
# ## ==== subsetting for specific country ====
# meta.f <- paste(git.dir, 'alignment/data/all_meta.tsv', sep='')
# meta <- read.delim(meta.f, as.is=TRUE, sep='\t', header=TRUE)
# 
# country <- 'USA'
# country.ind <- which(meta$country == country)
# fasta.out.2 <- paste(country, '.fasta', sep='')
# fasta.in.2 <- 'all_seq_reflen.fasta'
# subset.fasta(git.dir, country.ind, fasta.out.2, fasta.in.2)
# 
# ## ===== extracting single sequence also works! =====
# # Let's extract reference genome
# # GenBank: MN908947.3
# # GISAID: Wuhan-Hu-1/2019	EPI_ISL_402125
# ref.ind <- which(meta$gisaid_epi_isl == 'EPI_ISL_402125')
# fasta.out.3 <- 'ref.fasta'
# fasta.in.3 <- 'all_seq_reflen.fasta'
# subset.fasta(git.dir, ref.ind, fasta.out.3, fasta.in.3)
# 
# 
# 
# 
# 
# 
# 
# 
