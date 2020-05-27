## This script is to subset sequences for specific countries (or regions)
#  and filter sites for missing data within the subsetted sequeces.

rm(list=ls())

library(Rsamtools)
library(Biostrings)

# If you want to add reference genome to the output file, 
# do the following in the "Subset for country" part
# (you can order of the reference genome the way you want):
# headers <- headers.all[c(country.ind, ref.ind)]

# ================================
data.dir <- '~/Documents/Covid_Analysis/'

setwd(data.dir)

fasta.f <- 'all_seq.fasta'
# fasta.f <- 'all_seq_reflen.fasta'
meta.f <- 'all_meta.tsv'

country <- 'USA'

cat(paste('\n================ Processing', country, '================\n'))

# ================================
meta <- read.delim(meta.f, header=TRUE, as.is=TRUE)

if (!file.exists(sprintf("%s.fai", fasta.f))) {
    cat('\n Indexing all sequences. \n')
    indexFa(fasta.f)
    cat('\n Indexing all sequences completed. \n')
}

fa <- open(FaFile(fasta.f))
idx <- scanFaIndex(fa)
n.seq <- countFa(fa)
headers.all <- as.character(seqnames(idx))

ref.id <- 'hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125|2019-12-31'
ref.ind <- which(headers.all == ref.id)
ref.seq <- scanFa(fa, param=idx[ref.ind])[1][[1]]

tot.seq.len <- length(ref.seq)
stopifnot(all(width(ranges(idx)) == tot.seq.len)) # make sure all sequences are the same length
cat(paste('Total number of sequences in', fasta.f, ':', n.seq, '\n'))
cat(paste('Total length of per sequence in', fasta.f, ':', tot.seq.len, '\n\n'))

ref.seq.len <- 29903 # what it should be based on genbank MN908947
if (tot.seq.len == ref.seq.len) {
    cat('All sequences have the same length as the MN908947 reference sequence.\n')
    gap.pos <- unlist(gregexpr(pattern='-', ref.seq))
    stopifnot(gap.pos == -1)
}

## =======================================================
# Subset for country
cat('========== Subsetting for country ===========\n')
country.ind <- which(meta$country == country)
cat(paste('Total number of sequences from', country, 'is',
          length(country.ind), 'of', n.seq, 'sequences.\n\n'))

headers <- headers.all[country.ind]
# headers <- headers.all[c(country.ind, ref.ind)] #uncomment to include ref seq
country.seq <- scanFa(fa, param=GRanges(seqnames=headers,
                                        IRanges(start=1, width=tot.seq.len)))

# country.seq <- scanFa(fa)
# headers <- headers.all

## =======================================================
# Filtering 1
# Remove the sites that have missing or ambiguous values 
# in more than 20% of the samples.

cat('========== Filtering sites ===========\n')
site.th <- 0.2 # missing fraction threshold
sum.mat <- consensusMatrix(country.seq, as.prob=TRUE) 
to.include.site <- which(1 - colSums(sum.mat[1:4, ]) <= site.th)

cat(paste('# of sites with more than', site.th, 
          'of sequences missing or ambiguous are', 
          tot.seq.len - length(to.include.site), 'of',  
          tot.seq.len, 'sites. \n\n'))

# == Additional site-filtering and site subsetting can go here.


# == create site-filtered data
write.range <- reduce(IRanges(start=to.include.site, width=1))

stopifnot(length(country.seq) == length(headers))

filtered <- list()
for (i in 1:length(write.range)) {
    s.name <- headers[i]
    tmp <- subseq(country.seq, write.range[i])
    filtered[[i]] <- tmp
}

cat('Subsetting done. \n')
cat('Gathering sequences. \n\n')
filtered <- DNAStringSet(do.call(paste0, lapply(filtered, as.character)))
names(filtered) <- headers

## =======================================================
# Filtering 2
# Remove the samples that have more than 10% of 
# its sites missing or ambiguous.

cat('========== Filtering sequences ===========\n')
seq.th <- 0.1 # missing fraction threshold
site.stat <- alphabetFrequency(filtered, as.prob=TRUE)
to.include.seq <- which(1 - rowSums(site.stat[ ,1:4]) <= seq.th)

cat(paste('# of sequences from', country, 
          'with more than', seq.th, 'of sites missing or ambiguous are', 
          length(country.seq) - length(to.include.seq), 'of',  
          length(country.seq), 'sequences. \n\n'))

filtered <- filtered[to.include.seq]
headers <- headers[to.include.seq]


# =====================
# ====  write file ====
# =====================
cat('Writing file. \n')
save.f <- paste(country, '.fasta', sep='')
writeXStringSet(filtered, file=save.f, append=FALSE)

cat(paste('\nOutput is written to', data.dir, save.f, '\n\n', sep=''))
cat(paste(save.f, 'contains', length(to.include.seq), 'number of sequences,', 
          'each with length', length(to.include.site), '\n'))

close(fa)

