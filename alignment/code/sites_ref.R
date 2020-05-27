## This script removes the inserted gap in the Wuhan reference genome and
#  filter sites of each sequence in all_seq.fasta so that it retains
#  the sites only in the Wuhan reference genome. 
#  Meant to be ran after quality_check.R

rm(list=ls())

library(Rsamtools)
library(Biostrings)

data.dir <- '~/Documents/Covid_Analysis/'
setwd(data.dir)
fasta.f <- 'all_seq.fasta'

if (!file.exists(sprintf("%s.fai", fasta.f))) {
    # need to index the file only once.
    cat('\nIndexing all sequences. \n')
    indexFa(fasta.f)
    cat('Indexing all sequences completed. \n')
}

fa <- open(FaFile(fasta.f))
idx <- scanFaIndex(fa)
n.seq <- countFa(fa)
headers <- as.character(seqnames(idx))

ref.id <- 'hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125|2019-12-31'
ref.ind <- which(headers == ref.id)
ref.seq <- scanFa(fa, param=idx[ref.ind])[1][[1]]

tot.seq.len <- length(ref.seq)
stopifnot(all(width(ranges(idx)) == tot.seq.len)) # make sure all sequences are the same length
cat(paste('\nTotal number of sequences:', n.seq))
cat(paste('\nTotal length of per sequence:', tot.seq.len, '\n\n'))


ref.seq.len <- 29903 # what it should be based on genbank MN908947
gap.pos <- unlist(gregexpr(pattern='-', ref.seq))
to.include <- c(1:tot.seq.len)[-gap.pos]
stopifnot(length(to.include) == ref.seq.len)
write.range <- reduce(IRanges(start=to.include, end=to.include))

# =====================================
# ==== subset sites and write file ====
# =====================================

save.f <- 'all_seq_reflen.fasta'
if (file.exists(save.f)) {
    cat(paste('The file', save.f, 'already exists. Overwriting.\n\n'))
}

## ========================================================================
# Ver.1
# This is faster but more memory intensive. When the dataset gets larger
# or laptop has limited memory, comment this part out and use ver.2 below.

fa.seq <- getSeq(fa)
all.out <- list()
for (i in 1:length(write.range)) {
    s.name <- headers[i]
    tmp <- subseq(fa.seq, write.range[i])
    all.out[[i]] <- tmp
}

cat('Subsetting done. \n')
cat('Gathering sequences. \n')
all.out <- DNAStringSet(do.call(paste0, lapply(all.out, as.character)))
names(all.out) <- headers
cat('Writing file. \n')
writeXStringSet(all.out, file=save.f, append=FALSE)
cat(paste('\nDone! \nOutput is written to ', data.dir, save.f, '\n', sep=''))

## ========================================================================
# Ver.2
# This is a LOT slower than ver 1 (still a few mins), but memory efficient 
# as it doesn't load the whole data into R and just access the file using index.
# There is a more memory efficient version I tried, but it's a lot slower than ver 2.
# I can include the code if necessary.

# all.out <- list()
# for (i in 1:length(write.range)) {
#     tmp.idx <- GRanges(seqnames=headers, ranges=write.range[i])
#     tmp <- scanFa(fa, param=tmp.idx)
#     all.out[[i]] <- tmp
# }
# 
# cat('Subsetting done. \n')
# cat('Gathering sequences. \n')
# all.out <- DNAStringSet(do.call(paste0, lapply(all.out, as.character)))
# names(all.out) <- headers
# cat('Writing file. \n')
# writeXStringSet(all.out, file=save.f, append=FALSE)
# cat(paste('\n Done! \n Output is written to ', data.dir, save.f, '\n', sep=''))

## ========================================================================

close(fa)


