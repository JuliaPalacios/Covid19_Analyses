## This script is to subset gisaid by field's value.
# Note that the meta and msa fasta files from GASAID are NOT in the same order...
# Each sequence is in a single line in the current GISAID fasta file.
# Also, the bash script + R combo is suboptimal but this was a quick try.
#
# It assumes you Downloaded msa and metadata from Gisaid and placed them in /alignment/data
#
# This script creates pre-processed one master fasta and meta files for
# all sequences passing the quality check.

rm(list=ls())

## ==================================================
## Setting paths etc
## ==================================================
date <- "20200524"

git_root <- rprojroot::has_file(".git/index")
base.dir <- file.path(git_root$find_file("alignment", "data", date))
script.dir <- file.path(git_root$find_file("alignment", "code"))

setwd(base.dir)
fasta.f <- file.path(base.dir, "msa_0522.fasta")
meta.f <- file.path(base.dir, "metadata.tsv")
getmeta.sh <- file.path(script.dir, "get_meta_all.sh") 
getfasta.sh <- file.path(script.dir, "subset_fasta.sh")

my <- function(str){
  if(.Platform$OS.type == "windows") gsub("/", "//", str)
  else str
}

## ===================================================================
## get header info from fasta file and generate meta file based on it
## ===================================================================
getmeta.str <- paste('sh', my(getmeta.sh), my(base.dir), my(fasta.f), sep=' ')
system(getmeta.str)


## =========================================================
## quality check and get subset index for quality filtering
## =========================================================
fasta.meta <- read.delim(file.path(base.dir, 'meta_fasta.tsv'),
                         as.is=TRUE, header=FALSE)
colnames(fasta.meta) <- c('meta.id', 'gisaid_epi_isl', 'date', 'division')

# check all fasta.meta ID's are in
stopifnot(dim(fasta.meta)[1] == length(grep('EPI_ISL_', fasta.meta$gisaid_epi_isl)))

# check for duplicates (though gisaid said they already filtered for it)
stopifnot(length(fasta.meta$gisaid_epi_isl)
          == length(unique(fasta.meta$gisaid_epi_isl)))

## ==== nextstrain's exclude file ======
download.flag <- T
if (download.flag) {
  cat('\n\n downloading exclude.txt file from nextstrain \n')
  system('curl https://raw.githubusercontent.com/nextstrain/ncov/master/config/exclude.txt -o exclude.txt')
}
system(paste("sed '/#/d; /^$/d' exclude.txt > exclude_ns.txt", sep=''))
exclude.seq <- read.table(file='exclude_ns.txt',
                          header=FALSE, as.is=TRUE, col.names='strain')
to.include.1 <-  !(fasta.meta$gisaid_epi_isl %in% exclude.seq$strain)
cat(paste('\nNumber of sequences NOT in the exclude.txt is', sum(to.include.1), '\n'))


## === remove sequences not in meta file ===
gisaid.meta <- read.delim(meta.f, header=TRUE, as.is=TRUE)
to.include.2 <- fasta.meta$gisaid_epi_isl %in% gisaid.meta$gisaid_epi_isl
cat(paste('Number of sequences in the gisaid metadata is', sum(to.include.2), '\n'))


## === Finally, remove sequences with non-precise date. ===
to.include.3 <- !is.na(strptime(fasta.meta$date, format='%Y-%m-%d'))
cat(paste('Number of sequences with precise date is', sum(to.include.3), '\n'))


## === Final strains to be included ===
to.include <- which(to.include.1 & to.include.2 & to.include.3)
n.tot <- length(to.include)
cat(paste('\nFinal number of sequences after pre-processing is', n.tot, '\n'))
write.table(to.include, paste('tmp_ind.txt', sep=''), quote=FALSE,
            row.names=FALSE, col.names=FALSE)

## ===================================
## Write final processed files
## ===================================
## This part of code is not pretty but does the job and the fastest solution
#  for subsetting a large file for now. Will clean up later.
cat('\nWriting master fasta file and meta file, hang tight (a few mins max)! =)\n')
ind.fasta <- c(rbind(2*to.include-1, 2*to.include))
extract.ind.f <- file.path(base.dir, 'tmp_ind_fasta.txt')
write.table(ind.fasta, file=extract.ind.f,
            col.names=FALSE, row.names=FALSE, quote=FALSE)

getfasta.str <-  paste('sh', my(getfasta.sh), my(base.dir), my(extract.ind.f),
                       'tmp_beast.fasta', 'all_seq.fasta', sep=' ')
system(getfasta.str) # this takes some time (1-2 mins max) to run


## === create corresponding meta file from GISAID meta file ===
match.ind <- match(fasta.meta$gisaid_epi_isl[to.include],
                   gisaid.meta$gisaid_epi_isl)
stopifnot(!any(is.na(match.ind)))
field.include <- c("strain", "gisaid_epi_isl",  "date", "region", 
                   "country", "division", "location",
                   "region_exposure", "country_exposure", "division_exposure",
                   "length", "host")
to.write <- gisaid.meta[match.ind, field.include]
write.table(to.write, file='all_meta.tsv', sep="\t", quote=FALSE, 
            col.names=TRUE, row.names=FALSE)

system('rm tmp* exclude_ns.txt headers_all.txt meta_fasta.tsv')

cat('\n\n DONE!!! \n')






