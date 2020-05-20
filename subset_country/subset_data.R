## This script is to subset gisaid by field's value.
# Note that the meta and msa fasta files from GASAID are NOT in the same order...
# Each sequence is in a single line in the current GISAID fasta file. 
# Also, the bash script + R combo is suboptimal but this was a quick try. 

rm(list=ls())

## ==================================================
## Setting paths etc
## ==================================================
base.dir <- '~/Desktop/Coronavirus/data_GISAID/2020_05_19/msa_0519/'
setwd(base.dir)
fasta.f <- paste(base.dir, 'msa_0519.fasta', sep='')
meta.f <- paste(base.dir, 'metadata_2020-05-19_16-09.tsv', sep='')
country <- 'USA'
getmeta.sh <- paste(base.dir, 'get_meta_all.sh', sep='') # the scripts don't have to be in the base.dir
getfasta.sh <- paste(base.dir, 'get_country_fasta.sh', sep='')

## ===================================================================
## get header info from fasta file and generate meta file based on it
## ===================================================================
getmeta.str <- paste('sh', getmeta.sh, base.dir, fasta.f, country, sep=' ')
system(getmeta.str)

cat(paste('\nCountry chosen:', country, '\n\n'))

## ==================================================
## quality check and get subset index  
## ==================================================
fasta.meta <- read.delim(paste(base.dir, 'meta_fasta.tsv', sep=''), 
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


## === subset country by matching ID in both meta files ===
match.id <- match(fasta.meta$gisaid_epi_isl, gisaid.meta$gisaid_epi_isl)
match.id[which(is.na(match.id))] <- 1 
# it'll be filtered by to.include.2 but just to be sure
stopifnot(!(which(is.na(match.id)) %in% which(to.include.2)))
to.include.4 <- gisaid.meta$country[match.id] == country
cat(paste('Number of sequences from', country, 'is', sum(to.include.4), '\n'))

## Final strains to be included
to.include <- which(to.include.1 & to.include.2 & to.include.3 & to.include.4)
n.tot <- length(to.include)
cat(paste('\nFinal number of sequences after pre-processing and subsetting for',
            country, 'is', n.tot, '\n'))
write.table(to.include, paste('tmp_', country, '_ind.txt', sep=''), quote=FALSE, 
            row.names=FALSE, col.names=FALSE)

## ===================================
## Write final processed files
## ===================================
## This part of code is not pretty but does the job and the fastest solution 
#  for subsetting a large file for now. Will clean up later. 
cat(paste('\nWriting fasta file, meta file, and BEAST file for', country, 
          'hang tight (a few mins max)! =)\n'))
country.ind.fasta <- c(rbind(2*to.include-1, 2*to.include))
extract.ind.f <- paste(base.dir, 'tmp_', country, '_ind_fasta.txt', sep='')
write.table(country.ind.fasta, file=extract.ind.f,
            col.names=FALSE, row.names=FALSE, quote=FALSE)

getfasta.str <-  paste('sh', getfasta.sh, base.dir, extract.ind.f, 
                       fasta.f, country, sep=' ')
system(getfasta.str) # this takes some time (1-2 mins max) to run


## === create corresponding meta file from GISAID meta file ===
match.ind.country <- match(fasta.meta$gisaid_epi_isl[to.include],
                           gisaid.meta$gisaid_epi_isl)
stopifnot(!any(is.na(match.ind.country)))
write.table(gisaid.meta[match.ind.country, ], 
            file=paste(country, '_meta.tsv', sep=''),
            sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

cat('\n\n DONE!!! \n')








