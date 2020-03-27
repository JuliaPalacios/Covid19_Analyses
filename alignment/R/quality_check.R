## March 26, 2020
## This script is to read raw fasta file from GISAID and do quality check,
#  removes any duplicates or entries that don't meet the quality standard. 

## To do
# add in specific dates to final fasta file
# remove entries with ambiguous dates

rm(list=ls())

#base.dir <- '~/Documents/Covid_Analysis/'
base.dir <- '/GitHub/Covid_Analyses/'

data.dir <- paste(base.dir, 'alignment/data/', sep='')

## Automatically downloads necessary files from nextstrain repo
download.flag <- T

ns.dir <- paste(data.dir, 'nextstrain/', sep='') 
if (!dir.exists(ns.dir)) {
    dir.create(ns.dir)
}
if (download.flag) {
    setwd(ns.dir)
    system('curl https://raw.githubusercontent.com/nextstrain/ncov/master/data/metadata.tsv -o metadata.tsv')
    system('curl https://raw.githubusercontent.com/nextstrain/ncov/master/config/exclude.txt -o exclude.txt')
}

##Make sure you have these two files in the folder

gisaid.all.f <- paste(data.dir, 'gisaid_cov2020_sequences.fasta', sep='')
genbank.ref.f <- paste(data.dir, 'genbank_MN908947_ref_norm.fasta', sep='')

###Summary statistics of sequence length
library("ape")
alldata<-read.FASTA(gisaid.all.f)

lengthVec<-rep(0,length(alldata))
for (j in 1:length(alldata)){
  lengthVec[j]<-length(as.matrix(alldata[j]))
}

##Out of the original
length(alldata)
hist(lengthVec)
# This number have length<5000
sum(lengthVec<5000)
hist(lengthVec[lengthVec>5000])
plot(sort(lengthVec))

## Add in date field

## Eventually, all these preprocessing will be a single bash script
#  but for now, it's a sub-optimal hybrid of bash and R... =( 

## ===================================
## pre-processing: first pass
## ===================================
# 1. format the strain name the same as nextstrain
# 2. remove short strain (default cut off = 15000 following nextstrain)
# 3. remove obvious duplicates
setwd(data.dir)
min.len <- 15000
system(paste("bash ./preprocessing.sh ", gisaid.all.f, 
             ' ', min.len, sep=''))

# meta data after first pre-processing
meta.pp1 <- read.delim('meta_gisaid_tmp.txt', header=FALSE, sep='\t', as.is=TRUE,
                       col.names=c('strain', 'gisaid_epi_isl', 'date'))


## ===================================
## pre-processing: second pass
## ===================================
# 1. Remove strains in nextstrain's exclude.txt file
#    This part is done using "strain" field.
# 2. Remove strains not in nextstrain's metadata.tsv file
#    Note that the metadata portion is done with "gisaid_epi_isl" field 
#    instead of "strain" field. The "strain" field in nextstrain format can 
#    have duplicate entries even when "gisaid_epi_isl" is different. 
#    For example:
#    France/N1620/2020     EPI_ISL_414601  2020-02-27
#    France/N1620/2020     EPI_ISL_414624  2020-02-26
#    France/GE1583/2020    EPI_ISL_414600  2020-02-26
#    France/GE1583/2020    EPI_ISL_414623  2020-02-25

## First remove strains in nextstrain's exclude.txt file 
# nextstrain's exclude file, need to remove commented and empty lines
system(paste("sed '/#/d; /^$/d' ", ns.dir, 'exclude.txt', " > ", 
             data.dir, 'exclude_ns.txt', sep=''))
exclude.seq <- read.table(file=paste(data.dir, 'exclude_ns.txt', sep=''),
                          header=FALSE, as.is=TRUE, col.names='strain')
to.include.1 <-  !(meta.pp1$strain %in% exclude.seq$strain)
print(paste('Number of sequences NOT in the exclude.txt is', sum(to.include.1)))

## Next remove strains NOT in nextstrain's metadata.tsv file 
# nextstrain's metadata file
ns.meta.data <- read.delim(paste(ns.dir, 'metadata.tsv', sep=''), 
                           header=TRUE, as.is=TRUE, sep='\t')
to.include.2 <- meta.pp1$gisaid_epi_isl %in% ns.meta.data$gisaid_epi_isl
print(paste('Number of sequences in the metadata.tsv is', sum(to.include.2)))

## Finally, remove sequences with non-precise date.
to.include.3 <- !is.na(strptime(meta.pp1$date, format='%Y-%m-%d'))
print(paste('Number of sequences with precise date is', sum(to.include.3)))

## Final strains to be included
to.include <- to.include.1 & to.include.2 & to.include.3
n.tot <- sum(to.include)
print(paste('Final number of sequences after pre-processing is', n.tot))


## ===================================
## Write final processed files
## ===================================
seq.inc <- meta.pp1$gisaid_epi_isl[to.include]
ns.to.include <- ns.meta.data$gisaid_epi_isl %in% seq.inc
stopifnot(sum(ns.to.include) == n.tot)

## ==== Write metadata file ====
gisaid.meta.out.f <- paste(data.dir, 'gisaid_meta_pp.tsv', sep='')
meta.to.write <- ns.meta.data[ns.to.include, ]

# Insert field for beast strain name
meta.to.write$strain_beast <- paste(meta.to.write$strain, meta.to.write$date, sep='|')
meta.to.write <- meta.to.write[, c(1, 22, 2:21)]
write.table(meta.to.write, file=gisaid.meta.out.f,
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')


## ==== Write fasta file ====
gisaid.fasta.out.f <- paste(data.dir, 'gisaid_seq_pp.fasta', sep='') #regular
gisaid.fasta.out.beast.f <- paste(data.dir, 'gisaid_seq_pp_beast.fasta', sep='') #for beast


# Match index in gisaid fasta file to be in the order in the metadata file
fasta.ind <- rep(NA, n.tot) 
for (i in 1:n.tot) {
    tmp.id <- ns.meta.data$gisaid_epi_isl[which(ns.to.include)[i]]
    fasta.ind[i] <- which(meta.pp1$gisaid_epi_isl == tmp.id)
}
stopifnot(!any(is.na(fasta.ind)))

# Load pre-processed fasta file
pp.fasta.f <- 'seq_gisaid_tmp.fasta'
fasta <- readLines(pp.fasta.f)
header.lines <- which(grepl("^>", fasta, fixed=FALSE))

fasta.list <- list()
for (i in 1:length(header.lines)) {
    if (i == length(header.lines)) {
        tmp.ind <- header.lines[i]:length(fasta)
    } else {
        tmp.ind <- header.lines[i]:(header.lines[i+1] - 1)
    }
    fasta.list[[i]] <- fasta[tmp.ind]
}

fasta.list.order <- fasta.list[fasta.ind]

# check!
for (i in 1:length(fasta.list.order)) {
    tmp <- fasta.list.order[[i]]
    tmp.length <- sum(sapply(2:length(tmp), function(j) nchar(tmp[j])))
    stopifnot(ns.meta.data[ns.to.include, 'length'][i] == tmp.length)
}

fasta.list.order.beast <- list()
for (i in 1:length(fasta.list.order)) {
    tmp <- fasta.list.order[[i]]
    tmp[1] <- paste('>', meta.to.write$strain_beast[i], sep='')
    fasta.list.order.beast[[i]] <- tmp
}

# Write fasta file
# regular
fileConn<-file(gisaid.fasta.out.f)
writeLines(unlist(fasta.list.order), fileConn)
close(fileConn)

# beast
fileConn<-file(gisaid.fasta.out.beast.f)
writeLines(unlist(fasta.list.order.beast), fileConn)
close(fileConn)


## ===================================
## Remove temporary files
## ===================================
system('rm *tmp*')

## ===================================
## Prepare file for MAFFT alignment
## ===================================
# Add normalized GenBank reference sequence to the GISAID fasta file
# For the GenBank sequence normalization, see alignment.md in github.
system('cat genbank_MN908947_ref_norm.fasta gisaid_seq_pp_beast.fasta > mafft_in.fasta')


## ===================================
## MAFFT alignment
## ===================================
# --thread -1 option automatically detects number of cores in computer.
# uncomment below to run mafft

# default method FFT-NS-2 (fast; progressive method):
system('mafft --thread -1 mafft_in.fasta > mafft_out.fasta')

# FFT-NS-i (iterative refinement method; two cycles only):
# system('fftnsi --thread -1 mafft_in.fasta > mafft_out_fftnsi.fasta')

# FFT-NS-i (iterative refinement method; max. 1000 iterations):
# system('mafft --thread -1 --retree 2 --maxiterate 1000 mafft_in.fasta > mafft_out_fftnsi_2.fasta')

## =========================================================
## Remove reference sequence from MAFFT alignment for BEAST
## =========================================================
# The reference sequence occupies first 1,515 lines.
#For the whole alignment I will remove it because it is duplicated
system('sed 1,515d mafft_out.fasta > aligned.fasta')







