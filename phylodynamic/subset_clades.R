rm(list=ls())
library(phylodyn)
library(ape)
library(lubridate)
library(phangorn)
library(phytools)
library(ggplot2)
library(ggtree)
library(gridExtra)
setwd("~/Documents/Covid_Analysis/phylodynamic/")
#base.dir <- '~/Desktop/Coronavirus/R/tree_processing/'
#setwd(base.dir)
source("function_serial.R")
#Read fasta file
country<-"Californiadist"
data<-paste("~/Documents/Covid_Analysis/alignment/data/CaliforniaTest/",country,".RData",sep="")
#data<-paste(country,".RData",sep="")
load(data)
distList<-listout

# Extract sequences, sequence names, and sampling dates from the file
#Save name of the sequence and sampling times from the file
n<-distList$n
seq_names<-distList$seq_names
samp_times<-c()
for (r in seq_names){
    samp_times<-c(samp_times,paste(strsplit(r,"|")[[1]][(length(strsplit(r,"|")[[1]])-9):length(strsplit(r,"|")[[1]])],collapse =""))
}
samp_times<-decimal_date(date(samp_times))
lastdate<-max(samp_times)
samp_times2<-max(samp_times)-samp_times
name_samp<-cbind(samp_times2,seq_names)


# Ensure that there are at least two sequences at t=0, otherwise remove the sequence
while(table(samp_times2)[1]==1){
    idx<-which(samp_times2==0)
    #Remove the sequence observed only once
    #fastaformat <- fastaformat[-idx]
    seq_names<-seq_names[-idx]
    samp_times2<-samp_times2[-idx]
    name_samp<-name_samp[-idx,]
    #Make the the last sequence the new 0
    #samp_times<-samp_times-min(samp_times)
    distList$distGen<-distList$distGen[-idx,-idx]
}

# Compute mutation rate. Note: the reference sequence must be in the last spot
mu<-mu_linear_reg_inputDist(distList)

# Compute serial UPGMA tree and plot it
tree<-serial_upgma_inputDist(distList,mu, samp_times2, name_samp)
plot(tree,show.tip.label = FALSE,cex=.3)

# =================================================
# plot and store subtree, interactive
subtr <- subtreeplot(tree, wait=TRUE, show.tip.label=FALSE)

# =================================================
# get number of descendant tips from each internal node
n.tip <- Ntip(tree)
subtr.list <- ape::subtrees(tree, wait=TRUE)
child.ntip <- data.frame(node.id=(1:Nnode(tree))+n.tip,
                         n.tip=unlist(lapply(subtr.list, Ntip)))
n.tip.min <- 500 # the minimum number of tips in a subclade
n.tip.max <- 1000 # the maximum number of tips in a subclade
plt.ind <- which(child.ntip$n.tip > n.tip.min & child.ntip$n.tip < n.tip.max)
plt.node.id <- child.ntip$node.id[plt.ind] #internal node ids
n.plt <- length(plt.node.id)
plt.list <- vector("list", length=n.plt)
for (i in 1:n.plt) {
    print(paste('processing plot', i, 'out of', n.plt))
    plt.list[[i]] <- ggtree(tree) + 
        geom_hilight(node=plt.node.id[i], fill="steelblue", alpha=0.5) +
        theme(axis.title.x=element_text(size=15, face="bold")) +
        xlab(paste('node', plt.node.id[i]))
}

# plot all subtrees meeting the criteria
disp.plt <- marrangeGrob(plt.list, nrow=2, ncol=3)
print(disp.plt)

# Let's say we like internal node 1852.
# Extract tip labels in subclade with the selected node.
root.node <- 1852 # this is internal node index act as a root of a subtree
subtr <- subtr.list[[root.node - n.tip]]
seq.lab <- subtr$tip.label










