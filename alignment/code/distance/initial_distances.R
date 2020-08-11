##Coronavirus

##Preliminary molecular clock calculation and distance matrices. It assumes that 
##files per population have been created. We still subset when the file is too big


rm(list=ls())
library(phylodyn)
library(ape)
library(phangorn)
library(lubridate)

#git.dir <- '~/Documents/Covid_Analysis/'
#data.dir <- paste(git.dir, 'alignment/data/', sep='')

git.dir <- '/home/groups/juliapr/covid19/'
date <- '20200720'
data.dir <- paste(git.dir, 'alignment/data/', date, '/distance/', sep='')
code.dir <- paste(git.dir, 'alignment/code/distance/', sep='')
setwd(data.dir)
source(paste(code.dir, "subset_data.R",sep=""))


##Reads the case count files
system('curl https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv -o owid-covid-data.csv')
system('curl https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv -o confirmed.csv')
system('curl https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv -o cases.csv')


meta.f <- '../all_meta.tsv'
meta <- read.delim(meta.f, header=TRUE, as.is=TRUE)

meta$population<-meta$country
meta$population[meta$country=="USA"]<-meta$division[meta$country=="USA"]
meta$population[meta$country=="United Kingdom"]<-meta$division[meta$country=="United Kingdom"]

population.list<-unique(meta$population)

sizes<-table(meta$population)
dates_samp<-ymd(meta$date)
plot(table(dates_samp),cex.axis = 0.4,las=2,xlab="",ylab="Number of samples",main="Sequences by date")

hamminglist <- 0
datelist<-0
lenfast<-0

n.pop <- length(population.list)
cat(paste('\nTotal number of populations to go through is', n.pop, '\n'))
proc.time <- data.frame(time=rep(-1, n.pop), size=rep(-1, n.pop), pop=rep(NA, n.pop))

for (pop.ind in 1:n.pop){
    name <- population.list[pop.ind]
    sz <- sizes[names(sizes)==name]
    proc.time$size[pop.ind] <- sz
    proc.time$pop[pop.ind] <- name
    cat('\n======================================================= \n')
    cat(paste(name, 'has', sz, 'number of sequences'))
    cat('\n======================================================= \n')
    start_time <- Sys.time()
    
    if (sz<1000){#Do nothing special, read the whole dataset
        ##Sequence data (aligned already)
        gisaid.aligned <- paste(data.dir, gsub(" ", "", name, fixed=TRUE), ".fasta", sep='')
        gisaidall <- read.FASTA(gisaid.aligned)
        fastafile <- as.phyDat(gisaidall)
        ref<-length(gisaidall)
        lenfast<-max(lenfast,length(gisaidall[1]))
        hamming1<-as.matrix(dist.hamming(fastafile,ratio=TRUE))[ref,-ref]
        hamminglist<-c(hamminglist,hamming1)
        gendist<-as.matrix(dist.ml(fastafile,model="JC69"))[1:sz,1:sz]
        listout<-list(n=sz,seq_names=names(fastafile)[-ref],hamming=hamming1,distGen=gendist)
        save(listout,file=paste(gsub(" ", "", name, fixed=TRUE),"dist.RData",sep=""))
        datelist<-c(datelist,ymd(meta$date[meta$population==name]))
    }else{ #read by 500s
        gisaid.aligned <- paste(data.dir, gsub(" ", "", name, fixed=TRUE), ".fasta", sep='')
        subsets<-c(seq(0,sz,by=500),sz-1)
        gendist<-matrix(0,sz-1,sz-1)
        hamming1<-0
        for (j in 1:(length(subsets)-2)){
            for (k in (j+1):(length(subsets)-1)){
                if (k==j+1){
                    gisaid.aligned <- paste(data.dir, gsub(" ", "", name, fixed=TRUE), ".fasta", sep='')
                    subset.fasta(data.dir,c(seq(subsets[j]+1,subsets[j+1]),seq(subsets[k]+1,subsets[k+1])),fasta.out = "temp.fasta",fasta.in=gisaid.aligned,include.ref=TRUE)
                    gisaid.aligned <- paste(data.dir,"temp.fasta", sep='')
                    gisaidall <- read.FASTA(gisaid.aligned)
                    fastafile <- as.phyDat(gisaidall)
                    ref<-length(gisaidall)
                    if (j==1){seqnames<-names(fastafile)[-ref]}else{seqnames<-c(seqnames,names(fastafile))}
                    lenfast<-max(lenfast,length(gisaidall[1]))
                    hamming2<-as.matrix(dist.hamming(fastafile,ratio=TRUE))[ref,-ref]
                    hamminglist<-c(hamminglist,hamming2)
                    hamming1<-c(hamming1,hamming2)
                    gendist[c(seq(subsets[j]+1,subsets[j+1]),seq(subsets[k]+1,subsets[k+1])),c(seq(subsets[j]+1,subsets[j+1]),seq(subsets[k]+1,subsets[k+1]))]<-as.matrix(dist.ml(fastafile,model="JC69"))[-length(fastafile),-length(fastafile)]
                    lenfast<-max(lenfast,length(gisaidall[1]))
                }else{
                    gisaid.aligned <- paste(data.dir, gsub(" ", "", name, fixed=TRUE), ".fasta", sep='')
                    subset.fasta(data.dir,c(seq(subsets[j]+1,subsets[j+1]),seq(subsets[k]+1,subsets[k+1])),fasta.out = "temp.fasta",fasta.in=gisaid.aligned,include.ref=FALSE)
                    gisaid.aligned <- paste(data.dir,"temp.fasta", sep='')
                    gisaidall <- read.FASTA(gisaid.aligned)
                    fastafile <- as.phyDat(gisaidall)
                    gendist[c(seq(subsets[j]+1,subsets[j+1]),seq(subsets[k]+1,subsets[k+1])),c(seq(subsets[j]+1,subsets[j+1]),seq(subsets[k]+1,subsets[k+1]))]<-as.matrix(dist.ml(fastafile,model="JC69"))
                    
                }
            }
        }
        datelist<-c(datelist,ymd(meta$date[meta$population==name]))
        listout<-list(n=sz,seq_names=seqnames,hamming=hamming1,distGen=gendist)
        save(listout,file=paste(gsub(" ", "", name, fixed=TRUE),"dist.RData",sep=""))
    }
    end_time <- Sys.time()
    tmp.time <- difftime(end_time, start_time, units='mins')
    print(tmp.time)
    proc.time$time[pop.ind] <- tmp.time
}
list2<-list(hamming=hamminglist[-1],datelist=datelist[-1],lenfast=lenfast)
save(list2, file="hamming.RData")
save(proc.time, file='time.RData') 
# 
# 
# dates_samp<-ymd(meta$date)
# plot(sort(table(dates_samp)),cex.axis = 0.4,las=2,xlab="",ylab="Number of samples",main="Sequences by date")
# 
# fastafile <- as.phyDat(gisaidall)
# #Pairwise number of differences
# hamming<-as.matrix(dist.hamming(fastafile))*as.numeric(summary(gisaidall)[1])
# par(mfrow=c(1,2))
# plot(table(hamming)/2,xlab="Hamming Distance",ylab="Frequency")
# 
# root<-"Wuhan-Hu-1/2019"
# ref<-which(meta_fig$strain==root)
# plot(table(hamming[ref,]),ylab="Frequency",xlab="Hamming distance to reference")
# plot(hamming[ref,],dates_samp-dates_samp[ref],xlab="Distance to reference")
# 
# which.min(dates_samp)
# max.ham<-apply(hamming,1,max)
# tocheck<-which(max.ham>quantile(max.ham,.99))
# tocheck
# i<-9
# par(mfrow=c(2,2))
# plot(table(hamming[tocheck[i],]),ylab="Frequency",xlab="Hamming distance to reference",main=meta_fig[tocheck[i],1])
# plot(table(hamming[tocheck[i+1],]),ylab="Frequency",xlab="Hamming distance to reference",main=meta_fig[tocheck[i+1],1])
# plot(table(hamming[tocheck[i+2],]),ylab="Frequency",xlab="Hamming distance to reference",main=meta_fig[tocheck[i+2],1])
# plot(table(hamming[tocheck[i+3],]),ylab="Frequency",xlab="Hamming distance to reference",main=meta_fig[tocheck[i+3],1])
# 
# ##Hangzhou/ZJU-01/2020|2020-01-25  has the largest divergence to the root
# data.mat<-as.matrix(as.character(as.matrix.DNAbin(gisaidall)))
# differ<-seq(1,ncol(data.mat))[data.mat[tocheck[1],]!=data.mat[tocheck[2],]]
# data.mat[c(ref,tocheck),differ]
# hamming[ref,c(ref,tocheck)]
# #We remove it
# remove<-which(meta_fig$strain=="Hangzhou/ZJU-01/2020")
# hamming<-hamming[-remove,-remove]
# meta_fig<-meta_fig[-remove,]
# fastafile<-fastafile[-remove]
# gisaidall<-gisaidall[-remove]
# plot(table(hamming)/2,xlab="Hamming Distance",ylab="Frequency")
# ref<-which(meta_fig$strain==root)
# 
# 
# plot(table(hamming[ref,]),ylab="Frequency",xlab="Hamming distance to reference")
# max.ham<-apply(hamming,1,max)
# tocheck<-which(max.ham==63)
# plot(table(hamming[tocheck[1],]),ylab="Frequency",xlab="Hamming distance to reference",main=meta_fig[tocheck[1],1])
# plot(table(hamming[tocheck[2],]),ylab="Frequency",xlab="Hamming distance to reference",main=meta_fig[tocheck[2],1])
# data.mat<-as.matrix(as.character(as.matrix.DNAbin(gisaidall)))
# differ<-seq(1,ncol(data.mat))[data.mat[tocheck[1],]!=data.mat[tocheck[2],]]
# differ
# data.mat[c(ref,tocheck),differ]
# #The Denmark/SSI-102/202 seems also too divergent
# 
# 
# 
# 
# tocheck<-which(max.ham==63)
# plot(table(hamming[tocheck[1],]),ylab="Frequency",xlab="Hamming distance to reference",main=meta_fig[tocheck[1],1])
# which.max(hamming[tocheck,])
# 
# 
# data.mat2<-data.mat
# data.mat2<-data.mat[meta_fig$country=="USA",]
# pol<-rep(0,nrow(data.mat2))
# for (j in 1:nrow(data.mat2)){
#   pol[j]<-length(unique(data.mat2[j,]))
# }
# table(pol)
# seq(1,length(pol))[pol==12]
# meta_fig[748,]
# 
# bnp <- BNPR_PS(mcc.tree, prec_alpha=1)
# bnp$grid <- 2020.164 - bnp$grid
# bnp$x <- 2020.164 - bnp$x
# bnp$samp_times <- 2020.164 - bnp$samp_times
# bnp$coal_times <- 2020.164 - bnp$coal_times
# time<-2020.164-bnp$x
# plot_BNPR(bnp, xlim=c(2020,2020.16), heatmap_labels=FALSE)
# plot(time, bnp$effpop, type="l")
# 
# 
# loc <- usa_meta_fig
# plot(table(loc[loc>0]))
# ##distances between samples
# ##This is how I would estimate the rate of mutations per year
# fastafile <- as.phyDat(gisaidus)
# dh<-as.matrix(dist.hamming(fastafile),ratio=FALSE)
# 
# Dmat.l1 <- matrix(NA, nrow=length(dates_samp), ncol=length(dates_samp))
# 
# for (i in 1:length(dates_samp)) {
#   for(j in 1:i){
#     Dmat.l1[i,j] <- Dmat.l1[j,i] <- abs(decimal_date(dates_samp[j])-decimal_date(dates_samp[i]))
#   }
# }
# 
# rat<-dh/Dmat.l1
# mean(rat[Dmat.l1>0])
# # [1] 0.02075419
# ##this is mutation per year
# #Now the number of sites is
# 
# 
# ## ========================
# ## All samples
# ## ========================
# meta_fig <- read.delim(paste(data.dir, 'gisaid_meta_pp.tsv', sep=''), 
#                        header=TRUE, sep='\t', as.is=T)
# 
# ##pariwise hamming distance per time
# 
# 
# ##distances between samples
# ##This is how I would estimate the rate of mutations per year
# fastafile <- as.phyDat(gisaidall)
# dh <- as.matrix(dist.hamming(fastafile))
# 
# Dmat.l1 <- matrix(NA, nrow=length(dates_samp), ncol=length(dates_samp))
# 
# for (i in 1:length(dates_samp)) {
#   for(j in 1:i){
#     Dmat.l1[i,j] <- Dmat.l1[j,i] <- abs(decimal_date(dates_samp[j])-decimal_date(dates_samp[i]))
#   }
# }
# 
# rat<-dh/Dmat.l1
# mean(rat[Dmat.l1>0])
# # [1] 0.01787516  
# ##this is mutation per year
# #Now the number of sites is
# # 5.938986e-07
# # per site
# 
# dm <- dist.ml(fastafile)
# treeUPGMA <- upgma(dm)
# plot(treeUPGMA,show.tip.label = FALSE)
# 
# #Initial Tree
# 
# d.tree(text="((((((((((((((((((((Denmark/SSI-104/2020|2020-03-03,USA/WA-UW52/2020|2020-03-09),(Netherlands/NoordHolland_2/2020|2020-03-03,USA/UPHL-05/2020|2020-03-13)),Finland/FIN03032020C/2020|2020-03-03),Guangdong/20SF174/2020|2020-01-22),Beijing/233/2020|2020-01-28),(((HongKong/case52_VM20002582/2020|2020-02-12,USA/IL2/2020|2020-01-28),Jiangsu/IVDC-JS-001/2020|2020-01-19),USA/WA7-UW4/2020|2020-03-01)),((((((France/IDF0515/2020|2020-01-29,Japan/TY-WK-501/2020|2020-01-31),USA/TX1/2020|2020-02-11),USA/WA-UW60/2020|2020-03-09),((England/20102000506/2020|2020-03-01,Switzerland/1000477757/2020|2020-02-29),CzechRepublic/951/2020|2020-03-01)),((Australia/QLD02/2020|2020-01-30,Netherlands/NA/11/2020|2020-03-10),Switzerland/GR3043/2020|2020-02-27)),France/IDF0372/2020|2020-01-23)),(Wuhan/IVDC-HB-05/2019|2019-12-30,Wuhan/WH03/2020|2020-01-01)),((((((((Guangzhou/GZMU0031/2020|2020-02-25,Netherlands/NA/14/2020|2020-03-10),Shenzhen/HKU-SZ-002/2020|2020-01-10),Taiwan/NTU01/2020|2020-01-31),Guangdong/20SF201/2020|2020-01-23),Nepal/61/2020|2020-01-13),(England/20100022706/2020|2020-02-29,USA/CA4/2020|2020-01-29)),((England/201040081/2020|2020-03-02,USA/WA-UW29/2020|2020-03-08),Guangdong/GDSZ202015-P0019/2020|2020-02-05)),((((Switzerland/1000477796/2020|2020-02-29,Taiwan/4/2020|2020-01-28),Australia/NSW05/2020|2020-02-28),(Beijing/235/2020|2020-01-28,Wuhan/WIV05/2019|2019-12-30)),(Shandong/LY005/2020|2020-01-24,Sweden/01/2020|2020-02-07)))),France/IDF0373/2020|2020-01-23),((((((Brazil/BA-312/2020|2020-03-04,Netherlands/ZuidHolland_22/2020|2020-03-08),USA/WA-UW51/2020|2020-03-08),SouthKorea/KUMC05/2020|2020-02-27),(((England/Sheff01/2020|2020-03-04,Netherlands/NA/35/2020|2020-03-10),Australia/NSW06/2020|2020-02-29),Wuhan/IPBCAMS-WH-01/2019|2019-12-24)),(Italy/SPL1/2020|2020-01-29,USA/WA-UW20/2020|2020-03-05)),(((((Jiangsu/JS02/2020|2020-01-24,Wuhan/WIV07/2019|2019-12-30),France/HF2174/2020|2020-03-09),Switzerland/GE1402/2020|2020-02-28),(((Switzerland/VD0503/2020|2020-02-29,Zhejiang/WZ-02/2020|2020-01-17),Taiwan/2/2020|2020-01-23),Scotland/CVR07/2020|2020-03-09)),((Netherlands/Flevoland/1/2020|2020-03-09,Netherlands/Utrecht_1/2020|2020-03-03),USA/WA-UW32/2020|2020-03-07)))),((((((Netherlands/NoordBrabant_37/2020|2020-03-06,USA/CA9/2020|2020-02-23),Netherlands/NoordBrabant_25/2020|2020-03-06),Shandong/IVDC-SD-001/2020|2020-01-19),(((USA/UPHL-06/2020|2020-03-13,USA/WA-UW44/2020|2020-03-08),Netherlands/NA/2/2020|2020-03-10),(England/200990006/2020|2020-02-26,Netherlands/NA/28/2020|2020-03-12))),USA/CA2/2020|2020-01-22),Taiwan/CGMH-CGU-01/2020|2020-01-25)),(((((((((((Australia/NSW12/2020|2020-03-04,USA/UC-CDPH-UC11/2020|2020-03-05),Canada/BC_41851/2020|2020-03-02),USA/CA-CDPH-UC4/2020|2020-02-27),(((Switzerland/GE0199/2020|2020-02-28,Wuhan/IPBCAMS-WH-04/2019|2019-12-30),France/BFC2147/2020|2020-03-05),(Beijing/231/2020|2020-01-28,Switzerland/1000477377/2020|2020-02-27))),((((England/01/2020|2020-01-29,USA/WA-UW54/2020|2020-03-09),(England/201000003/2020|2020-03-01,Taiwan/CGMH-CGU-03/2020|2020-02-26)),Netherlands/ZuidHolland_20/2020|2020-03-03),Netherlands/Helmond_1363548/2020|2020-02-29)),(((France/PL1643/2020|2020-02-26,Netherlands/NoordHolland/3/2020|2020-03-12),Wuhan/WIV02/2019|2019-12-30),Chongqing/ZX01/2020|2020-01-23)),((China/IQTC02/2020|2020-01-29,USA/WA-UW69/2020|2020-03-10),USA/AZ1/2020|2020-01-22)),(((((Shandong/LY003/2020|2020-01-23,USA/WA-UW25/2020|2020-03-05),Australia/NSW01/2020|2020-01-24),Netherlands/ZuidHolland_13/2020|2020-03-06),(Netherlands/ZuidHolland/28/2020|2020-03-09,Wuhan/WIV06/2019|2019-12-30)),USA/WA-UW66/2020|2020-03-10)),(Netherlands/NoordBrabant_6/2020|2020-03-06,Netherlands/ZuidHolland_19/2020|2020-03-05)),Singapore/2/2020|2020-01-25),((((HongKong/case85_VM20002868/2020|2020-02-24,Netherlands/NoordBrabant/45/2020|2020-03-09),China/WHU02/2020|2020-01-02),(Netherlands/NoordBrabant_39/2020|2020-03-04,Singapore/5/2020|2020-02-06)),(((Singapore/10/2020|2020-02-04,USA/WA-UW49/2020|2020-03-09),USA/UPHL-01/2020|2020-03-13),(China/IQTC01/2020|2020-02-05,Guangzhou/GZMU0048/2020|2020-02-25))))),((((((((((Guangdong/20SF014/2020|2020-01-15,USA/UPHL-04/2020|2020-03-13),Australia/NSW03/2020|2020-01-25),Switzerland/1000477797/2020|2020-02-29),((Hangzhou/HZ-1/2020|2020-01-20,USA/WA-UW72/2020|2020-03-09),USA/CA3/2020|2020-01-29)),Switzerland/GR2988/2020|2020-02-27),((((Canada/BC_40860/2020|2020-03-03,Finland/FIN-455/2020|2020-03-08),England/200990002/2020|2020-02-28),Finland/FIN03032020A/2020|2020-03-03),((Canada/ON-PHL2445/2020|2020-01-25,USA/WA1/2020|2020-01-19),(Foshan/20SF211/2020|2020-01-22,Netherlands/NoordBrabant_20/2020|2020-03-04)))),(((France/HF1870/2020|2020-03-03,Netherlands/ZuidHolland/26/2020|2020-03-09),Switzerland/GE8102/2020|2020-03-01),(Shenzhen/SZTH-002/2020|2020-01-13,USA/CA-CDPH-UC3/2020|2020-02-27))),((((Belgium/VLM-03011/2020|2020-03-03,USA/CA5/2020|2020-01-29),HongKong/VB20024950/2020|2020-01-30),((Netherlands/NA/31/2020|2020-03-13,Netherlands/Utrecht_5/2020|2020-03-02),Switzerland/GE4984/2020|2020-03-07)),USA/WA-UW56/2020|2020-03-09)),(((((((USA/WA9-UW6/2020|2020-03-01,Vietnam/VR03-38142/2020|2020-01-24),France/GE1583/2020|2020-02-25),England/02/2020|2020-01-29),((Netherlands/NA/23/2020|2020-03-09,Switzerland/GE9586/2020|2020-02-27),Hangzhou/ZJU-01/2020|2020-01-25)),((Ireland/Limerick-19935/2020|2020-03-03,Netherlands/NoordBrabant_36/2020|2020-03-03),Australia/NSW02/2020|2020-01-22)),((((Chile/Talca-1/2020|2020-03-02,Wales/PHW37/2020|2020-03-12),Scotland/EDB004/2020|2020-03-04),Yunnan/IVDC-YN-003/2020|2020-01-17),Singapore/9/2020|2020-02-04)),(((Singapore/6/2020|2020-02-09,Wales/PHW1/2020|2020-02-27),USA/WA4-UW2/2020|2020-02-28),Netherlands/NoordBrabant_21/2020|2020-03-04))),((Brazil/SPBR-02/2020|2020-02-28,England/20100121007/2020|2020-02-29),Netherlands/ZuidHolland_24/2020|2020-03-06))),(((((((((Chile/Santiago_op4d1/2020|2020-03-08,Switzerland/GE1422/2020|2020-02-28),France/HF1795/2020|2020-03-02),USA/WA-UW70/2020|2020-03-10),SouthKorea/KCDC06/2020|2020-01-30),((Guangzhou/GZMU0014/2020|2020-02-25,Jiangxi/IVDC-JX-002/2020|2020-01-11),Wuhan/HBCDC-HB-06/2020|2020-02-07)),England/200960515/2020|2020-02-25),Hangzhou/HZCDC0001/2020|2020-01-19),Netherlands/NoordBrabant/58/2020|2020-03-11),((((Netherlands/NA/18/2020|2020-03-09,USA/WA-UW43/2020|2020-03-08),Japan/OS-20-07-1/2020|2020-01-23),(France/BFC2094/2020|2020-03-05,Germany/Baden-Wuerttemberg-1/2020|2020-02-25)),Belgium/VAG-03013/2020|2020-03-01))),(((Netherlands/NA/24/2020|2020-03-08,USA/WA-UW31/2020|2020-03-08),Wuhan/HBCDC-HB-03/2019|2019-12-30),HongKong/case90_VM20002907/2020|2020-02-25)),(((((((Netherlands/Gelderland/3/2020|2020-03-09,Panama/328677/2020|2020-03-06),Netherlands/Utrecht_11/2020|2020-03-03),Foshan/20SF207/2020|2020-01-22),Finland/FIN03032020B/2020|2020-03-03),USA/WA-S3/2020|2020-02-28),((((Australia/NSW08/2020|2020-02-28,Switzerland/AG0361/2020|2020-02-27),France/IDF0386-islP1/2020|2020-01-28),(SouthKorea/KCDC24/2020|2020-02-06,USA/WA13-UW9/2020|2020-03-02)),Shenzhen/HKU-SZ-005/2020|2020-01-11)),(((((HongKong/VM20001061/2020|2020-01-22,USA/WA6-UW3/2020|2020-02-29),England/09c/2020|2020-02-09),(Netherlands/Nootdorp_1364222/2020|2020-03-03,Netherlands/ZuidHolland_16/2020|2020-03-06)),((Netherlands/Tilburg_1364286/2020|2020-03-03,USA/WA-UW26/2020|2020-03-05),Wuhan/WIV04/2019|2019-12-30)),(((France/GE1973/2020|2020-03-04,Netherlands/Gelderland/2/2020|2020-03-09),Guangdong/20SF028/2020|2020-01-17),Taiwan/3/2020|2020-01-24)))),((((((((Belgium/BC-03016/2020|2020-03-01,France/IDF2075/2020|2020-03-02),Guangzhou/GZMU0042/2020|2020-02-25),Netherlands/NA/13/2020|2020-03-10),((Netherlands/ZuidHolland/25/2020|2020-03-09,USA/WA-UW50/2020|2020-03-08),USA/WA-UW76/2020|2020-03-10)),(Netherlands/ZuidHolland_17/2020|2020-03-07,Wuhan/HBCDC-HB-02/2020|2020-01-17)),((Netherlands/NoordBrabant_1/2020|2020-03-02,Netherlands/NoordBrabant_3/2020|2020-03-02),England/200940527/2020|2020-02-25)),Netherlands/ZuidHolland/29/2020|2020-03-09),(((((Netherlands/Haarlem_1363688/2020|2020-03-01,Switzerland/BS0914/2020|2020-03-02),Switzerland/BE2536/2020|2020-03-04),(USA/WA-UW48/2020|2020-03-09,USA/WA-UW68/2020|2020-03-09)),Australia/QLD03/2020|2020-02-05),Netherlands/Utrecht_15/2020|2020-03-08))),(((((((((((((Guangdong/GD2020246-P0028/2020|2020-02-09,Nonthaburi/61/2020|2020-01-08),France/HF1871/2020|2020-03-03),Netherlands/NoordHolland_1/2020|2020-03-03),Netherlands/NoordBrabant_26/2020|2020-03-06),Wuhan/HBCDC-HB-01/2019|2019-12-30),((((Netherlands/NA/16/2020|2020-03-11,Switzerland/GE3895/2020|2020-02-26),Netherlands/NoordBrabant_23/2020|2020-03-05),Switzerland/BL0902/2020|2020-02-27),(Guangdong/20SF012/2020|2020-01-14,Netherlands/NA/1/2020|2020-03-10))),Wales/PHW2/2020|2020-03-04),((England/20099038206/2020|2020-02-29,USA/WA12-UW8/2020|2020-03-03),USA/WA-S2/2020|2020-02-20)),Brazil/SPBR-05/2020|2020-02-29),((((((Belgium/SH-03014/2020|2020-03-01,USA/WA-UW27/2020|2020-03-04),Netherlands/Gelderland_1/2020|2020-03-02),((France/N1620/2020|2020-02-26,USA/WA-UW57/2020|2020-03-09),Chongqing/IVDC-CQ-001/2020|2020-01-18)),((France/IDF1980/2020|2020-03-04,Guangdong/20SF025/2020|2020-01-15),USA/MN2-MDH2/2020|2020-03-07)),((((Switzerland/GE4135/2020|2020-03-06,USA/CA-CDPH-UC9/2020|2020-03-05),USA/CA7/2020|2020-02-06),Wuhan/HBCDC-HB-02/2019|2019-12-30),((Belgium/GHB-03021/2020|2020-02-03,Guangdong/2020XN4459-P0041/2020|2020-01-30),USA/CA-CDPH-UC2/2020|2020-02-27))),(((((Finland/FIN-508/2020|2020-03-07,Netherlands/NA/19/2020|2020-03-12),USA/CA-CDPH-UC6/2020|2020-03-05),Singapore/7/2020|2020-01-27),(China/WHU01/2020|2020-01-02,Jiangsu/JS03/2020|2020-01-24)),(((Nonthaburi/74/2020|2020-01-13,USA/WA-UW18/2020|2020-03-05),Singapore/8/2020|2020-02-03),(Chile/Santiago-1/2020|2020-03-03,HongKong/VB20026565/2020|2020-02-01))))),(((((((Brazil/ES-225/2020|2020-02-29,USA/NY1-PV08001/2020|2020-02-29),Guangzhou/GZMU0044/2020|2020-02-25),(Switzerland/BE6651/2020|2020-02-29,Wuhan/IVDC-HB-04/2020|2020-01-01)),((Denmark/SSI-101/2020|2020-03-03,Switzerland/AG7120/2020|2020-02-29),USA/WA-UW33/2020|2020-03-08)),((((Jiangsu/JS01/2020|2020-01-23,SouthKorea/KCDC12/2020|2020-02-01),England/200990660/2020|2020-02-27),USA/CA1/2020|2020-01-23),USA/WA-UW73/2020|2020-03-10)),((((Australia/NSW11/2020|2020-03-02,Guangdong/GD2020087-P0008/2020|2020-02-01),Netherlands/NA/30/2020|2020-03-13),(Mexico/CDMX-InDRE_01/2020|2020-02-27,SouthKorea/KUMC03/2020|2020-02-27)),(((India/1-31/2020|2020-01-31,USA/WA-UW55/2020|2020-03-09),Wales/PHW28/2020|2020-03-12),(Netherlands/NA/15/2020|2020-03-11,USA/WA-UW22/2020|2020-03-06)))),(((Singapore/4/2020|2020-02-03,Wuhan/WH04/2020|2020-01-05),Wuhan/WH01/2019|2019-12-26),England/200991076/2020|2020-03-01))),(((((((((Netherlands/NA/6/2020|2020-03-10,Switzerland/TI2045/2020|2020-03-01),(Netherlands/NoordBrabant/56/2020|2020-03-09,Netherlands/ZuidHolland_10/2020|2020-03-03)),Singapore/11/2020|2020-02-02),((Finland/1/2020|2020-01-29,USA/WA-UW74/2020|2020-03-10),USA/WA-UW59/2020|2020-03-09)),Switzerland/GE062072020|2020-03-06),Guangzhou/20SF206/2020|2020-01-22),(England/20110003506/2020|2020-03-09,USA/WA-UW41/2020|2020-03-08)),((((Australia/NSW13/2020|2020-03-04,USA/WA2/2020|2020-02-24),France/IDF0626/2020|2020-01-29),Guangdong/2020XN4273-P0036/2020|2020-01-30),Guangdong/20SF040/2020|2020-01-18)),((((((Australia/NSW09/2020|2020-02-28,Georgia/Tb-82/2020|2020-02-28),Netherlands/Utrecht_16/2020|2020-03-08),Japan/TY-WK-521/2020|2020-01-31),((Netherlands/NA/26/2020|2020-03-09,Taiwan/NTU03/2020|2020-03-02),Australia/QLD04/2020|2020-02-05)),(((Italy/UniSR1/2020|2020-03-03,USA/IL1/2020|2020-01-21),(Netherlands/ZuidHolland_15/2020|2020-03-08,Scotland/CVR05/2020|2020-03-04)),Guangdong/2020XN4373-P0039/2020|2020-01-30)),Shandong/LY004/2020|2020-01-26))),(((((((((((Georgia/Tb-468/2020|2020-03-10,Georgia/Tb-477/2020|2020-03-10),Chile/Santiago_op2d1/2020|2020-03-06),(Guangdong/GD2020234-P0023/2020|2020-02-07,SouthKorea/KCDC07/2020|2020-01-31)),Netherlands/NoordBrabant_38/2020|2020-03-06),Wuhan/HBCDC-HB-05/2020|2020-01-18),((((Netherlands/Blaricum_1364780/2020|2020-03-02,Portugal/CV63/2020|2020-03-01),SouthKorea/KCDC05/2020|2020-01-30),Wuhan/IPBCAMS-WH-03/2019|2019-12-30),((Canada/BC_02421/2020|2020-03-01,Netherlands/NoordBrabant_22/2020|2020-03-04),Finland/FIN-266/2020|2020-03-04))),SouthKorea/KCDC03/2020|2020-01-25),(((Belgium/QKJ-03015/2020|2020-03-01,USA/WA-UW63/2020|2020-03-10),Hangzhou/ZJU-05/2020|2020-01-22),(Netherlands/Zeewolde_1365080/2020|2020-03-02,Scotland/CVR10/2020|2020-03-10))),((((((England/20100122106/2020|2020-03-02,Shenzhen/SZTH-003/2020|2020-01-16),Netherlands/Utrecht/17/2020|2020-03-10),(Germany/BavPat3/2020|2020-03-02,USA/MN3-MDH3/2020|2020-03-09)),(((Netherlands/Rotterdam_1364740/2020|2020-03-03,USA/WA-UW24/2020|2020-03-05),USA/WA-UW34/2020|2020-03-08),Foshan/20SF210/2020|2020-01-22)),(Germany/NRW-09/2020|2020-02-28,Netherlands/NA/7/2020|2020-03-09)),Chongqing/YC01/2020|2020-01-21)),((((((Guangdong/GD2020233-P0027/2020|2020-02-07,USA/WA-UW75/2020|2020-03-10),Switzerland/GE3121/2020|2020-02-27),Netherlands/ZuidHolland/27/2020|2020-03-09),Portugal/CV62/2020|2020-03-01),((((Netherlands/Utrecht/18/2020|2020-03-12,USA/WA-UW64/2020|2020-03-09),Belgium/BM-03012/2020|2020-03-01),Netherlands/ZuidHolland_8/2020|2020-03-06),Wuhan-Hu-1/2019|2019-12-26)),((Netherlands/NA/34/2020|2020-03-07,Netherlands/ZuidHolland_23/2020|2020-03-09),Japan/TY-WK-012/2020|2020-01-29))),(((((((((HongKong/case78_VM20002849/2020|2020-02-22,USA/WA-UW23/2020|2020-03-06),Netherlands/NA/21/2020|2020-03-08),(Chile/Santiago-2/2020|2020-03-05,Guangdong/GD2020086-P0021/2020|2020-02-01)),(France/GE1977/2020|2020-03-04,Germany/BavPat2/2020|2020-03-02)),(Germany/NRW-06/2020|2020-02-27,Guangdong/2020XN4448-P0002/2020|2020-01-31)),(((((Scotland/EDB003/2020|2020-03-08,USA/UPHL-03/2020|2020-03-13),Singapore/3/2020|2020-02-01),Shanghai/SH01/2020|2020-02-02),Australia/NSW07/2020|2020-02-29),((Finland/FIN-313/2020|2020-03-05,Georgia/Tb-54/2020|2020-02-27),Shandong/LY007/2020|2020-01-25))),((Fujian/13/2020|2020-01-22,Singapore/1/2020|2020-01-23),USA/WA-UW15/2020|2020-03-04)),Brazil-RJ/314/2020|2020-03-04),((England/200990725/2020|2020-02-28,Netherlands/Utrecht_1363628/2020|2020-03-01),Japan/Hu_DP_Kng_19-020/2020|2020-02-10))))),(((((((((((((((((Belgium/BA-02291/2020|2020-02-29,Chile/Talca-2/2020|2020-03-04),Beijing/105/2020|2020-01-26),Chile/Santiago_op3d1/2020|2020-03-07),Fujian/8/2020|2020-01-21),(Guangdong/2020XN4243-P0035/2020|2020-01-30,Switzerland/GE5373/2020|2020-02-27)),Taiwan/NTU02/2020|2020-02-05),((England/200981386/2020|2020-02-26,Switzerland/SZ1417/2020|2020-03-02),Guangdong/2020XN4475-P0042/2020|2020-01-30)),((Ireland/Limerick-19934/2020|2020-03-03,Netherlands/NA/8/2020|2020-03-09),France/HF1988/2020|2020-03-04)),(((France/HF1805/2020|2020-03-02,Netherlands/NA/12/2020|2020-03-10),USA/WA-UW62/2020|2020-03-09),Zhejiang/WZ-01/2020|2020-01-16)),((Australia/VIC01/2020|2020-01-25,Guangdong/GD2020080-P0010/2020|2020-02-01),SouthKorea/KUMC06/2020|2020-02-27)),((Netherlands/Eindhoven_1363782/2020|2020-03-02,Netherlands/NA/27/2020|2020-03-13),Netherlands/ZuidHolland_18/2020|2020-03-03)),((Netherlands/ZuidHolland/31/2020|2020-03-12,USA/WA3-UW1/2020|2020-02-27),Germany/BavPat1/2020|2020-01-28)),(((((France/HF2039/2020|2020-03-05,Japan/Hu_DP_Kng_19-027/2020|2020-02-10),USA/MA1/2020|2020-01-29),(USA/WA8-UW5/2020|2020-03-01,Wuhan/IPBCAMS-WH-02/2019|2019-12-30)),(Sichuan/IVDC-SC-001/2020|2020-01-15,SouthKorea/KUMC01/2020|2020-02-06)),(((Netherlands/Gelderland/1/2020|2020-03-10,USA/WA-UW28/2020|2020-03-04),Australia/NSW10/2020|2020-02-28),(Cambodia/0012/2020|2020-01-27,HongKong/case42_VM20002493/2020|2020-02-09)))),USA/NY2-PV08100/2020|2020-03-04),(((((((((((((Guangzhou/GZMU0030/2020|2020-02-27,Netherlands/ZuidHolland_2/2020|2020-03-06),USA/WA-UW21/2020|2020-03-05),SouthKorea/KUMC02/2020|2020-02-06),(SouthKorea/KUMC04/2020|2020-02-27,Switzerland/TI9486/2020|2020-02-24)),(((Germany/NRW-10/2020|2020-02-28,Netherlands/NA/9/2020|2020-03-09),(Netherlands/NoordBrabant_17/2020|2020-03-06,USA/WA-UW67/2020|2020-03-09)),France/HF1684/2020|2020-02-29)),(((Netherlands/NA/10/2020|2020-03-09,Netherlands/NA/17/2020|2020-03-09),China/WH-09/2020|2020-01-08),France/HF1993/2020|2020-03-04)),((Australia/QLD09/2020|2020-02-29,England/200990724/2020|2020-02-28),Australia/QLD01/2020|2020-01-28)),(USA/UPHL-02/2020|2020-03-13,Wuhan/HBCDC-HB-03/2020|2020-01-18)),(England/201040141/2020|2020-03-03,USA/WA-UW30/2020|2020-03-08)),((((((Canada/BC_83109/2020|2020-03-05,USA/WA-UW19/2020|2020-03-05),England/200990723/2020|2020-02-27),Tianmen/HBCDC-HB-07/2020|2020-02-08),(Netherlands/Utrecht_1363564/2020|2020-03-01,NewZealand/01/2020|2020-02-27)),((HongKong/case48_VM20002507/2020|2020-02-10,Switzerland/1000477806/2020|2020-02-29),(Luxembourg/Lux1/2020|2020-02-29,Netherlands/NA/25/2020|2020-03-09))),(Netherlands/NoordBrabant/59/2020|2020-03-11,Shandong/LY008/2020|2020-01-30))),((England/200960041/2020|2020-02-27,Japan/NA-20-05-1/2020|2020-01-25),Guangdong/20SF013/2020|2020-01-15)),(((Australia/NSW14/2020|2020-03-03,Japan/TK/20-31-3/2020|2020-02-20),Italy/CDG1/2020|2020-02-20),((Switzerland/GE6679/2020|2020-03-08,Wales/PHW27/2020|2020-03-12),Guangdong/2020XN4433-P0040/2020|2020-01-30))),(Germany/NRW-01/2020|2020-02-28,USA/WA-UW35/2020|2020-03-08))),((((((Netherlands/NA/29/2020|2020-03-13,USA/WA-UW46/2020|2020-03-09),USA/WA-UW71/2020|2020-03-09),France/HF1995/2020|2020-03-04),(((USA/CA-CDPH-UC7/2020|2020-03-05,Wuhan/HBCFailed to load native NucleotideLikelihoodCore in : /Users/juliapr/Documents/Covid19/alignment/BEAST
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
# Using Java nucleotide likelihood core because java.lang.UnsatisfiedLinkError: no NucleotideLikelihoodCore in java.library.path                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             DC-HB-04/2020|2020-01-18),USA/CA8/2020|2020-02-10),Wales/PHW33/2020|2020-03-12)),((Japan/KY-V-029/2020|2020-01-29,USA/CA-PC101P/2020|2020-03-11),Switzerland/VD5615/2020|2020-03-01)),(Netherlands/ZuidHolland_9/2020|2020-03-03,Wuhan/IVDC-HB-01/2019|2019-12-30))),Netherlands/NoordBrabant/55/2020|2020-03-09))
     # tree height = 200.48388627452053
