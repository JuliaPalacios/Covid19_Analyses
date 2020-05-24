###Covid Analysis -May 5


base.dir <- '~/Documents/Covid_Analysis/'
data.dir <- paste(base.dir, 'alignment/data/', sep='')

library(phylodyn)
library(ape)
library(phangorn)
library(lubridate)
source(paste(base.dir,"alignment/code/subset_data.R",sep=""))
source(paste(base.dir,"phylodynamic/function_serial.R",sep=""))

#gisaid.aligned <- paste(data.dir, 'aligned.fasta', sep='')
#gisaidall <- read.FASTA(gisaid.aligned)

meta_fig <- read.delim2(paste(data.dir, '/all_meta.tsv', sep=''), quote="",
                        header=TRUE, sep='\t', as.is=T)

#fastafile <- as.phyDat(gisaidall)

#meta_fig[meta_fig$strain_beast=="India/NCDC-3175/2020|2022-04-05",6]<-"2020-04-05"

countr<-length(unique(meta_fig$country))
total.sam<-nrow(meta_fig)
main1<-paste("A total of ",total.sam," from ",countr, "countries --",date(),collapse="")
plot(sort(table(meta_fig$country)),cex.axis = 0.35,las=2,xlab="",ylab="Number of samples",main="")


#fastafile<-fastafile[meta]
#Read fasta file, choose country. fastafile includes the ancestral reference.
country<-"Algeria"
subs<-c(seq(1,nrow(meta_fig))[meta_fig$country==country],which(meta_fig$strain=="Wuhan-Hu-1/2019"))

fasta.out <- 'test.fasta'
subset.fasta(base.dir,subs, fasta.out)
gisaid.aligned <- paste(data.dir, 'test.fasta', sep='')
gisaidall <- read.FASTA(gisaid.aligned)
fastafile <- as.phyDat(gisaidall)
ref<-which(names(fastafile)=="hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125|2019-12-31")
fastafile<-fastafile[c(seq(1,length(fastafile))[-ref],ref)]
fastafile2<-fastafile[-length(fastafile)]

# If you want to filter some sites for quality issues
# data.matrix<-as.character(as.matrix(fastafile))
# ref<-1
# ref.data<-data.matrix[ref,]
# 
# missing_data<-rep(0,ncol(data.matrix))
# tot.obs<-nrow(data.matrix)
# for (j in 1:ncol(data.matrix)){
#   missing_data[j]<-tot.obs-sum(data.matrix[,j]=="a")-sum(data.matrix[,j]=="c")-sum(data.matrix[,j]=="t")-sum(data.matrix[,j]=="g")
# }
# sum(missing_data>nrow(data.matrix)/5)
# tolerance<-nrow(data.matrix)/5
# ref.data<-ref.data[missing_data<=nrow(data.matrix)/5]
# data.matrix<-data.matrix[,missing_data<=nrow(data.matrix)/5]
# dim(data.matrix)
# 
# fasta1<-as.phyDat(data.matrix)
#Extract sequences, sequence names, and sampling dates from the file
#Save name of the sequence and sampling times from the file

seq_names<-names(fastafile2)
#seq_names[seq_names=="India/NCDC-3175/2020|2022-04-05"]<-"India/NCDC-3175/2020|2020-04-05"

samp_times<-c()
for (r in seq_names){
  samp_times<-c(samp_times,paste(strsplit(r,"|")[[1]][(length(strsplit(r,"|")[[1]])-9):length(strsplit(r,"|")[[1]])],collapse =""))
}
samp_times<-decimal_date(date(samp_times))
lastdate<-max(samp_times)
samp_times<-max(samp_times)-samp_times
name_samp<-cbind(samp_times,seq_names)

# root<-"Wuhan-Hu-1/2019"
# ref2<-which(meta_fig$strain==root)
# 
# meta.to.write <- names(fastafile_1)
# #meta.to.write <- meta.fig$strain_beast[meta.fig$country==chosen]
# #chosen<-"HongKong"
# seqinr::write.fasta(sequences=as.list(data.frame(t(data.matrix))), 
#                     names=meta.to.write,
#                     file.out=paste(base.dir, 'alignment/BEAST/fasta',country,'.fasta', sep=''),
#                     nbchar=60)
# library(seqinr)
# fastafile_1 = read.FASTA(paste(base.dir, 'alignment/BEAST/fasta',country,'.fasta', sep=''))
# fastafile_1<-as.phyDat(fastafile_1)
# 


mu<-mu_linear_reg(fastafile)
mu
mu<-.002
#fastafile_2<-fastafile_2[sample(seq(1,length(fastafile_2)),100)]
#Compute serial UPGMA tree and plot it
#```{r,echo=FALSE}
tree<-serial_upgma(fastafile2,mu, samp_times, name_samp, model="JC69")

par(mfrow=c(2,2))
plot(tree,show.tip.label = FALSE,cex=.3)
axisPhylo()

#Phylodynamic analysis and its plots.
bnp<-BNPR(tree)
bnp_ps<-BNPR_PS(tree)

axlabs<-axis_label(bnp,lastdate,byy=4/365)

country<-state
plot_BNPR2(bnp,axlabs = axlabs,log="",main=country)
plot_BNPR2(bnp_ps,axlabs = axlabs,log="",main=paste(country," PS",sep=""))
#plot_BNPR2(bnp,axlabs = axlabs,log="",main=country,xlim=c(.35,0))
#plot_BNPR2(bnp_ps,axlabs = axlabs,log="",main=paste(country," PS",sep=""),xlim=c(.35,0))

cases<-read.csv(file="~/Documents/Covid_Analysis/phylodynamic/owid-covid-data.csv",header=TRUE)
#cases<-read.csv(file="~/Documents/Covid_Analysis/phylodynamic/Covid19_JohnHopkins.csv",header=TRUE)
#cases_1a<-cases[cases$location=="United States",]
cases_1a<-cases[cases$location==country,]

plot(decimal_date(date(cases_1a$date)),cases_1a$new_cases,type="l",xaxt = "n",ylab="Count Data",xlab="",lwd=2)
mindate<-min(decimal_date(date(cases_1a$date)))
maxdate<-max(decimal_date(date(cases_1a$date)))
axlabs2 = list(x =seq(mindate,maxdate,by=4/365), labs = format(date_decimal(seq(mindate,maxdate,by=4/365)), "%b-%d"),cexlab=.1)
graphics::axis(1, at = axlabs2$x, labels = axlabs2$labs, cex.axis=1,
               las = 1)

abline(v=lastdate,lty=2)
points(decimal_date(date(cases_1a$date)),cases_1a$total_cases,type="l",xaxt = "n",col="blue",lwd=2)

points(decimal_date(date(cases_1a$date)),cases_1a$total_deaths,type="l",xaxt = "n",col="red",lwd=2)

legend(2020.15, y = 15000, bty = "n",c("New cases","Total", "Deaths"), col = c("black","blue","red"), lty=1,lwd=2)


# 
# 
# samp_times<-ymd(meta_fig$date)
# plot(table(decimal_date(date(samp_times))),xaxt = "n",ylab="Sampling Frequency",xlab="",lwd=2)
# mindate<-min(decimal_date(date(samp_times)))
# maxdate<-max(decimal_date(date(samp_times)))
# axlabs2 = list(x =seq(mindate,maxdate,by=4/365), labs = format(date_decimal(seq(mindate,maxdate,by=4/365)), "%b-%d"),cexlab=.1)
# graphics::axis(1, at = axlabs2$x, labels = axlabs2$labs, cex.axis=1,
#                las = 1)
# # 
# date2<-seq
# 
# 
# cases_1a<-cases[cases$province_state==country,]
# 
# date(cases_1a$date)
# cases_1a$date<-as.Date(cases_1a$date,format="%d/%m/%y")
# plot(decimal_date(date(as.Date(cases_1a$date))),cases_1a$new_cases,type="l",xaxt = "n",ylab="Count Data",xlab="",lwd=2)
# mindate<-min(decimal_date(date(cases_1a$date)))
# maxdate<-max(decimal_date(date(cases_1a$date)))
# axlabs2 = list(x =seq(mindate,maxdate,by=4/365), labs = format(date_decimal(seq(mindate,maxdate,by=4/365)), "%b-%d"),cexlab=.1)
# graphics::axis(1, at = axlabs2$x, labels = axlabs2$labs, cex.axis=1,
#                las = 1)
# 
# points(decimal_date(date(cases_1a$date)),cases_1a$total_cases,type="l",xaxt = "n",col="blue",lwd=2)
# 
# points(decimal_date(date(cases_1a$date)),cases_1a$total_deaths,type="l",xaxt = "n",col="red",lwd=2)

