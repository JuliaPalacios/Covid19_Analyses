###Covid Analysis -May 5


base.dir <- '~/Documents/Covid_Analysis/'
data.dir <- paste(base.dir, 'alignment/data/', sep='')

library(phylodyn)
library(ape)
library(phangorn)
library(lubridate)
source(paste(base.dir,"phylodynamic/function_serial.R",sep=""))

gisaid.aligned <- paste(data.dir, 'aligned.fasta', sep='')
gisaidall <- read.FASTA(gisaid.aligned)

meta_fig <- read.delim(paste(data.dir, '/gisaid_meta_pp.tsv', sep=''), 
                       header=TRUE, sep='\t', as.is=T)

fastafile <- as.phyDat(gisaidall)


countr<-length(unique(meta_fig$country))
total.sam<-nrow(meta_fig)
main1<-paste("A total of ",total.sam," from ",countr, "countries --",date(),collapse="")
plot(sort(table(meta_fig$country)),cex.axis = 0.35,las=2,xlab="",ylab="Number of samples",main="")


#Read fasta file
country<-"Italy"
subs<-c(seq(1,length(fastafile))[meta_fig$country==country],seq(1,length(fastafile))[meta_fig$strain=="Wuhan-Hu-1/2019"])
fastafile_1<-fastafile[subs]
fastafile_2<-fastafile[meta_fig$country==country]

#Extract sequences, sequence names, and sampling dates from the file
#Save name of the sequence and sampling times from the file
seq_names<-names(fastafile_2)
samp_times<-c()
for (r in seq_names){
  samp_times<-c(samp_times,paste(strsplit(r,"|")[[1]][(length(strsplit(r,"|")[[1]])-9):length(strsplit(r,"|")[[1]])],collapse =""))
}
samp_times<-decimal_date(date(samp_times))
lastdate<-max(samp_times)
samp_times<-max(samp_times)-samp_times
name_samp<-cbind(samp_times,seq_names)
#Compute mutation rate

mu<-mu_linear_reg(fastafile_1)

#Compute serial UPGMA tree and plot it
#```{r,echo=FALSE}
tree<-serial_upgma(fastafile_2,mu, samp_times, name_samp, model="F81")

par(mfrow=c(2,2))
plot(tree,show.tip.label = TRUE,cex=.3)

#Phylodynamic analysis and its plots.
bnp<-BNPR(tree)
bnp_ps<-BNPR_PS(tree)

axlabs<-axis_label(bnp,lastdate,byy=4/365)

plot_BNPR2(bnp,axlabs = axlabs,log="",main=country)
plot_BNPR2(bnp_ps,axlabs = axlabs,log="",main=paste(country," PS",sep=""))

cases<-read.csv(file="~/Documents/Covid_Analysis/phylodynamic/owid-covid-data.csv",header=TRUE)
cases_1a<-cases[cases$location==country,]

plot(decimal_date(date(cases_1a$date)),cases_1a$new_cases,type="l",xaxt = "n",ylab="Count Data",xlab="")
mindate<-min(decimal_date(date(cases_1a$date)))
maxdate<-max(decimal_date(date(cases_1a$date)))
axlabs2 = list(x =seq(mindate,maxdate,by=4/365), labs = format(date_decimal(seq(mindate,maxdate,by=4/365)), "%b-%d"),cexlab=.1)
graphics::axis(1, at = axlabs2$x, labels = axlabs2$labs, cex.axis=1,
               las = 1)

abline(v=lastdate,lty=2)
points(decimal_date(date(cases_1a$date)),cases_1a$total_cases,type="l",xaxt = "n",col="blue")

points(decimal_date(date(cases_1a$date)),cases_1a$total_deaths,type="l",xaxt = "n",col="red")
