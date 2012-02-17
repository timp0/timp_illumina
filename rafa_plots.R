source("leekasuccs_ordinal.R")
figdir="figs"
#source("preprocessing.R")

# load dataset
load("~hcorrada/barcode/exps/succs/methyl/data/fulldataset_QN.rda")

source("wrangling.R")

figdir <- "figs1_thyroid"
if (!file.exists(figdir))
  dir.create(figdir)

set.seed(2357)
##ii <- which(dat$ipd2$tissue2=="Colon" & dat$ipd2$subclass!="Adenoma")
pdf("prelim.pdf")
Tissues=unique(dat$ipd2$tissue2)[c(2,1,3,4,5)]
for(h in Tissues){
  mypar(2,2,oma=c(0,0,2,0))
  ii <- which(dat$ipd2$tissue2==h)
  
  y=dat$beta[,ii]
  x=as.numeric(dat$ipd2[ii,"subclass"]!="Normal")
  mds=cmdscale(dist(t(y)))
  plot(mds,col=x+1,main="distance")
  legend("topleft",c("normal","tumor"),pch=16,col=c(1,2))
         
  tt=genefilter::rowttests(y,as.factor(x))
  ms=sapply(split(1:ncol(y),x), function(ind) rowMeans(y[,ind]))
  Ns=sapply(split(1:ncol(y),x),length)
  sds=sapply(split(1:ncol(y),x), function(ind) genefilter::rowSds(y[,ind]))
  f=(sds[,2]/sds[,1])^2
  fp=1-pf(f,Ns[1],Ns[2])
  ps=cbind(tt$p.value,fp)
  colnames(ps)<-c("Mean","Variance")
  barplot(colMeans(ps<0.05)*100,ylim=c(0,100),main="Percent significance")
  plot(ms[,1],ms[,2],col=as.numeric(ps[,1]<0.05)+3,xlab="Normal",ylab="Cancer",main="Mean")
  legend("topleft",c("not signif","signif"),pch=1,col=c(3,4))
  abline(0,1)
  plot(sds[,1],sds[,2],col=as.numeric(ps[,2]<0.05)+3,xlab="Normal",ylab="Cancer",main="Variance",xlim=c(0,0.3),ylim=c(0,0.3))
  abline(0,1)
  mtext(h,side=3,outer=TRUE)
}
dev.off()
