##Get path for code on this machine
if (!exists("codedir")) {
  codedir=getwd()
}

source(paste(codedir,"450k_sum_stats.R", sep="/"))
source(paste(codedir,"OU_model.R", sep="/"))

library(minfiLocal)
library(grid)


##load in probes hector sent as csv
setwd("~/Data/Infinium/022412_analysis")

hector.select=read.csv("winstonTab.csv", stringsAsFactors=F)

##select.beta=beta[match(hector.select$X, rownames(beta)),]



##Read in plates my way

##Experimental file location
expdatapath="/thumper2/feinbergLab/core/arrays/illumina/"

plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))

##Read in data
RGset=read.450k.exp(base=expdatapath, targets=plates)

##Get out just cancer samples
tis.samp=(pData(RGset)$Tissue %in% c("colon", "lung", "breast", "thyroid",
  "kidney", "pancreas", "esophagus"))
##Alternatively pData(RGset) gets out the plates variable again

tis.data=preprocessRaw(RGset[,(tis.samp)])

p.list=rownames(getMeth(tis.data))

hector.probes=match(hector.select$X, p.list)

##Just Icr.probes
sel.data=tis.data[hector.probes,]

##Number of samples
samp.tab=tis.pheno(pData(sel.data))

norm.vals=numeric()

##Do per probe tests for each tissue
for (i in 1:dim(samp.tab)[1]) {

  norm.stat=per.stats(sel.data, tissue=rownames(samp.tab)[i], pheno="normal")
  
  norm.vals=cbind(norm.vals, norm.stat$meds)
  
  colnames(norm.vals)[i]=rownames(samp.tab)[i]

}
  
rownames(norm.vals)=p.list[hector.probes]

overall.norm=per.stats(sel.data, tissue=rownames(samp.tab), pheno="normal")$meds

##All tissue seem to be about the same at these probes
tis.var=apply(norm.vals,1, var)



##ok - current idea - run each probe(for whatever reason) lots of times, then generate a hexbin
##of the monte carlo



pdf("Plots/probe_simul.pdf", width=11, height=8.5)
for (i in 1:length(hector.probes)) {
  grid.newpage()
  grid.text(hector.select$X[i], x=0.5, y=.98)
  CpG.plot(sel.data[i,], panel=T, loc=c(0.75, 0.75, 0.4, 0.4))
  ou.hex(loc=c(0.25, 0.75, 0.4, 0.4), mu=overall.norm[i], red=T)
  ou.density(mu=overall.norm[i], panel=T)
  CpG.pheno.density(sel.data[i,], panel=T, loc=c(0.75, 0.25, 0.4, 0.4))
}
dev.off()





