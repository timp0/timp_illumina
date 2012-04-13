source("~/Code/timp_illumina/450k_cancer_loadin.R")
setwd("~/Data/Infinium/041012_analysis")

all.tis=c("colon", "lung", "breast", "thyroid", "kidney", "pancreas",
  "esophagus", "liver")


##Get out just cancer samples
tis.samp=(pData(RGset)$Tissue %in% all.tis)

tis.data=preprocessIllumina(RGset[,(tis.samp)])

tis.data=tis.data[good.probes,]

norm.data=tis.data[,tis.data$Phenotype=="normal"]

##Number of samples
samp.tab=tis.pheno(pData(tis.data))
                                   
tdmr=list()
tvmr=list()

comp=combn(dim(samp.tab)[1], 2)

for (i in 1:dim(comp)[2]) {
##for (i in 1:1) {
  set1=norm.data[,(norm.data$Tissue==rownames(samp.tab)[comp[1,i]])]

  set2=norm.data[,(norm.data$Tissue==rownames(samp.tab)[comp[2,i]])]
  
  dmr=mean.ttest(set1, set2)
  vmr.u=incvar.ftest(set1, set2)
  vmr.d=incvar.ftest(set2, set1)


  pdf(paste(rownames(samp.tab)[comp[1,i]], rownames(samp.tab)[comp[2,i]],
            "dmr","pdf", sep="."))
  MDS.CpG(norm.data[dmr$idx[1:1e3],], col.p=F)
  for (j in 1:50) {
    CpG.plot(norm.data[dmr$idx[j],], col.p=F)
  }
  dev.off()
  
  pdf(paste(rownames(samp.tab)[comp[1,i]], rownames(samp.tab)[comp[2,i]],
            "vmr.u","pdf", sep="."))
  MDS.CpG(norm.data[vmr.u$idx[1:1e3],], col.p=F)
  for (j in 1:50) {
    CpG.plot(norm.data[vmr.u$idx[j],], col.p=F)
  }
  dev.off()

  pdf(paste(rownames(samp.tab)[comp[1,i]], rownames(samp.tab)[comp[2,i]],
            "vmr.d","pdf", sep="."))
  MDS.CpG(norm.data[vmr.d$idx[1:1e3],], col.p=F)
  for (j in 1:50) {
    CpG.plot(norm.data[vmr.d$idx[j],], col.p=F)
  }
  dev.off()
  
}


