source("~/Code/timp_illumina/450k_cancer_loadin.R")
setwd("~/Data/Infinium/030812_analysis")

all.tis=c("colon", "lung", "breast", "thyroid", "kidney", "pancreas",
  "esophagus", "liver")


##Get out just cancer samples
tis.samp=(pData(RGset)$Tissue %in% all.tis)

tis.data=preprocessMinfi(RGset[,(tis.samp)])


##Number of samples
samp.tab=tis.pheno(pData(tis.data))
                                   
probe.meds=numeric()              

pdf("norms.pdf")

##Do per probe tests for each tissue
for (i in 1:dim(samp.tab)[1]) {

  norm.stat=per.stats(tis.data, tissue=rownames(samp.tab)[i], pheno="normal")

  probe.meds=cbind(probe.meds, norm.stat$meds)

  colnames(probe.meds)[i]=rownames(samp.tab)[i]

  hist((probe.meds[!is.na(probe.meds)]), main=rownames(samp.tab)[i], plot=T)

}

dev.off()

full.meds=per.stats(tis.data, tissue=all.tis, pheno="normal")$meds

##
middle.probes=which((full.meds<.65)&(full.meds>.35))

all.norm=tis.data[,pData(tis.data)$Phenotype=="normal"]
all.canc=tis.data[,pData(tis.data)$Phenotype=="cancer"]

univ.vinc=incvar.ftest(all.canc, all.norm)
univ.vinc.middle=univ.vinc[(univ.vinc$idx %in% middle.probes),]

fiddy.sig=univ.vinc.middle$idx[univ.vinc.middle$var.pval<1e-15]

pdf("sigs.pdf")


for (i in fiddy.sig) {

##  CpG.plot(tis.data[i,])
  
}

dev.off()







  
