source("~/Code/timp_illumina/450k_cancer_loadin.R")
setwd("~/Data/Infinium/030812_analysis")

all.tis=c("colon", "lung", "breast", "thyroid", "kidney", "pancreas",
  "esophagus", "liver")


##Get out just cancer samples
tis.samp=(pData(RGset)$Tissue %in% all.tis)

tis.data=preprocessIllumina(RGset[,(tis.samp)])

tis.data=tis.data[good.probes,]

##Number of samples
samp.tab=tis.pheno(pData(tis.data))
                                   
probe.meds=numeric()              

##Do per probe tests for each tissue
for (i in 1:dim(samp.tab)[1]) {

  norm.stat=per.stats(tis.data, tissue=rownames(samp.tab)[i], pheno="normal")

  probe.meds=cbind(probe.meds, norm.stat$meds)

  colnames(probe.meds)[i]=rownames(samp.tab)[i]
}

##
middle.probes=which(rowSums((probe.meds>.3)&(probe.meds<.7))==8)

mid.data=tis.data[middle.probes,]

all.norm=mid.data[,pData(mid.data)$Phenotype=="normal"]
all.canc=mid.data[,pData(mid.data)$Phenotype=="cancer"]

univ.vinc=incvar.ftest(all.canc, all.norm)

fiddy.sig=univ.vinc[univ.vinc$var.pval<1e-15,]

pdf("sigs.pdf")

for (i in 1:50) {

  CpG.plot(mid.data[fiddy.sig$idx[i],])
  
}

dev.off()

anno=gprobes[match(fiddy.sig$probe.name, values(gprobes)$name)]


load("~/Data/Genetics/072111_blocks/lg_regions2.rda")


save(file="fifty.rda", list=c("anno.top1000", "fiddy.sig"))
