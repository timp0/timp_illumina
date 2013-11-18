source("450k_cancer_loadin.R")

##Put us in the right dir for output
setwd("~/Data/Genetics/Infinium/020312_analysis/")

isl.probes=values(gprobes)$minfi.idx[values(gprobes)$dist.island==0]

tis.samp=(pData(RGset)$Tissue %in% c("colon", "lung", "breast", "thyroid","kidney", "pancreas"))

tis.data=preprocessIllumina(RGset[,tis.samp])

##Find Bad probes



##For loaf

##chosen.probes=good.probes[!(good.probes %in% isl.probes)]
chosen.probes=good.probes

##Get out all loaf
loaf.data=tis.data[chosen.probes,(pData(tis.data)$Tissue %in% c("kidney", "colon"))]
loaf.data$Phenotype[loaf.data$Plate=="IL005"]="loaf"

##Get out all Wilms normal and cancer *that aren't loaf
wilms.data=tis.data[chosen.probes,(pData(tis.data)$Tissue=="kidney")]
wilms.data$Phenotype[wilms.data$Plate=="IL005"]="loaf"

wilms.norm=wilms.data[, pData(wilms.data)$Phenotype=="normal"]
wilms.canc=wilms.data[, pData(wilms.data)$Phenotype=="cancer"]

wilms.m.probes=mean.ttest(wilms.canc, wilms.norm)
wilms.v.probes=incvar.ftest(wilms.canc,wilms.norm)
##wilms.madup=mad.ftest(wilms.canc, wilms.norm)
##wilms.maddown=mad.ftest(wilms.norm, wilms.canc)

##Best difference in variation probes
##Get out all Colon normal and cancer *that aren't loaf

colon.data=tis.data[chosen.probes,pData(tis.data)$Tissue=="colon"]
colon.data$Phenotype[colon.data$Plate=="IL005"]="loaf"

colon.norm=colon.data[,pData(colon.data)$Phenotype=="normal"]
colon.canc=colon.data[,pData(colon.data)$Phenotype=="cancer"]

colon.m.probes=mean.ttest(colon.canc, colon.norm)
colon.v.probes=incvar.ftest(colon.canc,colon.norm)
##colon.v.lower=incvar.ftest(colon.norm, colon.canc)
##colon.madup=mad.ftest(colon.canc, colon.norm)


pdf("colon_more.pdf")
##MDS of colon
MDS.CpG(colon.data[colon.v.probes$idx[1:1e3],])
for (i in 1:50) {
  CpG.plot(colon.data[colon.v.probes$idx[i],])
}
dev.off()
pdf("wilms_more.pdf")
##MDS of wilms
MDS.CpG(wilms.data[wilms.v.probes$idx[1:1e3],])
for (i in 1:50) {
  CpG.plot(wilms.data[wilms.v.probes$idx[i],])
}
dev.off()


loaf.v.probes=cbind(wilms.v.probes,
  colon.v.probes$var.pval[match(colon.v.probes$idx, wilms.v.probes$idx)])

colnames(loaf.v.probes)[3:4]=c("wilms.p", "colon.p")

##Trim off 0 p-values
loaf.v.probes=loaf.v.probes[(loaf.v.probes$wilms.p!=0)&(loaf.v.probes$colon.p!=0),]

##As an estimate let's multiply
loaf.v.probes$both.p=loaf.v.probes$wilms.p*loaf.v.probes$colon.p

loaf.v.probes=loaf.v.probes[order(loaf.v.probes$both.p),]


pdf("whole_loaf.pdf")

##MDS of wilms
MDS.CpG(loaf.data[loaf.v.probes$idx[1:1e3],])
for (i in 1:50) {
  CpG.plot(loaf.data[loaf.v.probes$idx[i],])
}
dev.off()



