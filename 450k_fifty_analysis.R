source("450k_cancer_loadin.R")

##Get out just cancer samples
tis.samp=(pData(RGset)$Tissue %in% c("colon", "lung", "breast", "thyroid",
  "kidney", "pancreas", "esophagus", "liver"))

tis.data=preprocessIllumina(RGset[,(tis.samp)])


##Number of samples
samp.tab=tis.pheno(pData(tis.data))
                                   
probe.meds=numeric()              
                                   
##Do per probe tests for each tissue
for (i in 1:dim(samp.tab)[1]) {

  norm.stat=per.stats(tis.data, tissue=rownames(samp.tab)[i], pheno="normal")

  probe.meds=cbind(probe.meds, norm.stat$meds)

  colnames(probe.meds)[i]=rownames(samp.tab)[i]

}
  
  
  
