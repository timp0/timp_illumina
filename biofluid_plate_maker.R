##This is simply to make a list of the samples for Biofluids, in a balanced
##Easy way, should be faster than excel

setwd("~/Work/Notes/Genomics/Infinium/Biofluid_Mut")

samp.list=data.frame(tissue=c(rep("oral", 32),rep("blood", 16),
                       rep("urine", 30), rep("saliva", 40),
                       rep("stomach", 30),
                       rep("liver", 11),rep("pancreas", 9),
                       rep("cingulate", 4), rep("frontal.cortex", 4),
                       rep("hippocampus", 4), rep("spleen", 3),
                       rep("liver", 3)),
  state=c(rep("normal", 16), rep("cancer", 32),
    rep("normal", 10), rep("cancer", 20),
    rep("normal", 20), rep("cancer", 20),
    rep("gastritis", 20), rep("cancer", 10),
    rep("cancer", 7), rep("normal", 4),
    rep("cancer", 5), rep("normal", 4),
    rep("normal", 18)),
  phenotype=c(rep("epi.v.gen", 48),
    rep("biofluid", 70), rep("endoscopy.gastritis", 10),
    rep("endoscopy.miscalled", 10), rep("endoscopy.cancer", 10),
    rep("", 11), rep("NET", 5), rep("", 4),
    rep("BSSeq", 18)),
  num=c(1:16, 1:16, 1:16,
    1:10, 1:20,
    1:20, 1:20,
    1:10, 1:10, 1:10,
    1:7, 1:4,
    1:5, 1:4,
    1:4, 1:4, 1:4, 1:3, 1:3),
  un.id=sample(1:1000, 186),
  id=c(rep("UPPP", 16), rep("OralT", 16), rep("O.PBL", 16),
    rep("UrineN", 10), rep("UrineT", 20),
    rep("SalivaN", 20),rep("SalivaT", 20),
    rep("Gastritis", 10), rep("Mis_Gastritis", 10),rep("StomachT", 10),
    rep("LiverT", 7), rep("LiverN", 4),
    rep("PancreasNET", 5), rep("PancreasN", 4),
    rep("CingulateN", 4), rep("FrontalN", 4),
    rep("HippocampusN", 4), rep("SpleenN", 3),
    rep("LiverN", 3)),
  ori=c(rep("RafaGP", 148),
    rep("Timp", 20),
    rep("Sarven", 18))
  )


samp.list$idnum=paste(samp.list$un.id,samp.list$id, samp.list$num, sep="")

samp.list=samp.list[order(samp.list$num),]

write.csv(samp.list, file="listy.csv",row.names=F)
