if (!exists("codedir")) {
  codedir=getwd()
}
source(file.path(codedir,"450k_general_init.R"))

##Read in plates
plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))
plates=rbind(plates, read.450k.sheet(expdatapath, "IL04[56].csv$", recursive=T))

##Read in data
RGset=read.450k.exp(base=expdatapath, targets=plates)

##Add something indexing gprobes vs actual data order
values(gprobes)$minfi.idx=match(values(gprobes)$name,
                 rownames(getM(RGset[,1])))


  
  
  
