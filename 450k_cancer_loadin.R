if (!exists("codedir")) {
  codedir=getwd()
}
source(file.path(codedir,"450k_general_init.R"))

##Read in plate filenames
plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))
plates=rbind(plates, read.450k.sheet(expdatapath, "IL04[56].csv$", recursive=T))



  
  
  
