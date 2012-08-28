if (!exists("codedir")) {
  codedir=getwd()
}
codedir="/thumper2/feinbergLab/personal/aunterma/repos/timp_illumina"
source(file.path("/thumper2/feinbergLab/personal/aunterma/repos/timp_illumina/450k_general_init.R"))

##Read in plates-- Amy's plates
plates=read.450k.sheet("/thumper2/feinbergLab/personal/aunterma/IL058", pattern = "_v2.csv$", recursive = TRUE)
##Read in data
##RGset=read.450k.exp(base=expdatapath, targets=plates)

##Add something indexing gprobes vs actual data order
##values(gprobes)$minfi.idx=match(values(gprobes)$name,
##                 rownames(getM(RGset[,1])))


  
  
  
