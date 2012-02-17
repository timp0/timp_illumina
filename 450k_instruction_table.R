##Put us in the right dir for output
setwd("~/Work/Notes/Genomics/Infinium/Instructions/")

##Experimental file location
expdatapath="/thumper2/feinbergLab/core/iScan/Experiments/"

##For Latex Tables
library(xtable)

##Current 450k package name
require(minfi)
##This will list commands
ls("package:minfi") 

##Needs this data object apparently?  Let's test!!
##data(IlluminaHumanMethylation450kmanifest)


##Read in plates
plate.id=c("02", "05", "07","08", "09", "10")

targy=NULL


for (i in 1:length(plate.id)) {
  platename=paste("IL0", plate.id[i], sep="")

  ## Read in plate - data is in this long path
  basepath=paste(expdatapath, platename, sep="")
  
  targy=rbind(targy,(read.csv(file.path(basepath, paste(platename, "_v2.csv", sep="")), stringsAsFactors=F, skip=5)))
}

##Number of columns
coly=4


tissue.types=levels(factor(targy$Tissue))

tissue.types=c(tissue.types, rep('', coly-(length(tissue.types) %% coly)))

dim(tissue.types)=c(coly, length(tissue.types)/coly)

z=xtable(tissue.types)
print(z, file="try1.tex", table.placement="!htb")

