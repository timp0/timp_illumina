##Put us in the right dir for output
setwd("~/Work/Notes/Genomics/Infinium/Instructions/")

##Experimental file location
##expdatapath="/thumper2/feinbergLab/core/iScan/Experiments/"
expdatapath="~/Work/Notes/Genomics/Infinium/Instructions/"

##For Latex Tables
library(xtable)

##Current 450k package name
##require(minfi)
##This will list commands
##ls("package:minfi") 

##Needs this data object apparently?  Let's test!!
##data(IlluminaHumanMethylation450kmanifest)


##Read in plates
plate.id=c("02", "05", "07","08", "09", "10")

targy=NULL


for (i in 1:length(plate.id)) {
  platename=paste("IL0", plate.id[i], sep="")

  ## Read in plate - data is in this long path
  ##basepath=paste(expdatapath, platename, sep="")
  basepath=expdatapath
  targy=rbind(targy,(read.csv(file.path(basepath, paste(platename, "_v2.csv", sep="")), stringsAsFactors=F, skip=5)))
}

##Number of columns
coly=4

#To get vertical lines in 
vlines=rep("r|", coly+1)
vlines[2]="|r|"

##Find different tissue levels
tissue.types=levels(factor(targy$Tissue))
##Fill in multiple of 4 with blanks
tissue.types=c(tissue.types, rep('', coly*ceiling(length(tissue.types)/coly)-length(tissue.types)))
##Make into 4 columns
dim(tissue.types)=c(length(tissue.types)/coly, coly)

##Make table
z=xtable(tissue.types, align=vlines)
print(z, file="tissue1.tex", table.placement="h!", include.rownames=F, include.colnames=F, hline.after=0:nrow(tissue.types), floating=F)

##Find different status levels
status.types=levels(factor(targy$Status))
##Fill in multiple of 4 with blanks
status.types=c(status.types, rep('', coly*ceiling(length(status.types)/coly)-length(status.types)))
##Make into 4 columns
dim(status.types)=c(length(status.types)/coly,coly)

##Make table
z=xtable(status.types, align=vlines)
print(z, file="status1.tex", table.placement="h!", include.rownames=F, include.colnames=F, hline.after=0:nrow(status.types), floating=F)
##Find different Purification levels
Purification.types=levels(factor(targy$Purification))
##Fill in multiple of 4 with blanks
Purification.types=c(Purification.types, rep('', coly*ceiling(length(Purification.types)/coly)-length(Purification.types)))
##Make into 4 columns
dim(Purification.types)=c(length(Purification.types)/coly, coly)

##Make table
z=xtable(Purification.types, align=vlines)
print(z, file="Purification1.tex", table.placement="h!", include.rownames=F,
      include.colnames=F, hline.after=0:nrow(Purification.types), floating=F)






