##Put us in the right dir for output
setwd("~/Data/Infinium/080111_analysis/")

##Experimental file location
expdatapath="/thumper2/feinbergLab/core/iScan/Experiments/"

##Current 450k package name
require(minfi)
##This will list commands
ls("package:minfi") 

##Needs this data object apparently?  Let's test!!
##data(IlluminaHumanMethylation450kmanifest)


##Read in plates
plate.id=c("02", "05", "07","08", "09", "10")

targy=NULL

rgdata=list()
targets=list()
preraw=list()
prenorm=list()


for (i in 1:length(plate.id)) {
  platename=paste("IL0", plate.id[i], sep="")

  ## Read in plate - data is in this long path
  basepath=paste(expdatapath, platename, sep="")
  
  ## Read in sample IDs
  ## For this you need rbind
  ##targets[[i]]=read.csv(file.path(basepath, paste(platename, "_v2.csv", sep="")), stringsAsFactors=F, skip=5)

  targy=rbind(targy,(read.csv(file.path(basepath, paste(platename, "_v2.csv", sep="")), stringsAsFactors=F, skip=5)))
  ##targy=rbind(targy, targets[[i]])
  ##load data
  ##rgdata[[i]]=read.450k.exp(basedir=basepath, targets[[i]])

  ##Preprocess Raw
  ##preraw[[i]]=preprocessRaw(rgdata[[i]])

  ##Preprocess Normalized(Background subtracted
  ##prenorm[[i]]=preprocessIllumina(rgdata[[i]], bg.correct=T, normalize="controls", reference=8)
}


##Get out just dilution samples

##nin.dil=nine.targets$Sample_Group %in% "Cells_DNA_test"
##ten.dil=ten.targets$Sample_Group %in% "Cells_DNA_test"

##Ok - on nine, 36, is 10^5, 75 is 10^4, 37 is 500ng, 76 is 250ng
##on ten, 10-^3 is 16, 10^2 is 53, 10^1 is 87
##on ten 100ng is 17, 50ng is 54, 20ng is 88

##Conc
##conc=cbind(Beta(nine.ISet.raw[,c(37, 76)], type="Illumina"), Beta(ten.ISet.raw[,c(17, 54, 88)]))

##pdf("concentration.pdf")

##smoothScatter(conc[,5], conc[,1])
##smoothScatter(conc[,4], conc[,1])
##smoothScatter(conc[,3], conc[,1])
##smoothScatter(conc[,2], conc[,1])

##dev.off()                               


