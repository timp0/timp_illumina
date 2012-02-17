##Put us in the right dir
setwd("~/Data/Infinium/080111_analysis/")


##Current 450k package name
require(minfi)



## Read in IL009
basepath="/thumper2/feinbergLab/core/iScan/Experiments/IL009"
nine.targets=read.csv(file.path(basepath, "IL009.csv"), stringsAsFactors=F, skip=7)

##load data
nine.RGset=read.450k.exp(basedir=basepath, nine.targets)


basepath="/thumper2/feinbergLab/core/iScan/Experiments/IL010"
ten.targets=read.csv(file.path(basepath, "IL010.csv"), stringsAsFactors=F, skip=7)
##load data
ten.RGset=read.450k.exp(basedir=basepath, ten.targets)


##Needs this data object apparently?
data(IlluminaHumanMethylation450kmanifest)

##preprocesing
nine.ISet.raw=preprocess.raw(nine.RGset)
ten.ISet.raw=preprocess.raw(ten.RGset)

##nine.Beta(ISet.raw, type="Illumina")[1:4, 1:3]

##normalized(background subtracted)

nine.ISet.norm=preprocess.illumina(nine.RGset, bg.correct=T, normalize="controls", reference=8)
ten.ISet.norm=preprocess.illumina(ten.RGset, bg.correct=T, normalize="controls", reference=8)

##Get out just dilution samples

nin.dil=nine.targets$Sample_Group %in% "Cells_DNA_test"
ten.dil=ten.targets$Sample_Group %in% "Cells_DNA_test"

##Ok - on nine, 36, is 10^5, 75 is 10^4, 37 is 500ng, 76 is 250ng
##on ten, 10-^3 is 16, 10^2 is 53, 10^1 is 87
##on ten 100ng is 17, 50ng is 54, 20ng is 88

##Conc
conc=cbind(Beta(nine.ISet.raw[,c(37, 76)], type="Illumina"), Beta(ten.ISet.raw[,c(17, 54, 88)]))

pdf("concentration.pdf")

smoothScatter(conc[,5], conc[,1])
smoothScatter(conc[,4], conc[,1])
smoothScatter(conc[,3], conc[,1])
smoothScatter(conc[,2], conc[,1])

dev.off()


##Get out betas for conc of DNA
conc=cbind(Beta(nine.ISet.raw[,c(36, 75)], type="Illumina"), Beta(ten.ISet.raw[,c(16, 53, 87)]))

pdf("cells.pdf")

smoothScatter(conc[,5], conc[,1])
smoothScatter(conc[,4], conc[,1])
smoothScatter(conc[,3], conc[,1])
smoothScatter(conc[,2], conc[,1])

dev.off()
