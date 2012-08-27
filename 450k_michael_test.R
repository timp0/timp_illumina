##Stuff from Amy, but altered by me
##To find blocks and regions
if (!exists("codedir")) {
  codedir=getwd()
}

plotdir="~/Dropbox/Data/Genetics/Infinium/082712_analysis"
filedir="~/Dropbox/Data/Genetics/Infinium/082712_analysis"

source(file.path(codedir, "450k_general_init.R"))


if (file.exists(file.path(filedir, "ast.rda"))) {
  load(file.path(filedir, "ast.rda"))
} else {
  target=read.450k.sheet(base=file.path(expdatapath,"IL066"),pattern=".csv$",verbose=TRUE)

  ##Now we subset Micahel's data out of Yun's.
  Index=which(target[,"Experimenter"]=="Michael") 
  pd=target[Index,]

  dat=dat.preload(plates=pd,plotdir=plotdir,expdatapath=expdatapath)

  save(list="dat", file=file.path(filedir, "ast.rda"))
}

