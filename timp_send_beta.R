codedir="~/Code/timp_illumina"
##plotdir="~/Dropbox/Data/Genetics/Infinium/102212_analysis"
filedir="~/LData/Genetics/Infinium/111712_betas"

source(file.path(codedir,"450k_general_init.R"))

require(doMC)
registerDoMC()
cores=4
options(cores=cores)


##Gets only plate files and loads in basic stuff
source(file.path(codedir,"450k_cancer_loadin.R"))

thyroid.plates=plates[plates$Tissue %in% c("thyroid"),]
thyroid.dat <- read.450k.exp(targets=thyroid.plates, verbose=TRUE)
thyroid.dat=preprocessMinfi(thyroid.dat)
thyroid.beta=getBeta(thyroid.dat)
write.csv(thyroid.beta, file=gzfile(file.path(filedir, "thyroid_beta.csv.gz")))

esophagus.plates=plates[plates$Tissue %in% c("esophagus"),]
esophagus.dat <- read.450k.exp(targets=esophagus.plates, verbose=TRUE)
esophagus.dat=preprocessMinfi(esophagus.dat)
esophagus.beta=getBeta(esophagus.dat)
write.csv(esophagus.beta, file=gzfile(file.path(filedir, "esophagus_beta.csv.gz")))

breast.plates=plates[plates$Tissue %in% c("breast"),]
breast.dat <- read.450k.exp(targets=breast.plates, verbose=TRUE)
breast.dat=preprocessMinfi(breast.dat)
breast.beta=getBeta(breast.dat)
write.csv(breast.beta, file=gzfile(file.path(filedir, "breast_beta.csv.gz")))


