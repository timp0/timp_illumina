codedir="~/Code/timp_illumina"
plotdir="~/Dropbox/Data/Genetics/Infinium/102212_analysis"
filedir="~/LData/Genetics/Infinium/102212_analysis"

#####NEW RAFA EMAIL INFO

source(file.path(codedir,"450k_general_init.R"))

require(doMC)
registerDoMC()
cores=4
options(cores=cores)

if (file.exists(file.path(filedir, "cancer.rda"))) {
  load(file.path(filedir, "cancer.rda"))
} else {
  
  ##Gets only plate files and loads in basic stuff
  source(file.path(codedir,"450k_cancer_loadin.R"))
  plates=plates[!(plates$Tissue %in% c("cell.line", "mix", "urine", "saliva", "blood", "placenta", "na")),] 
  plates=plates[plates$Phenotype!="bad",]
  
  dat <- read.450k.exp(targets=plates, verbose=TRUE)
  dat=preprocessMinfi(dat)

  ##Add CpG Island anno info
  dat=timp.probeanno(dat)
  
  save(list="dat", file=file.path(filedir, "cancer.rda"), compress="gzip")
}

pd=colData(dat)

Index=which(pd$Tissue=="pancreas" &pd$Phenotype %in% c("normal", "canncer"))

sub=dat[,Index]
subpd=colData(sub)

panc=dat[,which(pd$Tissue=="pancreas" & pd$Phenotype!="bad")]

design=model.matrix(~factor(subpd$Phenotype)+factor(subpd$predictedSex))

res=bumphunter(sub,design,B=100,smooth=FALSE)
sres=bumphunter(sub, design, B=100, smooth=T)

cobj=cpgCollapse(sub)
blocks=blockFinder(cobj$object,design,B=100)


##Plot dmrs
##Change to GRange first
dmr=bump2grange(res$table)

pdf(file.path(plotdir, paste0("dmrggplot.pdf")), width=11, height=8.5)
range.plot(panc, dmr, grp="Phenotype", logit=F)
dev.off()
 
##Gviz plots!
pdf(file.path(plotdir, paste0("dmrgviz.pdf")), width=11, height=8.5)
anno.region.plot(panc, dmr, grp="Phenotype", logit=F)
dev.off()

blocky=blocks$tab

##Plot blocks now
pdf(file.path(plotdir, "blockgg.pdf"), width=11, height=8.5)
range.plot(panc, blocky, grp="Phenotype", logit=F)
dev.off()

pdf(file.path(plotdir, "blockgviz.pdf"), width=11, height=8.5)
anno.region.plot(panc, blocky, grp="Phenotype", logit=F)
dev.off()
