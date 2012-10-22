codedir="~/Code/timp_illumina"
plotdir="~/Dropbox/Data/Genetics/Infinium/071312_analysis"
filedir="~/Dropbox/Data/Genetics/Infinium/071312_analysis"

if (file.exists(file.path(filedir, "thy.rda"))) {
  source(file.path(codedir, "450k_general_init.R"))
  load(file.path(filedir, "thy.rda"))
} else {

  ##Gets only plate files and loads in basic stuff
  source(file.path(codedir,"450k_cancer_loadin.R"))

  ##Get just thy
  dat=dat.preload(plates=plates[plates$Tissue %in% "thyroid",],filt.thresh=11)
  ##TAKE OUT SNPS- Chris
  x=load("/home/bst/faculty/ririzarr/projects/cegs/450/snps_chris.rda")
  ##Winston's way: x=load(file.path("~/Dropbox/Data/Genetics/Infinium/071312_analysis", "snps_chris.rda"))
  snps1=get(x)
  snps1=snps1[match(rownames(dat$meth),snps1$IlmnID),]
  keepIndex=which(snps1$SBEsnp_RefSeqID=="FALSE"&snps1$CGsnp_RefSeqID=="FALSE")
  ##Ok - because dat is meth, unmeth, locs, everything (arrays of probes) get rid of annotation at same time
  for(i in 1:4) dat[[i]]=dat[[i]][keepIndex,]
    
  
  ##Timp single descriptor
  dat$pd$desc=paste(dat$pd$Tissue, dat$pd$Status, dat$pd$Phenotype, sep="-")
  
  ##save dat file
  save(list=c("dat"),file=file.path(filedir, "thy.rda"))
}


dat=dat.init(dat)


if (file.exists(file.path(filedir, "reg.rda"))) {
  load(file.path(filedir, "reg.rda"))
} else {

  block.hyper=block.finding(dat, grps=c("normal", "hyperplastic"), permute.num=100, cores=4)
  block.cancer=block.finding(dat, grps=c("normal", "cancer"), permute.num=100, cores=4)
  

  dmr=dmr.find(dat, grps=c("normal", "cancer"))
  vmr=vmr.find(dat, grps=c("normal", "cancer"))
  
  save(list=c("block.hyper", "block.cancer","dmr", "vmr"), file=file.path(filedir, "reg.rda"), compress="gzip")

}

##Plot blocks now
pdf(file.path(plotdir, "hyperblock.pdf"), width=11, height=8.5)
range.plot(dat, block.hyper, grp="pheno", logit=F)
dev.off()

pdf(file.path(plotdir, "hyperblockgviz.pdf"), width=11, height=8.5)
anno.region.plot(dat, block.hyper, grp="pheno", logit=F)
dev.off()

pdf(file.path(plotdir, "cancerblock.pdf"), width=11, height=8.5)
range.plot(dat, block.cancer, grp="pheno", logit=F)
dev.off()

pdf(file.path(plotdir, "cancerblockgviz.pdf"), width=11, height=8.5)
anno.region.plot(dat, block.cancer, grp="pheno", logit=F)
dev.off()



sel=c("normal", "cancer")

pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()

##Plot dmrs

pdf(file.path(plotdir, paste0(sel[1], sel[2], "dmrggplotd.pdf")), width=11, height=8.5)
range.plot(dat, dmr, grp="pheno")
dev.off()

   

##Gviz plots!
pdf(file.path(plotdir, paste0(sel[1], sel[2], "dmrgviz.pdf")), width=11, height=8.5)
anno.region.plot(dat, dmr)
dev.off()

##Make plots more flexible(colored lines for different sample types, colored ribbon for iqr range instead of loess 95%?  

#Do a bunch of different tests, ignore permute for now

#####NEW RAFA EMAIL INFO

library(minfiLocal)

require(doMC)
registerDoMC()
cores=4
options(cores=cores)


##Gets only plate files and loads in basic stuff
source(file.path(codedir,"450k_cancer_loadin.R"))

##plates=plates[plates$Tissue %in% "thyroid",]
##plates=plates[plates$Tissue %in% "breast",]

plates=plates[!(plates$Tissue %in% c("cell.lines", "mix", "urine", "saliva", "blood", "placenta", "na")),] 

RGset <- read.450k.exp(targets=plates, verbose=TRUE)
object=preprocessMinfi(RGset)
pd=pData(object)

Index=which(pd$Tissue=="thyroid" & pd$Phenotype %in% c("normal", "cancer"))

sub=object[,Index]
subpd=colData(sub)

design=model.matrix(~factor(subpd$Phenotype)+factor(subpd$predictedSex))

res=bumphunter(sub,design,parallel=TRUE,B=25,smooth=FALSE)

cobj=cpgCollapse(sub)

blocks=blockFinder(cobj$object,design,B=100)
