##Stuff from Amy, but altered by me
##To find blocks and regions
if (!exists("codedir")) {
  codedir=getwd()
}

plotdir="~/Dropbox/Data/Genetics/Infinium/071312_analysis"
filedir="~/LData/Genetics/Infinium/072512_analysis"

##Gets only plate files and loads in basic stuff
source(file.path(codedir,"450k_cancer_loadin.R"))

library(limma)
library(DNAcopy)

##Rafa functions - God I hate these - prob should move to git and 
##source("/home/bst/faculty/ririzarr/projects/cegs/functions.R")
##source("/home/bst/faculty/ririzarr/projects/cegs/450/functions.R")

##source("/home/bst/faculty/ririzarr/projects/cegs/cis-genome.R")
##source("/home/bst/faculty/ririzarr/projects/cegs/functions.R")

##library(genefilter) #for row Sds()
##library(ggplot2)
##library(RColorBrewer)
##library(matrixStats)


##to get betas                 
##Get just thyroid
data.preload(plates=plates[plates$Tissue %in% "thyroid",],filt.thresh=11)


##save dat file
save(list=c("dat"),file=file.path(filedir, "thy.rda"))

#######################################################
##TAKE OUT SNPS- Chris
x=load(file.path(filedir, "snps_chris.rda"))
snps1=get(x)
snps1=snps1[match(rownames(dat$meth),snps1$IlmnID),]
keepIndex=which(snps1$SBEsnp_RefSeqID=="FALSE"&snps1$CGsnp_RefSeqID=="FALSE")

##Ok - because dat is meth, unmeth, locs, everything (arrays of probes) get rid of annotation at same time
for(i in 1:4) dat[[i]]=dat[[i]][keepIndex,]



#######################################################
# BLOCKS
#######################################################

MG=500 ##max gap for groups to be like charm.. we choose 500 for now.
MNP=3 ##min number of probes per group
BMNP=15


##Concatenate to make single thing
dat$pd$Cat2<-paste(dat$pd$Phenotype,dat$pd$Status,sep="-")
outcome=dat$pd$Cat2

##Timp single descriptor
dat$pd$desc=paste(dat$pd$Tissue, dat$pd$Status, dat$pd$Phenotype, sep="-")



sel=c("normal", "cancer")
block.finding(dat, grps=sel)


pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()
pdf(file.path(plotdir, paste0(sel[1], sel[2],"dmr.pdf")), width=11, height=8.5)
dmr.find(dat, grps=sel)
dev.off()


sel=c("hyperplastic", "cancer")
pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()


sel=c("normal", "hyperplastic")
pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()

sel=c("hyperplastic", "adenoma")
pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()
