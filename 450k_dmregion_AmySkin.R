codedir="~/thumper/repos/timp_illumina"
plotdir="/thumper2/feinbergLab/personal/aunterma/SkinAging/Winston"
filedir="/thumper2/feinbergLab/personal/aunterma/SkinAging/Winston"

if (file.exists(file.path(filedir, "epi.rda"))) {
  source(file.path(codedir, "450k_general_init.R"))
  load(file.path(filedir, "epi.rda"))
} else {

  ##Gets only plate files and loads in basic stuff
  source(file.path(codedir,"450k_cancer_loadin.R"))
  

  ##to get betas                 

  ##Get just epi
  dat=dat.preload(plates=plates[plates$Tissue %in% "epi",],filt.thresh=11)
  ##TAKE OUT SNPS- Chris
  x=load("/home/bst/faculty/ririzarr/projects/cegs/450/snps_chris.rda")
  ##Winston's way: x=load(file.path("~/Dropbox/Data/Genetics/Infinium/071312_analysis", "snps_chris.rda"))
  snps1=get(x)
  snps1=snps1[match(rownames(dat$meth),snps1$IlmnID),]
  keepIndex=which(snps1$SBEsnp_RefSeqID=="FALSE"&snps1$CGsnp_RefSeqID=="FALSE")
  ##Ok - because dat is meth, unmeth, locs, everything (arrays of probes) get rid of annotation at same time
  for(i in 1:4) dat[[i]]=dat[[i]][keepIndex,]
      
  
  ##Concatenate to make single thing
  dat$pd$Cat2<-paste(dat$pd$Category,dat$pd$Status,sep="-")
  outcome=dat$pd$Cat2
  
  ##Timp single descriptor
  ##dat$pd$desc=paste(dat$pd$Tissue, dat$pd$Status, dat$pd$Phenotype, sep="-")
  
  ##save dat file
  save(list=c("dat"),file=file.path(filedir, "epi.rda"))
}
source("~/thumper/repos/timp_illumina/450k_timp_functions.R")
dat=dat.init(dat,refdir="~wtimp/Dropbox/Data/Genetics/Infinium/121311_analysis")

sel=c("Old-exposed", "Young-exposed")
dat$pd$Phenotype=dat$pd$Cat2
block=block.finding(dat, grps=sel, permute.num=0)

dmr1=dmr.find(dat, grps=sel)


pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()

#Plot dmrs

pdf(file.path(plotdir, paste0(sel[1], sel[2], "dmrggplot.pdf")), width=11, height=8.5)
region.plot(dat, dmr, var=note)
dev.off()
write.csv(dmr1,file="/thumper2/feinbergLab/personal/aunterma/SkinAging/Winston/ExposedAgeDMRs.csv")
#Plot blocks now
pdf(file.path(plotdir, paste0(sel[1], sel[2], "blockggplot.pdf")), width=11, height=8.5)
block.plot(dat, block)
dev.off()
#Figure out vmr, vblock
#Do a bunch of different tests, ignore permute for now

