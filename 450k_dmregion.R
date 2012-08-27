##Stuff from Amy, but altered by me
##To find blocks and regions
if (!exists("codedir")) {
  codedir=getwd()
}

plotdir="~/Dropbox/Data/Genetics/Infinium/071312_analysis"
filedir="~/Dropbox/Data/Genetics/Infinium/071312_analysis"

if (file.exists(file.path(filedir, "thy.rda"))) {
  source(file.path(codedir, "450k_general_init.R"))
  load(file.path(filedir, "thy.rda"))
} else {

  ##Gets only plate files and loads in basic stuff
  source(file.path(codedir,"450k_cancer_loadin.R"))
  
  
  ##to get betas                 
  ##Get just thyroid
  dat=dat.preload(plates=plates[plates$Tissue %in% "thyroid",],filt.thresh=11)
  
  ##TAKE OUT SNPS- Chris
  x=load(file.path(filedir, "snps_chris.rda"))
  snps1=get(x)
  snps1=snps1[match(rownames(dat$meth),snps1$IlmnID),]
  keepIndex=which(snps1$SBEsnp_RefSeqID=="FALSE"&snps1$CGsnp_RefSeqID=="FALSE")

  ##Ok - because dat is meth, unmeth, locs, everything (arrays of probes) get rid of annotation at same time
  for(i in 1:4) dat[[i]]=dat[[i]][keepIndex,]
      
  
  ##Concatenate to make single thing
  dat$pd$Cat2<-paste(dat$pd$Phenotype,dat$pd$Status,sep="-")
  outcome=dat$pd$Cat2
  
  ##Timp single descriptor
  dat$pd$desc=paste(dat$pd$Tissue, dat$pd$Status, dat$pd$Phenotype, sep="-")
  
  ##save dat file
  save(list=c("dat"),file=file.path(filedir, "thy.rda"))
}

sel=c("normal", "cancer")
z=block.finding(dat, grps=sel, permute.num=0)

pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()

pdf(file.path(plotdir, paste0(sel[1], sel[2],"dmr.pdf")), width=11, height=8.5)
a=dmr.find(dat, grps=sel)

dev.off()


pdf(file.path(plotdir, paste0(sel[1], sel[2], "dmrggplot.pdf")), width=11, height=8.5)
region.plot(dat, a)
dev.off()

#Plot blocks now
#Figure out vmr, vblock
#Do a bunch of different tests, ignore permute for now
