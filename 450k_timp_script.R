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

sel=c("normal", "cancer")

block=block.finding(dat, grps=sel, permute.num=0)
dmr=dmr.find(dat, grps=sel)
vmr=vmr.find(dat, grps=sel)

pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()

##Plot dmrs

pdf(file.path(plotdir, paste0(sel[1], sel[2], "dmrggplotd.pdf")), width=11, height=8.5)
range.plot(dat, dmr, grp="pheno")
dev.off()

##Plot blocks now
pdf(file.path(plotdir, paste0(sel[1], sel[2], "blockggplot2.pdf")), width=11, height=8.5)
range.plot(dat, block, grp="pheno")
dev.off()

##Gviz plots!
pdf(file.path(plotdir, paste0(sel[1], sel[2], "dmrgviz.pdf")), width=11, height=8.5)
anno.region.plot(dat, dmr)
dev.off()

##Make plots more flexible(colored lines for different sample types, colored ribbon for iqr range instead of loess 95%?  

#Do a bunch of different tests, ignore permute for now

