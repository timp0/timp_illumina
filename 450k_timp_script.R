codedir="~/Code/timp_illumina"
plotdir="~/Dropbox/Data/Genetics/Infinium/071313_analysis"
filedir="/mithril/homes/timp/LData/Genetics/Infinium/071313_analysis"
expdatapath="/mithril/Data/Infinium"

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
  ##Filter data we don't have
  plates=plates[plates$Basename!="character(0)",]
  plates=plates[plates$Phenotype!="bad",]
  
  dat <- read.450k.exp(targets=plates, verbose=TRUE) 
  dat=preprocessQuantile(dat)

  ##Add CpG Island anno info
  dat=timp.probeanno(dat)
  
  save(list="dat", file=file.path(filedir, "cancer.rda"), compress="gzip")
}

pd=colData(dat)

##Let's get each of the tissue blocks and dmrs

colon.reg=canc.dmrblock(dat, tis="colon")

save(file=file.path(filedir, "colon_cancer.rda"), compress="gzip", list=c("colon.reg"))

breast.reg=canc.dmrblock(dat, tis="breast")

save(file=file.path(filedir, "breast_cancer.rda"), compress="gzip", list=c("breast.reg"))

lung.reg=canc.dmrblock(dat, tis="lung")

save(file=file.path(filedir, "lung_cancer.rda"), compress="gzip", list=c("lung.reg"))

kidney.reg=canc.dmrblock(dat, tis="kidney")

save(file=file.path(filedir, "kidney_cancer.rda"), compress="gzip", list=c("kidney.reg"))


##Make this into a plotting function, so as not to call it constantly

##Colon as a start
coly=dat[,(pd$Tissue=="colon"|((pd$Tissue%in%c("liver","lung"))&(pd$Phenotype!="cancer")))]
colypd=colData(coly)

colypd$anno=NA
colypd$anno[colypd$Tissue=="colon" & colypd$Phenotype=="normal"]="colon.normal"
colypd$anno[colypd$Tissue=="lung" & colypd$Phenotype=="normal"]="lung.normal"
colypd$anno[colypd$Tissue=="liver" & colypd$Phenotype=="normal"]="liver.normal"
colypd$anno[colypd$Tissue=="colon" & colypd$Phenotype=="adenoma"]="colon.adenoma"
colypd$anno[colypd$Tissue=="colon" & colypd$Phenotype=="cancer"]="colon.cancer"
colypd$anno[colypd$Notes=="colon.liver.metastasis"]="colon.liver.metastasis"
colypd$anno[colypd$Notes=="colon.lung.metastasis"]="colon.lung.metastasis"

coly=coly[,!is.na(colypd$anno)]
colData(coly)=colypd[!is.na(colypd$anno),]

compname="colon"

##Plot clusters
pdf(file.path(plotdir, paste0("mds_", compname, ".pdf")), width=11, height=8.5)
cg.cluster(coly, ccomp="anno", grps=c("colon.normal", "colon.cancer"))
dev.off()

pdf(file.path(plotdir, paste0("linkage_", compname, ".pdf")), width=11, height=8.5)
cg.dendro(coly, ccomp="anno", grps=c("colon.normal", "colon.cancer"))
dev.off()

##Plot dmrs

pdf(file.path(plotdir, paste0("dmrggplot_", compname, ".pdf")), width=11, height=8.5)
range.plot(coly, colon.dmr, grp="anno", logit=F)
dev.off()
 
pdf(file.path(plotdir, paste0("dmrgviz_", compname, ".pdf")), width=11, height=8.5)
anno.region.plot(coly, colon.dmr, grp="anno", logit=F)
dev.off()

pdf(file.path(plotdir, paste0("dmrmds_", compname, ".pdf")), width=11, height=8.5)
reg.cluster(coly, colon.dmr, ccomp="anno",namey="phenotype")
dev.off()

##Plot blocks


pdf(file.path(plotdir, paste0("blockmds_", compname, ".pdf")), width=11, height=8.5)
reg.cluster(coly, colon.block, ccomp="anno") 
dev.off()

pdf(file.path(plotdir, paste0("blockgviz_", compname, ".pdf")), width=11, height=8.5)
anno.region.plot(coly, colon.block, grp="anno", logit=F)
dev.off()

pdf(file.path(plotdir, paste0("blockggplot_", compname, ".pdf")), width=11, height=8.5)
range.plot(coly, colon.block, grp="anno", logit=F)
dev.off()

     
Index=which(pd$Tissue=="pancreas" &
      (pd$Phenotype=="normal"|pd$Notes=="adenocarcinoma"))

sub=dat[,Index]
subpd=colData(sub)

panc=dat[,which(pd$Tissue=="pancreas"|
  grepl("pancreas", pd$Notes)|
  (pd$Phenotype=="normal"&pd$Tissue %in% c("lung", "liver")))]

pancpd=colData(panc)

pancpd$anno=NA
##pancpd$anno[pancpd$Notes=="NET"]="NET"
pancpd$anno[pancpd$Notes=="adenocarcinoma"]="adenocarcinoma"
pancpd$anno[pancpd$Phenotype=="metastasis"]=pancpd$Notes[pancpd$Phenotype=="metastasis"]
##pancpd$anno[pancpd$Phenotype=="metastasis"]="metastasis"
pancpd$anno[pancpd$Tissue=="lung"&pancpd$Phenotype=="normal"]="lung.normal"
pancpd$anno[pancpd$Tissue=="liver"&pancpd$Phenotype=="normal"]="liver.normal"
pancpd$anno[pancpd$Phenotype=="normal"&pancpd$Tissue=="pancreas"]="pancreas.normal"
pancpd$anno[pancpd$Phenotype=="adenoma"]="IPMN"

panc=panc[,!is.na(pancpd$anno)]
colData(panc)=pancpd[!is.na(pancpd$anno),]



design=model.matrix(~factor(subpd$Phenotype)+factor(subpd$predictedSex))

res=bumphunter(sub,design,B=100,smooth=FALSE, pickCutoff=T)
sres=bumphunter(sub, design, B=100, pickCutoff=T)

cobj=cpgCollapse(sub)
blocks=blockFinder(cobj$object,design,B=100, pickCutoff=T)


##Plot clusters
pdf(file.path(plotdir, paste0("mds.pdf")), width=11, height=8.5)
cg.cluster(panc, ccomp="anno", grps=c("pancreas.normal", "adenocarcinoma"))
dev.off()

pdf(file.path(plotdir, paste0("linkage.pdf")), width=11, height=8.5)
cg.dendro(panc, ccomp="anno", grps=c("pancreas.normal", "adenocarcinoma"))
dev.off()

##Plot dmrs
##Change to GRange first
dmr=bump2grange(res$table)
sdmr=bump2grange(sres$table)

pdf(file.path(plotdir, paste0("sdmrggplot.pdf")), width=11, height=8.5)
range.plot(panc, sdmr, grp="anno", logit=F)
dev.off()
 
##Gviz plots!
pdf(file.path(plotdir, paste0("sdmrgviz.pdf")), width=11, height=8.5)
anno.region.plot(panc, sdmr, grp="anno", logit=F)
dev.off()

pdf(file.path(plotdir, paste0("sdmrmds.pdf")), width=11, height=8.5)
reg.cluster(panc, sdmr, ccomp="anno",namey="phenotype")
dev.off()


blocky=bump2grange(blocks$tab)

pdf(file.path(plotdir, paste0("blockmds.pdf")), width=11, height=8.5)
reg.cluster(panc, blocky, ccomp="anno") 
dev.off()

pdf(file.path(plotdir, "blockgviz.pdf"), width=11, height=8.5)
anno.region.plot(panc, blocky, grp="anno", logit=F)
dev.off()

mutdir="~/Dropbox/Data/Genetics/Mutations/092512_dl"
#Mutation tables
load(file.path(mutdir, "muts.rda"))

panmut=cosmic.mutation$pancreas

z=values(panmut)$gene.name %in% c("KRAS", "TP53", "CTNNB1", "CDKN2A", "MEN1", "SMAD4",
    "GNAS", "APC", "VHL", "MLL3", "DAXX")

mutblock=subsetByOverlaps(blocky,panmut[z])

pdf(file.path(plotdir, "mutblock.pdf"), width=11, height=8.5)
range.plot(panc, mutblock, grp="anno", logit=F)
anno.region.plot(panc, mutblock, grp="anno", logit=F)
dev.off()


panc=dat[,which(pd$Tissue=="pancreas"|
  grepl("pancreas", pd$Notes)|
  (pd$Phenotype=="normal"&pd$Tissue %in% c("lung", "liver")))]

pancpd=colData(panc)

pancpd$anno=NA
##pancpd$anno[pancpd$Notes=="NET"]="NET"
pancpd$anno[pancpd$Notes=="adenocarcinoma"]="adenocarcinoma"
##pancpd$anno[pancpd$Phenotype=="metastasis"]=pancpd$Notes[pancpd$Phenotype=="metastasis"]
pancpd$anno[pancpd$Phenotype=="metastasis"]="metastasis"
##pancpd$anno[pancpd$Tissue=="lung"&pancpd$Phenotype=="normal"]="lung.normal"
##pancpd$anno[pancpd$Tissue=="liver"&pancpd$Phenotype=="normal"]="liver.normal"
pancpd$anno[pancpd$Phenotype=="normal"&pancpd$Tissue=="pancreas"]="pancreas.normal"
pancpd$anno[pancpd$Phenotype=="adenoma"]="IPMN"

panc=panc[,!is.na(pancpd$anno)]
colData(panc)=pancpd[!is.na(pancpd$anno),]


##Sort blocks and sdmrs on area

osdmr=sdmr[order(-values(sdmr)$area )]

barea=values(blocky)$seg.mean*values(blocky)$num.mark

normy=colData(panc)$anno=="pancreas.normal"
cancy=colData(panc)$anno=="adenocarcinoma"

bdiff=numeric(length(blocky))

for (i in seq(along=blocky)) {
  pprobes=rowData(panc) %in% blocky[i]
  bdiff[i]=median(rowMedians(getBeta(panc[pprobes,cancy]))-
    rowMedians(getBeta(panc[pprobes,normy])))
}

oblocky=blocky[order(-abs(bdiff))]

#oblocky=oblocky[(width(oblocky)>1e4)&(oblocky$num.mark>50)]



##Plot blocks now

pdf(file.path(plotdir, "oblockgg.pdf"), width=11, height=8.5)
range.plot(panc, oblocky[13], grp="anno", logit=F, num.plot=100)
dev.off()


##Plot dmrs now
pdf(file.path(plotdir, "osdmrgg2.pdf"), width=11, height=8.5)

range.plot(panc, osdmr, grp="anno", logit=F, num.plot=100)
dev.off()
