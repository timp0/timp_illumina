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
##pancpd$anno[pancpd$Phenotype=="metastasis"]=pancpd$Notes[pancpd$Phenotype=="metastasis"]
pancpd$anno[pancpd$Phenotype=="metastasis"]="metastasis"
##pancpd$anno[pancpd$Tissue=="lung"&pancpd$Phenotype=="normal"]="lung.normal"
##pancpd$anno[pancpd$Tissue=="liver"&pancpd$Phenotype=="normal"]="liver.normal"
pancpd$anno[pancpd$Phenotype=="normal"&pancpd$Tissue=="pancreas"]="pancreas.normal"
pancpd$anno[pancpd$Phenotype=="adenoma"]="IPMN"

panc=panc[,!is.na(pancpd$anno)]
colData(panc)=pancpd[!is.na(pancpd$anno),]



design=model.matrix(~factor(subpd$Phenotype)+factor(subpd$predictedSex))

res=bumphunter(sub,design,B=100,smooth=FALSE)
sres=bumphunter(sub, design, B=100, smooth=T)

cobj=cpgCollapse(sub)
blocks=blockFinder(cobj$object,design,B=100)

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

pdf(file.path(plotdir, paste0("dmrggplot.pdf")), width=11, height=8.5)
range.plot(panc, dmr, grp="anno", logit=F)
dev.off()
 
##Gviz plots!
pdf(file.path(plotdir, paste0("dmrgviz.pdf")), width=11, height=8.5)
anno.region.plot(panc, dmr, grp="anno", logit=F)
dev.off()

pdf(file.path(plotdir, paste0("sdmrggplot.pdf")), width=11, height=8.5)
range.plot(panc, sdmr, grp="anno", logit=F)
dev.off()
 
##Gviz plots!
pdf(file.path(plotdir, paste0("sdmrgviz.pdf")), width=11, height=8.5)
anno.region.plot(panc, sdmr, grp="anno", logit=F)
dev.off()

blocky=blocks$tab

##Plot blocks now
pdf(file.path(plotdir, "blockgg.pdf"), width=11, height=8.5)
range.plot(panc, blocky, grp="anno", logit=F)
dev.off()

pdf(file.path(plotdir, "blockgviz.pdf"), width=11, height=8.5)
anno.region.plot(panc, blocky, grp="anno", logit=F)
dev.off()

mutdir="~/Dropbox/Data/Genetics/Mutations/092512_dl"

#Mutation tables
load(file.path(mutdir, "muts.rda"))

panmut=cosmic.mutation$pancreas

z=match(c("KRAS", "TP53", "CTNNB1", "CDKN2A", "MEN1", "SMAD4",
  "GNAS", "APC", "VHL", "MLL3", "DAXX"), values(panmut)$gene.name)

mutblock=subsetByOverlaps(blocky,panmut[z])

pdf(file.path(plotdir, "mutblock.pdf"), width=11, height=8.5)
range.plot(panc, mutblock, grp="anno", logit=F)
anno.region.plot(panc, mutblock, grp="anno", logit=F)
dev.off()
