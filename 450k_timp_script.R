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

if (file.exists(file.path(filedir, "colon_cancer.rda"))) {
    load(file=file.path(filedir,"colon_cancer.rda"))
} else {
    ##Let's get each of the tissue blocks and dmrs
    colon.reg=canc.dmrblock(dat, tis="colon")
    save(file=file.path(filedir, "colon_cancer.rda"), compress="gzip", list=c("colon.reg"))
}

if (file.exists(file.path(filedir, "breast_cancer.rda"))) {
    load(file=file.path(filedir,"breast_cancer.rda"))
} else {
    ##Let's get each of the tissue blocks and dmrs
    breast.reg=canc.dmrblock(dat, tis="breast")
    save(file=file.path(filedir, "breast_cancer.rda"), compress="gzip", list=c("breast.reg"))
}

if (file.exists(file.path(filedir, "pancreas_cancer.rda"))) {
    load(file=file.path(filedir,"pancreas_cancer.rda"))
} else {
    ##Let's get each of the tissue blocks and dmrs
    pancreas.reg=canc.dmrblock(dat, tis="pancreas")
    save(file=file.path(filedir, "pancreas_cancer.rda"), compress="gzip", list=c("pancreas.reg"))
}

if (file.exists(file.path(filedir, "lung_cancer.rda"))) {
    load(file=file.path(filedir,"lung_cancer.rda"))
} else {
    ##Let's get each of the tissue blocks and dmrs
    lung.reg=canc.dmrblock(dat, tis="lung")
    save(file=file.path(filedir, "lung_cancer.rda"), compress="gzip", list=c("lung.reg"))
}

if (file.exists(file.path(filedir, "kidney_cancer.rda"))) {
    load(file=file.path(filedir,"kidney_cancer.rda"))
} else {
    ##Let's get each of the tissue blocks and dmrs
    kidney.reg=canc.dmrblock(dat, tis="kidney")
    save(file=file.path(filedir, "kidney_cancer.rda"), compress="gzip", list=c("kidney.reg"))
}



compname="colon"
##Colon as a start
sub=dat[,(pd$Tissue=="colon"|((pd$Tissue%in%c("liver","lung"))&(pd$Phenotype!="cancer")))]
subpd=colData(sub)

subpd$anno=NA
subpd$anno[subpd$Tissue=="colon" & subpd$Phenotype=="normal"]="colon.normal"
subpd$anno[subpd$Tissue=="lung" & subpd$Phenotype=="normal"]="lung.normal"
subpd$anno[subpd$Tissue=="liver" & subpd$Phenotype=="normal"]="liver.normal"
subpd$anno[subpd$Tissue=="colon" & subpd$Phenotype=="adenoma"]="colon.adenoma"
subpd$anno[subpd$Tissue=="colon" & subpd$Phenotype=="cancer"]="colon.cancer"
subpd$anno[subpd$Notes=="colon.liver.metastasis"]="colon.liver.metastasis"
subpd$anno[subpd$Notes=="colon.lung.metastasis"]="colon.lung.metastasis"

sub=sub[,!is.na(subpd$anno)]
colData(sub)=subpd[!is.na(subpd$anno),]

##Make linkage plot with sig probes
sig.probe.plot(sub, plotdir=plotdir, ccomp="anno", grps=c("colon.normal", "colon.cancer"), compname=compname)

##Plot DMRs
variety.reg.plot(sub, colon.reg$dmrs, grp="anno", compname=paste0(compname, "_dmr")
                 ,plotdir=plotdir)
##Plot blocks
variety.reg.plot(sub, colon.reg$blocks, grp="anno", compname=paste0(compname, "_block")
                 ,plotdir=plotdir)


compname="simpcolon"
##Simple Colon
sub=dat[,(pd$Tissue=="colon")]
subpd=colData(sub)

subpd$anno=NA
subpd$anno[subpd$Tissue=="colon" & subpd$Phenotype=="normal"]="colon.normal"
subpd$anno[subpd$Tissue=="colon" & subpd$Phenotype=="adenoma"]="colon.adenoma"
subpd$anno[subpd$Tissue=="colon" & subpd$Phenotype=="cancer"]="colon.cancer"
subpd$anno[grepl("metastasis", subpd$Notes)]="colon.metastasis"

sub=sub[,!is.na(subpd$anno)]
colData(sub)=subpd[!is.na(subpd$anno),]

##Make linkage plot with sig probes
sig.probe.plot(sub, plotdir=plotdir, ccomp="anno",
               grps=c("colon.normal", "colon.cancer"), compname=compname)


##Plot DMRs
variety.reg.plot(sub, colon.reg$dmrs, grp="anno", compname=paste0(compname, "_dmr")
                 ,plotdir=plotdir)
##Plot blocks
variety.reg.plot(sub, colon.reg$blocks, grp="anno", compname=paste0(compname, "_block")
                 ,plotdir=plotdir)



compname="pancreas"
sub=dat[,which(pd$Tissue=="pancreas"|grepl("pancreas", pd$Notes)|
    (pd$Phenotype=="normal"&pd$Tissue %in% c("lung", "liver")))]

subpd=colData(sub)

subpd$anno=NA
subpd$anno[subpd$Notes=="NET"]="NET"
subpd$anno[subpd$Notes=="adenocarcinoma"]="pancreas.cancer"
subpd$anno[subpd$Phenotype=="metastasis"]=subpd$Notes[subpd$Phenotype=="metastasis"]
subpd$anno[subpd$Tissue=="lung"&subpd$Phenotype=="normal"]="lung.normal"
subpd$anno[subpd$Tissue=="liver"&subpd$Phenotype=="normal"]="liver.normal"
subpd$anno[subpd$Phenotype=="normal"&subpd$Tissue=="pancreas"]="pancreas.normal"
subpd$anno[subpd$Phenotype=="adenoma"]="IPMN"

sub=sub[,!is.na(subpd$anno)]
colData(sub)=subpd[!is.na(subpd$anno),]

##Make linkage plot with sig probes
sig.probe.plot(sub, plotdir=plotdir, ccomp="anno", grps=c("pancreas.normal", "pancreas.cancer"), compname=compname)



##Plot DMRs
variety.reg.plot(sub, pancreas.reg$dmrs, grp="anno", compname=paste0(compname, "_dmr"), plotdir=plotdir)
##Plot blocks
variety.reg.plot(sub, pancreas.reg$blocks, grp="anno", compname=paste0(compname, "_block"), plotdir=plotdir)


compname="pancreassimple"
sub=dat[,which(pd$Tissue=="pancreas"|grepl("pancreas", pd$Notes))]

subpd=colData(sub)

subpd$anno=NA
subpd$anno[subpd$Notes=="adenocarcinoma"]="pancreas.cancer"
subpd$anno[subpd$Phenotype=="metastasis"]="metastasis"
subpd$anno[subpd$Phenotype=="normal"&subpd$Tissue=="pancreas"]="pancreas.normal"
subpd$anno[subpd$Phenotype=="adenoma"]="IPMN"


sub=sub[,!is.na(subpd$anno)]
colData(sub)=subpd[!is.na(subpd$anno),]

##Make linkage plot with sig probes
sig.probe.plot(sub, plotdir=plotdir, ccomp="anno", grps=c("pancreas.normal", "pancreas.cancer"), compname=compname)


##Plot DMRs
variety.reg.plot(sub, pancreas.reg$dmrs, grp="anno", compname=paste0(compname, "_dmr"), plotdir=plotdir)
##Plot blocks
variety.reg.plot(sub, pancreas.reg$blocks, grp="anno", compname=paste0(compname, "_block"), plotdir=plotdir)



compname="breast"
sub=dat[,pd$Tissue=="breast"]

##Make linkage plot with sig probes
sig.probe.plot(sub, plotdir=plotdir, ccomp="Phenotype", grps=c("normal", "cancer"), compname=compname)

##Plot DMRs
variety.reg.plot(sub, breast.reg$dmrs, grp="Phenotype", compname=paste0(compname, "_dmr"), plotdir=plotdir)
##Plot blocks
variety.reg.plot(sub, breast.reg$blocks, grp="Phenotype", compname=paste0(compname, "_block"), plotdir=plotdir)



compname="lung"
sub=dat[,pd$Tissue=="lung"]

##Make linkage plot with sig probes
sig.probe.plot(sub, plotdir=plotdir, ccomp="Phenotype", grps=c("normal", "cancer"), compname=compname)

##Plot DMRs
variety.reg.plot(sub, lung.reg$dmrs, grp="Phenotype", compname=paste0(compname, "_dmr"), plotdir=plotdir)
##Plot blocks
variety.reg.plot(sub, lung.reg$blocks, grp="Phenotype", compname=paste0(compname, "_block"), plotdir=plotdir)


compname="kidney"
sub=dat[,pd$Tissue=="kidney"]

##Make linkage plot with sig probes
sig.probe.plot(sub, plotdir=plotdir, ccomp="Phenotype", grps=c("normal", "cancer"), compname=compname)

##Plot DMRs
variety.reg.plot(sub, kidney.reg$dmrs, grp="Phenotype", compname=paste0(compname, "_dmr"), plotdir=plotdir)
##Plot blocks
variety.reg.plot(sub, kidney.reg$blocks, grp="Phenotype", compname=paste0(compname, "_block"), plotdir=plotdir)




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
