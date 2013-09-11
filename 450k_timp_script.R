codedir="~/Code/timp_illumina"
plotdir="~/Dropbox/Data/Genetics/Infinium/081313_analysis"
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


if (file.exists(file.path(filedir, "thyroid_cancer.rda"))) {
    load(file=file.path(filedir,"thyroid_cancer.rda"))
} else {
    ##Let's get each of the tissue blocks and dmrs
    thyroid.reg=canc.dmrblock(dat, tis="thyroid")
    save(file=file.path(filedir, "thyroid_cancer.rda"), compress="gzip", list=c("thyroid.reg"))
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


compname="thyroid"
sub=dat[,pd$Tissue=="thyroid"]

##Make linkage plot with sig probes
sig.probe.plot(sub, plotdir=plotdir, ccomp="Phenotype", grps=c("normal", "cancer"), compname=compname)

##Plot DMRs
variety.reg.plot(sub, thyroid.reg$dmrs, grp="Phenotype", compname=paste0(compname, "_dmr"), plotdir=plotdir)
##Plot blocks
variety.reg.plot(sub, thyroid.reg$blocks, grp="Phenotype", compname=paste0(compname, "_block"), plotdir=plotdir)


mutdir="~/Dropbox/Data/Genetics/Mutations/092512_dl"
#Mutation tables
load(file.path(mutdir, "muts.rda"))

#panmut=cosmic.mutation$pancreas

#z=values(panmut)$gene.name %in% c("KRAS", "TP53", "CTNNB1", "CDKN2A", "MEN1", "SMAD4",
#    "GNAS", "APC", "VHL", "MLL3", "DAXX")

#mutblock=subsetByOverlaps(pancreas.reg$blocks,panmut[z])

#pdf(file.path(plotdir, "pancmutblock.pdf"), width=11, height=8.5)
#range.plot(panc, mutblock, grp="anno", logit=F)
#anno.region.plot(panc, mutblock, grp="anno", logit=F)
#dev.off()

##colon, breast, lung, thyroid, pancreas

mut.tab=data.frame(row.names=c("colon.mut", "breast.mut", "kidney.mut", "lung.mut", "pancreas.mut", "thyroid.mut"))

mut.tab$colon.block=c(sum(cosmic.mutation$colon %over% colon.reg$blocks),
    sum(cosmic.mutation$breast %over% colon.reg$blocks),
    sum(cosmic.mutation$wilms %over% colon.reg$blocks),
    sum(cosmic.mutation$lung %over% colon.reg$blocks),
    sum(cosmic.mutation$pancreas %over% colon.reg$blocks),
    sum(cosmic.mutation$thyroid %over% colon.reg$blocks))

mut.tab$breast.block=c(sum(cosmic.mutation$colon %over% breast.reg$blocks),
    sum(cosmic.mutation$breast %over% breast.reg$blocks),
    sum(cosmic.mutation$wilms %over% breast.reg$blocks),
    sum(cosmic.mutation$lung %over% breast.reg$blocks),
    sum(cosmic.mutation$pancreas %over% breast.reg$blocks),
    sum(cosmic.mutation$thyroid %over% breast.reg$blocks))

mut.tab$kidney.block=c(sum(cosmic.mutation$colon %over% kidney.reg$blocks),
    sum(cosmic.mutation$breast %over% kidney.reg$blocks),
    sum(cosmic.mutation$wilms %over% kidney.reg$blocks),
    sum(cosmic.mutation$lung %over% kidney.reg$blocks),
    sum(cosmic.mutation$pancreas %over% kidney.reg$blocks),
    sum(cosmic.mutation$thyroid %over% kidney.reg$blocks))

mut.tab$lung.block=c(sum(cosmic.mutation$colon %over% lung.reg$blocks),
    sum(cosmic.mutation$breast %over% lung.reg$blocks),
    sum(cosmic.mutation$wilms %over% lung.reg$blocks),
    sum(cosmic.mutation$lung %over% lung.reg$blocks),
    sum(cosmic.mutation$pancreas %over% lung.reg$blocks),
    sum(cosmic.mutation$thyroid %over% lung.reg$blocks))

mut.tab$pancreas.block=c(sum(cosmic.mutation$colon %over% pancreas.reg$blocks),
    sum(cosmic.mutation$breast %over% pancreas.reg$blocks),
    sum(cosmic.mutation$wilms %over% pancreas.reg$blocks),
    sum(cosmic.mutation$lung %over% pancreas.reg$blocks),
    sum(cosmic.mutation$pancreas %over% pancreas.reg$blocks),
    sum(cosmic.mutation$thyroid %over% pancreas.reg$blocks))

mut.tab$thyroid.block=c(sum(cosmic.mutation$colon %over% thyroid.reg$blocks),
    sum(cosmic.mutation$breast %over% thyroid.reg$blocks),
    sum(cosmic.mutation$wilms %over% thyroid.reg$blocks),
    sum(cosmic.mutation$lung %over% thyroid.reg$blocks),
    sum(cosmic.mutation$pancreas %over% thyroid.reg$blocks),
    sum(cosmic.mutation$thyroid %over% thyroid.reg$blocks))

mut.tab[1,]=mut.tab[1,]/length(cosmic.mutation$colon)
mut.tab[2,]=mut.tab[2,]/length(cosmic.mutation$breast)
mut.tab[3,]=mut.tab[3,]/length(cosmic.mutation$wilms)
mut.tab[4,]=mut.tab[4,]/length(cosmic.mutation$lung)
mut.tab[5,]=mut.tab[5,]/length(cosmic.mutation$pancreas)
mut.tab[6,]=mut.tab[6,]/length(cosmic.mutation$thyroid)


