##Stuff from Amy, but altered by me
##To find blocks and regions
if (!exists("codedir")) {
  codedir=getwd()
}

plotdir="~/Dropbox/Data/Genetics/Infinium/082712_analysis"
filedir="~/Dropbox/Data/Genetics/Infinium/082712_analysis"

source(file.path(codedir, "450k_general_init.R"))


if (file.exists(file.path(filedir, "ast.rda"))) {
  load(file.path(filedir, "ast.rda"))
} else {
  target=read.450k.sheet(base=file.path(expdatapath,"IL066"),pattern=".csv$",verbose=TRUE)

  ##Now we subset Micahel's data out of Yun's.
  Index=which(target[,"Experimenter"]=="Michael") 
  target=target[Index,]

  dat=dat.preload(plates=target,plotdir=plotdir,expdatapath=expdatapath, plotter=T,sex=F)

  save(list="dat", file=file.path(filedir, "ast.rda"))
}

##TAKE OUT SNPS- Chris
x=load(file.path("~/Dropbox/Data/Genetics/Infinium/071312_analysis", "snps_chris.rda"))
snps1=get(x)
snps1=snps1[match(rownames(dat$meth),snps1$IlmnID),]
keepIndex=which(snps1$SBEsnp_RefSeqID=="FALSE"&snps1$CGsnp_RefSeqID=="FALSE")

##Ok - because dat is meth, unmeth, locs, everything (arrays of probes) get rid of annotation at same time
for(i in 1:4) dat[[i]]=dat[[i]][keepIndex,]

##Rafa code (according to Michael)

pd=dat$pd
y=log2(dat$meth/dat$unmeth)
tt=factor(pd$Phenotype,c("Asthmatic","Non-asthmatic"))
sex=factor(pd$Sex,c("m","f"))
X=model.matrix(~tt+sex)
fit=lmFit(y,X)
eb=ebayes(fit)
MG=500 ##max gap for groups to be like charm.. we choose 500 for now.
MNP=3 ##min number of probes per group
pns=clusterMaker(dat$locs$chr,dat$locs$pos,maxGap=MG)
Ns=tapply(seq(along=pns),pns,length)
pnsind=which(pns%in%as.numeric(names(which((Ns>MNP)))))
dm=rowMeans(ilogit(y[,tt=="Asthmatic"]))-rowMeans(ilogit(y[,tt=="Non-asthmatic"]))

###I have no idea what cutoff=5 does.
rafa.tab=regionFinder(eb$t[,2],pns,dat$locs$chr,dat$locs$pos,y=dm,cutoff=5,ind=pnsind)

##Winston code

##dmr.find <- function(dat, ccomp="Phenotype", grps=c("normal", "cancer"), MG=500, MNP=3, cutoff=0.5,
sel=c("Asthmatic", "Non-asthmatic")
pdf(file.path(plotdir, paste0(sel[1], sel[2],"dmr.pdf")), width=11, height=8.5)
winston.tab=dmr.find(dat, grps=sel, cutoff=5)

dev.off()
