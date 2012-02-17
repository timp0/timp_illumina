##Put us in the right dir for output
setwd("~/Data/Infinium/122711_analysis/")

##Experimental file location
expdatapath="/thumper2/feinbergLab/core/arrays/illumina/"

##Current 450k package name
##Kasper says ignore warnings about "contains no R code"
##This is the local(expanded) version of minfi
library(minfi)
library(minfiLocal)

##This will list commands
##ls("package:minfi")


##Read in plates
plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))

##Read in data
RGset=read.450k.exp(base=expdatapath, targets=plates)


##Get methyl and unmethyl signal

##ID lung and liver samples
lsq.samp=(((pData(RGset)$Tissue=="lung")|(pData(RGset)$Tissue=="liver")) &
          (pData(RGset)$Status=="normal"))


##Get out just pancreas samples and mets
pancreas.samp=pData(RGset)$Tissue=="pancreas" ##Alternatively pData(RGset) gets out the plates variable again
pan.met.samp=grepl(pattern="pancreas",pData(RGset)$Notes)&pData(RGset)$Phenotype=="metastasis"

pancreas.data=RGset[,(pancreas.samp|pan.met.samp|lsq.samp)]

##Get out just colon samples and mets
colon.samp=pData(RGset)$Tissue=="colon" ##Alternatively pData(RGset) gets out the plates variable again
col.met.samp=grepl(pattern="colon",pData(RGset)$Notes)&pData(RGset)$Phenotype=="metastasis"

colon.data=RGset[,(colon.samp|col.met.samp|lsq.samp)]


pan.MSet=preprocessMinfi(pancreas.data)
col.MSet=preprocessMinfi(colon.data)

load("~/Data/Infinium/121311_analysis/probe_obj_final.rda")

##Add something indexing gprobes vs actual data order
values(gprobes)$minfi.idx=match(values(gprobes)$name, rownames(getM(pan.MSet[,1])))

##Keep just probes with no SNPs
good.probes=values(gprobes)$minfi.idx[ (!values(gprobes)$sbe.snp.boo)&
  (!values(gprobes)$boo.snps)]


pan.label=ifelse(pData(pancreas.data)$Phenotype=="metastasis",
  pData(pancreas.data)$Notes, paste(pData(pancreas.data)$Tissue,
                                    pData(pancreas.data)$Phenotype, sep="."))

col.label=ifelse(pData(colon.data)$Phenotype=="metastasis",
  pData(colon.data)$Notes, paste(pData(colon.data)$Tissue,
                                    pData(colon.data)$Phenotype, sep="."))



pdf("Plots/met1.pdf")


mdsPlot(pan.MSet[good.probes,], numPositions=1000,
        sampGroups=pan.label)
mdsPlot(col.MSet[good.probes,], numPositions=1000,
        sampGroups=col.label)

dev.off()                               



##PCA

p=prcomp(t(data$fqbeta[top25,]))

##MDS

z=cmdscale( dist( t( data$fqbeta[top25,]) ), k=2)

norms=data$fsamp$Progression==0
carcy=data$fsamp$Progression>2
ady=data$fsamp$Progression>0&data$fsamp$Progressoin<3

##Now thyroid progression plot
colly=c("green", "orange", "blue", "red","brown")

type_col=factor(data$fsamp$Class)

types=levels(type_col)

levels(type_col)=colly

pnorms=p$x[norms,]
pn_med=apply(pnorms,2,median)
pn_mad=apply(pnorms,2,mad)

nprof=list(xcoor=rep(0,5), ycoor=rep(0,5), xrad=rep(0,5), yrad=rep(0,5), col=rep(as.character("black"),5))

for (i in 1:5) {
  type_norms=data$fsamp$Class==types[i]&data$fsamp$Progression==0
  ztn=z[type_norms,]
  nprof$xcoor[i]=median(ztn[,1])
  nprof$ycoor[i]=median(ztn[,2])
  nprof$xrad[i]=mad(ztn[,1])*3
  nprof$yrad[i]=mad(ztn[,2])*3
  nprof$col[i]=colly[i]
}

pca_range_x=range(p$x[,1])
pca_range_y=range(p$x[,2])

mds_range_x=range(z[,1])
mds_range_y=range(z[,2])

pdf("Movie/norms_mds.pdf")

plot(z[norms,1], z[norms,2], bg=as.character(type_col[norms]), pch=21, xlim=mds_range_x, ylim=mds_range_y)

draw.ellipse(x=nprof$xcoor, y=nprof$ycoor, a=nprof$xrad, b=nprof$yrad, lty=2, lwd=2,border=nprof$col)

legend("topright", c("Breast", "Colon", "Lung", "Wilms", "Thyroid"), col=as.character(levels(type_col)), pch=16)

dev.off()

pdf("Movie/carc_mds.pdf")

plot(z[carcy,1], z[carcy,2], bg=as.character(type_col[norms]), pch=21, xlim=mds_range_x, ylim=mds_range_y)

draw.ellipse(x=nprof$xcoor, y=nprof$ycoor, a=nprof$xrad, b=nprof$yrad, lty=2, lwd=2,border=nprof$col)

legend("topright", c("Breast", "Colon", "Lung", "Wilms", "Thyroid"), col=as.character(levels(type_col)), pch=16)

dev.off()

#Colon adenoma/tumor sample stuff
colons=data$fsamp$Class==3

colony=list(fqbeta=data$fqbeta[,colons], fsamp=data$fsamp[colons,])

##MDS

z=cmdscale( dist( t( colony$fqbeta[top25,]) ), k=2)

norms=colony$fsamp$Progression==0
carcy=colony$fsamp$Progression>2
ady=colony$fsamp$Progression>0&colony$fsamp$Progression<3

##Now progression plot
pcolly=c("blue", "green", "red")

prog_col=factor(ifelse(colony$fsamp$Progression>0,1,0)+ifelse(colony$fsamp$Progression>2,1,0))

levels(prog_col)=pcolly

pnorms=p$x[norms,]
pn_med=apply(pnorms,2,median)
pn_mad=apply(pnorms,2,mad)

ztn=z[norms,]
c_x=median(ztn[,1])
c_y=median(ztn[,2])
c_a=mad(ztn[,1])*3
c_b=mad(ztn[,2])*3


mds_range_x=range(z[,1])
mds_range_y=range(z[,2])

pdf("Movie/colons_prog_mds.pdf")

plot(z[,1], z[,2], bg=as.character(prog_col), pch=21, xlim=mds_range_x, ylim=mds_range_y)

draw.ellipse(x=rep(c_x,2), y=rep(c_y,2), a=rep(c_a,2), b=rep(c_b,2), lty=2, lwd=2, border="green")

legend("topleft", c("Normal", "Adenoma", "Carcinoma"), col=as.character(levels(prog_col)), pch=16)

dev.off()

