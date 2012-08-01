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

##Take just thyroid
RGset=read.450k.exp(base=expdatapath, targets=plates[plates$Tissue %in% "thyroid",])
MSet <- preprocessRaw(RGset)

##Minfilocal function
dat=orderByChromosome(MSet)##,everything = TRUE)
dat$pd=pData(MSet)

##QUICK PREPROC
##RafaFunction - winsorization
dat$meth=fixUMoutliers(dat$meth)
dat$unmeth=fixUMoutliers(dat$unmeth)

###LOOK for bad arrays
Umeds=log2(colMedians(dat$unmeth))
Mmeds=log2(colMedians(dat$meth))

##Use plots to look
pdf(file.path(plotdir, "badwolf.pdf"))
plot(Umeds,Mmeds)
dev.off()

##data exploration shows meds should be above 10.5. THIS SHOULD NOT
###BE AUTOMATIC.. EXPLORATION NEEDED TO CHOSE 10.5
##Timp-thyroid choice is ~11
thresh=11

meds=(Umeds+Mmeds)/2
###filter the obviously bad ones
##Did nothing to Amy's Data
dat$meth=dat$meth[,meds>thresh]
dat$unmeth=dat$unmeth[,meds>thresh]
dat$pd=dat$pd[meds>thresh,]


dat=qnorm.subset(dat)

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

###Ready to go...
Y=log2(dat$meth/dat$unmeth)
A=log2(dat$meth+dat$unmeth)

#######################################################
# BLOCKS
#######################################################

MG=500 ##max gap for groups to be like charm.. we choose 500 for now.
MNP=3 ##min number of probes per group
BMNP=15



##Finds blocks of probes
blocks<-collapse450(dat)

island=as.numeric(dat$everything$Relation_to_UCSC_CpG_Island=="Island")

##Concatenate to make single thing
dat$pd$Cat2<-paste(dat$pd$Phenotype,dat$pd$Status,sep="-")
outcome=dat$pd$Cat2

##Timp single descriptor
dat$pd$desc=paste(dat$pd$Tissue, dat$pd$Status, dat$pd$Phenotype, sep="-")

##Do all the comparisons?
#http://www.quickmeme.com/meme/3q8jmd/
comps=combn(unique(outcome),2)

Tab=vector("list",ncol(comps))
Blocktab=vector("list",ncol(comps))


for(h in 1:ncol(comps)){
  
  T1=comps[1,h];T2=comps[2,h]
  compname=paste(T1,T2,sep="V")
  pdfname=paste(compname,"pdf",sep=".")
  csvname1=paste(compname,"csv",sep=".")
  csvname2=paste0(compname,"-blocks.csv")

  keep=outcome%in%c(T1,T2)
  y=Y[,keep]

  ##Drop samples not used in generating the comparison
  tt=factor(outcome[keep],c(T1,T2))

  X=model.matrix(~tt)
  fit=lmFit(y,X)

  ##True blocks
  ##
  ##coef[,2] of fit is (in this case due to fit from linear model) the difference between
  ##the two sample groups
  ##Then take the mean of all differences per cluster
  value=tapply(fit$coef[,2],blocks$pns,mean)

  ##This was backwards for some reason - you want all opensea clusters
  ind=which(blocks$anno$type=="OpenSea")
  ##Create the cna object, giving is locations
  ##Rafa segments chromosome up for areas that have too big a gap can't trust
  ##Blocks.
  cna=CNA(matrix(value[ind],ncol=1),
    chrom=paste(blocks$anno$chr[ind],blocks$anno$blockgroup[ind],sep="_"),
    maploc=blocks$anno$pos[ind],data.type="logratio",presorted=TRUE)
  cbs=segment(cna)
  
  blocktab=cbs$output;names(blocktab)[c(2,3,4)]<-c("pns","start","end")
  blocktab$chr=sub("_.*","",cbs$output$chrom)
  blocktab$pns=sub(".*_","",cbs$output$chrom)
  blocktab$indexStart=cbs$segRows[,1]
  blocktab$indexEnd=cbs$segRows[,2]
  ##area is number of probes times mean difference
  blocktab$area=abs(cbs$output$seg.mean)*cbs$output$num.mark
  blocktab$calc=apply(blocktab[,8:9], 1, function(x) mean((value[ind])[x[1]:x[2]]))
  #*(abs(cbs$output$seg.mean)>0.1)
  blocktab=blocktab[order(-blocktab$area),-1]
  
  
  ##DMRs
  eb=ebayes(fit)
  ##Get the differential methylation(?)
  dm=rowMeans(ilogit(y[,tt==T2]))-rowMeans(ilogit(y[,tt==T1]))
  ##Cluster the probes
  pns=clusterMaker(dat$locs$chr,dat$locs$pos,maxGap=MG)
  ##ss=loessByGroup(fit$coef[,2],pns,dat$locs$pos,se=sqrt(eb$s2.post))
  ## log2 difference
  ss=fit$coef[,2]
  ##Number of probes per cluster
  Ns=tapply(seq(along=pns),pns,length)
  ##Find clusters with at least mininmum probe number
  pnsind=which(pns%in%as.numeric(names(which((Ns>MNP)))))

  ##Just finds areas of consistant difference, based on a simple cutoff of difference
  tab=regionFinder(ss,pns,dat$locs$chr,dat$locs$pos,y=dm,cutoff=0.5,ind=pnsind)
  ##Area is amount of difference*number of probes(clusters in this case)
  tab=tab[tab$area>0.5,]


  ##Take mean value per cluster for each sample
  pp=t(apply(tab[,7:8],1,function(i) colMeans(y[i[1]:i[2],,drop=FALSE])))
  ipp=ilogit(pp)

  
  pdf(file.path(plotdir,pdfname),width=8,height=11)
  plot(rowMedians(pp[,tt==T1]),rowMedians(pp[,tt==T2]),xlim=c(0,3),ylim=c(0,3),
       xlab=T1,ylab=T2,main="log2 Medians comparison")
  abline(0,1)
  plot(rowSds(pp[,tt==T1]),rowSds(pp[,tt==T2]),xlim=c(0,3),ylim=c(0,3),
         xlab=T1,ylab=T2,main="log2 Variance comparison")
  abline(0,1)
  plot(rowMedians(ipp[,tt==T1]),rowMedians(ipp[,tt==T2]),xlim=c(0,1),ylim=c(0,1),
       xlab=T1,ylab=T2,main="Medians comparison")
  abline(0,1)
  plot(rowSds(ipp[,tt==T1]),rowSds(ipp[,tt==T2]),xlim=c(0,.5),ylim=c(0,.5),
         xlab=T1,ylab=T2,main="Variance comparison")
  abline(0,1)
  dev.off()
  
  
  Blocktab[[h]]=blocktab
  Tab[[h]]=tab

  ##Plot top 200 or all dmrs, whichever is less
  M=min(nrow(tab),200)
  ##Plot this far (in bp) on either side
  ADD=2000

  

  pdf(file.path(plotdir, paste0("dmrs", pdfname)), width=8, height=11)
  
  for(i in 1:M){

    ##Take probes within this defined region
    Index=which(dat$locs$chr==tab$chr[i] &
      dat$locs$pos >= tab$start[i] -ADD &
      dat$locs$pos <= tab$end[i] + ADD)
    ##x is genomic position
    x=dat$locs$pos[Index]
    #log2 ratio, just these probes
    yy=y[Index,]
    matplot(jitter(x),ilogit(yy),col=as.fumeric(tt),ylim=c(0,1),
            main=paste(tab$region[i],":",tab$name[i]),
            pch=16,cex=0.75,xlab=paste("location on",tab$chr[i]),ylab="Beta")

    ##Rug with color according to island status
    for(j in seq(along=x)) rug(x[j],col=island[Index][j]+3,lwd=3)

    ##Make a mean line for the regions for each sample (tt)
    tmpIndexes=split(1:ncol(yy),tt)
    for(j in seq(along=tmpIndexes)){
      yyy=rowMeans(yy[,tmpIndexes[[j]]])
      lfit=loess(yyy~x,span=0.75)
      lines(x,ilogit(lfit$fitted),col=j)
    }
    legend("bottomleft",levels(tt),col=1:2,lty=1)
  }
  dev.off()
 
  ##write.csv(tab,paste(outpath,csvname1))
  ##write.csv(blocktab,paste(outpath,csvname2))
}


sel=c("normal", "cancer")
pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()


sel=c("hyperplastic", "cancer")
pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()


sel=c("normal", "hyperplastic")
pdf(file.path(plotdir, paste0(sel[1], sel[2],"mds.pdf")), width=11, height=8.5)
cg.cluster(dat, grps=sel)
dev.off()
