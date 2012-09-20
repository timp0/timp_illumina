###This next commang gets i from SGE
##it is the number of the comparison we will make
##the list of comparisons is in compTable
library(limma)
library(GenomicFeatures)
library(matrixStats)
library(DNAcopy)
source("../functions.R")
source("../bumphunter.R")
source("../450/450functions.R")

cat("Loading and preping data.")
load("rdas/data.rda")
cat(".")

###Just thyroid to start
pd2=read.csv("FeinbergCohortPath.csv",as.is=TRUE)
ind=which(pd$Tissue=="thyroid" & !pd$Sample.ID%in%c("Thyroid_Normal_14","Thyroid_Normal_6"))
tpd=pd[ind,]
theind=match(tpd$Sample.ID,pd2$minfi.id)
tpd=cbind(tpd,pd2[theind,])
tmpind=which(is.na(tpd$LcT) & tpd$X.Tumor=="")
tpd=tpd[tmpind,]
tIndex=ind[tmpind]

##harthle - grep("HC|HA",tpd$Sample.ID)]
#######################################################
### DMR finding
#######################################################
keep=tpd$DignityScore%in%c(0,6)
type=as.fumeric(tpd$DignityScore[keep])
sex=factor(tpd$sex[keep],c("M","F"))
X=model.matrix(~type+sex)
res=bumphunter(Y[,tIndex[keep]],X,chr,pos,group,cutoff=1,smooth=FALSE)



##for testing
###design=X;object=Y[,tIndex[keep]];group=NULL;cutoff=1;maxGap=500;smooth=TRUE;B=10;coef=2


fit=lmFit(Y[,tIndex[keep]],X)
eb=ebayes(fit)

pns=clusterMaker(chr,pos,maxGap=750)
tmp<-loessByGroup(fit$coef[,2],pos,pns)
#tmp2<-loessByGroup(fit$coef[,2],pos,pns,weight=1/sqrt(eb$s2.post))
#tmp3<-runmedByGroup(fit$coef[,2],pns,k=5)
tmp4=regionFinder(tmp$fitted,chr,pos,pns)

  
mypar(2,1)
for(i in tab$group){
  ind=which(pns==i)
  if(!is.na(tmp$spans[i])){
    plot(pos[ind],fit$coef[ind,2],main=tmp$span[i],ylim=c(-1,1))
    lines(pos[ind],tmp$fitted[ind],col=2)
    lines(pos[ind],tmp2$fitted[ind],col=3)
    lines(pos[ind],tmp3$fitted[ind],col=4)
    plot(pos[ind],sqrt(eb$s2.post)[ind])
  }
}
#######################################################
## START BLOCK FINDING
#######################################################

load("rdas/blocks.rda")

keep=tpd$DignityScore%in%c(1,3,4,5,6)
type=as.fumeric(tpd$DignityScore[keep])
sex=factor(tpd$sex[keep],c("M","F"))
X=model.matrix(~type+sex)

b=blockFinder(Y[,tIndex[keep]],X,dat$locs$chr,dat$locs$pos,
  dat$everything$Relation_to_UCSC_CpG_Island,
  dat$everything$UCSC_CpG_Islands_Name,
  coef=2,
  returnBlocks=TRUE,blocks=blocks)



##GET NULL DISTRIBUTION
B<-100
cat("Performing",B,"permutations")
L <- vector("list",B)
V <- vector("list",B)
for(j in 1:B){
  cat(j,"")
  sX=model.matrix(~sample(type)+sex)
  nb=blockFinder(Y[,keep],sX,dat$locs$chr,dat$locs$pos,
    dat$everything$Relation_to_UCSC_CpG_Island,
    dat$everything$UCSC_CpG_Islands_Name,
    blocks=blocks,
    returnBlocks=TRUE)
  
  L[[j]]<-elementMetadata(nb$tab)$num.mark
  V[[j]]<-elementMetadata(nb$tab)$seg.mean
  cat(length(L[[j]]),"\n")
}

l=unlist(L)
cutoffs <- c(0,unique(quantile(l,seq(0,1,len=round(length(l)/10000)))))
cutoffs[length(cutoffs)]<-Inf
l=cut(l,cutoffs,include.lower=TRUE)
v=unlist(V)
nulldist=split(v,l)
obsv <- elementMetadata(b$tab)$seg.mean
obsl <-cut(elementMetadata(b$tab)$num.mark,cutoffs,include.lower=TRUE)
Indexes=split(seq(along=obsv),obsl)
pvals=vector("numeric",length(obsv))
for(j in seq(along=Indexes)){
  tmpind=Indexes[[j]]
  null<-abs(nulldist[[names(Indexes)[j]]])
  pvals[tmpind]=sapply(abs(obsv[tmpind]),function(x) mean(x<null))
}

elementMetadata(b$tab)$pvals=pvals
usedIndex<-keepIndex[keep]
save(b,usedIndex,keepIndex,pd,file="rdas/progression-blocks-thyroid.rda")
system("rm rdas/tmp-progression-blocks-thyroid.rda")
     

     
