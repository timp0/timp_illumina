preprocessMinfi<- function(object){

  mset <- mapToGenome(RGset, mergeManifest=TRUE)
  
  
}
bumphunter450<-function(mat,design,coef=2,
                        regionName,chr,pos,pns,cutoff=NULL,cutoffQ=0.95,
                        smooth=FALSE,numProbes=5,minNumProbes=5,
                        maxGap=500,minGroupSize=5){
  
  ##maxGap is max gap for groups to be like charm.. we choose 500 for now.
  ##minGroupSize - min number of probes for a group to be considered
  pns=clusterMaker(chr,pos,maxGap=MG)
  Ns=tapply(seq(along=pns),pns,length)
  pnsind=which(pns%in%as.numeric(names(which((Ns>MNP)))))
  
  

collapse450 <- function(chr,pos,relationToIsland,islandName,
                        maxGap=400,blockMaxGap=10^5,
                        locSummary=mean,verbose=TRUE){
  
  ##this function aggregates probes in islands
  ##and the rest are aggregated into clusters
  ##then annotation is created for the new aggregated data
  ##the indexes of the original probes along with the group ids are
  ##returned indexes and pns respectively. note these two are redudant
  ## but they save me work
  
  ##make sure chr is a factor.. Kasper, should this be done from get-go?
  if(is.character(chr)){
    chr=factor(as.character(chr),
      levels=paste("chr",c(1:22,"X","Y"),sep=""))
  }
  if(verbose) cat("Clustering islands and close probes.\n")
  ##block tab has the block group ids
  blocktab=clusterMaker4Blocks(chr=chr,pos,relationToIsland,islandName,
    maxGap=maxGap, FUN=locSummary)
  
  ###split rows by group id
  groupIndexes<-split(seq(along=blocktab$pns),blocktab$pns)

  ##make an object with results
  if(verbose) cat("Computing new annotation.\n")
  res=list(anno=data.frame(id=as.numeric(names(groupIndexes)),
             chr=tapply(as.character(dat$loc$chr),blocktab$pns,first),
             pos=round(tapply(dat$loc$pos,blocktab$pns,mean)),
             type=tapply(as.character(blocktab$type),blocktab$pns,first)),
    indexes=groupIndexes)
  
  if(verbose) cat("Defining blocks.\n")
  ind=res$anno$type=="OpenSea"
  pns<-rep(NA,nrow(res$anno))
  pns[ind] =clusterMaker(res$anno$chr[ind],res$anno$pos[ind],maxGap=blockMaxGap)
  res$anno$blockgroup=pns
  res$pns=blocktab$pns

  return(res)
}

blockFinder <- function(mat,design,chr,pos,relationToIsland,islandName,
                        coef=2,blocks=NULL,
                        maxGap=400,
                        blockMaxGap=10^5,
                        summary=mean,
                        returnBlocks=FALSE,verbose=TRUE,returnGRanges=TRUE){
  require(DNAcopy)
  if(is.null(blocks)){
    blocks<- collapse450(chr,pos,relationToIsland,islandName,
                         maxGap=400,blockMaxGap=10^5,verbose=verbose)
  }
  ind=which(blocks$ann$type=="OpenSea")
  
  fit=lmFit(mat,design)
  value=tapply(fit$coefficients[,coef],blocks$pns,summary)

  cna=CNA(matrix(value[ind],ncol=1),
    chrom=paste(blocks$anno$chr[ind],blocks$anno$blockgroup[ind],sep="_"),
    maploc=blocks$anno$pos[ind],data.type="logratio",presorted=TRUE)
  cbs=segment(cna,verbose=ifelse(verbose,1,0))

  blocktab=cbs$output;names(blocktab)[c(2,3,4)]<-c("pns","start","end")
  blocktab$chr=sub("_.*","",cbs$output$chrom)
  blocktab$pns=sub(".*_","",cbs$output$chrom)
  blocktab$indexStart=cbs$segRows[,1]
  blocktab$indexEnd=cbs$segRows[,2]

  if(returnGRanges){
    tmp<-blocktab[,-c(1,3,4,7)]
    blocktab<-GRanges(seqnames=Rle(blocktab$chr),
                      ranges=IRanges(start=blocktab$start,end=blocktab$end))
    elementMetadata(blocktab)<-tmp
  }
  
  if(returnBlocks) return(list(tab=blocktab,dat=cna,blockIndex=ind,blocks=blocks)) else return(list(tab=blocktab,dat=cna,blockIndex=ind))
}



## blockTable <- function(,fac,blocks=NULL,nulldist=NULL,
##                        ){
##   if(is.null(blocks)) blocks<- collapse450(dat)
##   Index=splitit(blocks$pns)
##   ind=which(blocks$ann$type=="OpenSea")

##   outcome=paste0(tpd$Status,tpd$DignityScore)
##   T1="cancer1";T2="cancer6"
## keep=outcome%in%c(T1,T2)
## y=(Y[,tIndex])[,keep]
## tt=factor(outcome[keep],c(T1,T2))
## sex=factor(tpd$sex[keep],c("M","F"))
## X=model.matrix(~tt+sex)
## fit=lmFit(y,X)
## value=tapply(fit$coef[,2],blocks$pns,mean)
## cna=CNA(matrix(value[ind],ncol=1),
##   chrom=paste(blocks$anno$chr[ind],blocks$anno$blockgroup[ind],sep="_"),
##   maploc=blocks$anno$pos[ind],data.type="logratio",presorted=TRUE)
## cbs=segment(cna)

## obsv <- cbs$output$seg.mean
## obsl <-cut(cbs$output$num.mark,cutoffs,include.lower=TRUE)
## Indexes=split(seq(along=obsv),obsl)
## pvals=vector("numeric",length(obsv))
## for(i in seq(along=Indexes)){
##   tmpind=Indexes[[i]]
##   null<-abs(nulldist[[names(Indexes)[i]]])
##   pvals[tmpind]=sapply(abs(obsv[tmpind]),function(x) mean(x<null))
## }

## blocktab=cbs$output;names(blocktab)[c(2,3,4)]<-c("pns","start","end")
## blocktab$chr=sub("_.*","",cbs$output$chrom)
## blocktab$pns=sub(".*_","",cbs$output$chrom)
## blocktab$indexStart=cbs$segRows[,1]
## blocktab$indexEnd=cbs$segRows[,2]
## blocktab$pval=pvals

## }

fixUMoutliers <- function(mat, K=-3) {
  require(matrixStats)
  TMP <- log2(mat + 0.5)
  MU <- colMedians(TMP)
  SD <- rowMads(t(TMP))         # more to explain
  for(j in 1:ncol(mat)) {
    tmp <- TMP[,j]
    mu <- MU[j]
    sd <- SD[j]
    mat[(tmp - mu)/sd < K,j] <- 2^(K*sd + mu)
  }
  return(mat)
}

qnorm450 <- function(mat,auIndex,xIndex,yIndex,sex){
  mat[auIndex,]=preprocessCore::normalize.quantiles(mat[auIndex,])
  sexIndexes=split(1:ncol(mat),sex)
  for(Index in sexIndexes)
    mat[c(xIndex,yIndex),Index]=preprocessCore::normalize.quantiles(mat[c(xIndex,yIndex),Index])
  mat
}

clusterMaker4Blocks <- function(chr, pos, relationToIsland,
                                islandName, maxGap,FUN=mean) {
  ## get middle position of each island
  tmpName<-paste0(islandName,relationToIsland)
  
  ##these are the "average positions" on shelfs, shores, and islands
  ipos=round(tapply(pos[tmpName!=""], tmpName[tmpName!=""],FUN))
  ## make non-island/shore/shelf positions their own island
  islFactor=ifelse(tmpName!="",tmpName, seq_along(pos))
  ## check which of the new positions correspond to
  ## islands/shores/shelfs, i.e. open sea
  isNotSea=tmpName!=""
  ##assign probes in islands/shores/shelves just one position
  pos[isNotSea]=ipos[islFactor[isNotSea]]
  ## make first pass clusters with new positions
  pns=clusterMaker(chr,pos,maxGap=maxGap)
  ##but islands must be on their own so we change those
  ##beginning and ends of each genomic region type:
  ##island,shore,shelf)
  types=unique(relationToIsland)
  ##take out open sea
  types=types[types!=""]
  for(i in types){
    ind=relationToIsland==i
      StartEnd=abs(diff(c(0,ind)!=0))
    add2pns <- cumsum(StartEnd)
    pns <- pns+ add2pns
  }    
  ## where are the islands
  pns=as.numeric(as.factor(pns))
  
  ##annotation
  type=rep("OpenSea",length(pos))
  type[relationToIsland=="Island"]<-"Island"
  type[grep("Shore",relationToIsland)] <- "Shore"
  type[grep("Shelf",relationToIsland)]<- "Shelf"
  
  return(data.frame(pns=pns,type=type))
}



oldcollapse450 <- function(dat,maxGap=400,blockMaxGap=2*10^5,
                        locSummary=mean,verbose=TRUE){
  
  ##this function aggregates probes in islands
  ##and the rest are aggregated into clusters
  ##then annotation is created for the new aggregated data
  ##the indexes of the original probes along with the group ids are
  ##returned indexes and pns respectively. note these two are redudant
  ## but they save me work
  
  ##THIS NEEDS TO CHANGE
  if(!exists("IlluminaHumanMethylation450kannotationFeatures")){
      if(verbose) cat("Loading probe annotation data.\n")
      load("/Users/ririzarr/tmp/cegs/450/IlluminaHumanMethylation450kannotationFeatures.rda")
    }
  if(verbose) cat("Processing annotation data.\n")
  tmpind=match(rownames(dat$everything),rownames(IlluminaHumanMethylation450kannotationFeatures))
  islandAnn=IlluminaHumanMethylation450kannotationFeatures[tmpind,"HMM_Island"]

  ##make sure chr is a factor.. Kasper, should this be done from get-go?
  if(is.character(dat$loc$chr)){
    dat$loc$chr=factor(as.character(dat$loc$chr),
      levels=paste("chr",c(1:22,"X","Y"),sep=""))
  }
  if(verbose) cat("Clustering islands and close probes.\n")
  ##block tab has the block group ids
  blocktab=clusterMaker4Blocks(chr=dat$loc$chr,
    dat$loc$pos, islandAnn, maxGap=maxGap, FUN=locSummary)

  ###split rows by group id
  groupIndexes<-split(seq(along=blocktab$pns),blocktab$pns)

  ##make an object with results
  if(verbose) cat("Computing new annotation.\n")
  res=list(anno=data.frame(id=as.numeric(names(groupIndexes)),
             chr=tapply(as.character(dat$loc$chr),blocktab$pns,first),
             pos=round(tapply(dat$loc$pos,blocktab$pns,mean)),
             island=tapply(blocktab$isIsland,blocktab$pns,first)),
    indexes=groupIndexes)
  
  if(verbose) cat("Defining blocks.\n")
  ind=!res$anno$island
  pns<-rep(NA,nrow(res$anno))
  pns[ind] =clusterMaker(res$anno$chr[ind],res$anno$pos[ind],maxGap=blockMaxGap)
  res$anno$blockgroup=pns
  res$pns=blocktab$pns

  return(res)
}


oldclusterMaker4Blocks <- function(chr, pos, islandAnn, maxGap,FUN=mean)
  {
    ## get middle position of each island
    islandPos=round(tapply(pos[islandAnn!=""], islandAnn[islandAnn!=""],FUN))
    ## make non-island positions their own island
    islFactor=ifelse(islandAnn!="",islandAnn, seq_along(pos))
    ## check which of the new positions correspond to islands
    isIsland=islandAnn!=""
    ##islands have the same location
    pos[isIsland]=islandPos[islFactor[isIsland]]
    ## make first pass clusters with now positions
    pns=clusterMaker(chr,pos,maxGap=maxGap)
    ##but islands must be on their own so we change those
    ##beginning and ends of islands
    islandStartEnds=abs(diff(c(0,isIsland)!=0))
    ##add to posistions
    add2pns <- cumsum(islandStartEnds)
    ## where are the islands
    pns=as.numeric(as.factor(pns+add2pns))
    return(list(pns=pns,isIsland=isIsland))
  }
}
