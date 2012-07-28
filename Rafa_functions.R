##Convenience functions:
##Get first index of array
first<- function(x,k=1) x[k]
ilogit=function(x) 1/(1+exp(-x))
as.fumeric<- function(x,...) as.numeric(as.factor(x,...))


qnorm450 <- function(mat,auIndex,xIndex,yIndex,sex){
  mat[auIndex,]=preprocessCore::normalize.quantiles(mat[auIndex,])
  sexIndexes=split(1:ncol(mat),sex)
  for(Index in sexIndexes)
    mat[c(xIndex,yIndex),Index]=preprocessCore::normalize.quantiles(mat[c(xIndex,yIndex),Index])
  mat
}

clusterMaker <- function(chr,pos,order.it=TRUE,maxGap=300){
  nonaIndex=which(!is.na(chr) & !is.na(pos))
  ##split into list on chromosomes
  Indexes=split(nonaIndex,chr[nonaIndex])
  clusterIDs=rep(NA,length(chr))
  LAST=0
  ##Loop through Chromosomes
  for(i in seq(along=Indexes)){
    Index=Indexes[[i]]
    x=pos[Index]

    if(order.it){ Index=Index[order(x)];x=pos[Index] }

    y=as.numeric(diff(x)>maxGap)
    z=cumsum(c(1,y))
    clusterIDs[Index]=z+LAST
    LAST=max(z)+LAST
  }
  return(clusterIDs)
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
  ##It has N_Shore and Shelf?
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


collapse450 <- function(dat,maxGap=400,blockMaxGap=10^5,
                        locSummary=mean,verbose=TRUE){
  
  ##this function aggregates probes in islands
  ##and the rest are aggregated into clusters
  ##then annotation is created for the new aggregated data
  ##the indexes of the original probes along with the group ids are
  ##returned indexes and pns respectively. note these two are redudant
  ## but they save me work
  
  ##make sure chr is a factor.. Kasper, should this be done from get-go?
  if(is.character(dat$loc$chr)){
    dat$loc$chr=factor(as.character(dat$loc$chr),
      levels=paste("chr",c(1:22,"X","Y"),sep=""))
  }
  if(verbose) cat("Clustering islands and close probes.\n")
  ##block tab has the block group ids
  blocktab=clusterMaker4Blocks(chr=dat$loc$chr,
    dat$loc$pos,
    dat$everything$Relation_to_UCSC_CpG_Island,
    dat$everything$UCSC_CpG_Islands_Name,
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
  pns[ind]=clusterMaker(res$anno$chr[ind],res$anno$pos[ind],maxGap=blockMaxGap)
  res$anno$blockgroup=pns
  res$pns=blocktab$pns

  return(res)
}



fixUMoutliers <- function(mat, K=-3) {
  ##Relatively efficient winsorization function

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

boyorgirl <- function(A,xIndex,yIndex,plot=FALSE,id=1:ncol(A)){
  ##A is average of log int
  x=Biobase::rowMedians(t(A[xIndex,]),na.rm=TRUE)
  y=Biobase::rowMedians(t(A[yIndex,]),na.rm=TRUE)
  k=kmeans(y-x,c(min(y-x),max(y-x)))
  sex=factor(ifelse(k$cluster==which.min(k$centers),"F","M"),levels=c("M","F"))
  if(plot){
    plot(x,y,type="n")
    text(x,y,id,col=as.numeric(sex))
    legend("bottomleft",c("M","F"),col=c(1,2),pch=1)
  }
  return(sex)
}

regionFinder<-function(x,regionNames,chr,position,y=x,
                       summary=mean,ind=seq(along=x),order=TRUE,oneTable=TRUE,
                       ...){

  Indexes=getSegments(x[ind],regionNames[ind],...)
  
  res=vector("list",2)
  for(i in 1:2){
    res[[i]]=data.frame(chr=sapply(Indexes[[i]],function(Index) chr[ind[Index[1]]]),
         start=sapply(Indexes[[i]],function(Index) min(position[ind[Index]])),
         end=sapply(Indexes[[i]],function(Index) max(position[ind[Index]])),
         value=sapply(Indexes[[i]],function(Index) summary(y[ind[Index]])),
         area=sapply(Indexes[[i]],function(Index) abs(sum(y[ind[Index]]))),
         pns=sapply(Indexes[[i]],function(Index) regionNames[ind[Index]][1]),
         indexStart=sapply(Indexes[[i]],function(Index) min(ind[Index])),
         indexEnd=sapply(Indexes[[i]],function(Index) max(ind[Index])))
    res[[i]]$L=res[[i]]$indexEnd-res[[i]]$indexStart+1
  }

  names(res)=c("up","dn")
  if(order & !oneTable){
    if(nrow(res$up)>0) res$up=res$up[order(-res$up$area),]
    if(nrow(res$dn)>0) res$dn=res$dn[order(-res$dn$area),]
  }
  if(oneTable){
    res=rbind(res$up,res$dn)
    if(order & nrow(res)>0) res=res[order(-res$area),]
  }
  return(res)

}


getSegments <- function(x,factor,cutoff=quantile(abs(x),0.99),verbose=TRUE){

  Indexes=split(seq(along=x),factor)
  regionID=vector("numeric",length(x))
  LAST = 0
  
  segmentation = vector("numeric", length(x))
  type = vector("numeric", length(x))

  for (i in seq(along = Indexes)) {
    if (verbose) if (i%%1000 == 0) cat(".")
    Index = Indexes[[i]]
    y = x[Index]
    z = sign(y) * as.numeric(abs(y) > cutoff)
    w = cumsum(c(1, diff(z) != 0)) + LAST
    segmentation[Index] = w
    type[Index] = z
    LAST = max(w)
  }
  ##add a vector of the pns
  res=list(upIndex=split(which(type>0),segmentation[type>0]),
    dnIndex=split(which(type<0),segmentation[type<0]),
    zeroIndex=split(which(type==0),segmentation[type==0]))
  names(res[[1]])<-NULL
  names(res[[2]])<-NULL
  names(res[[3]])<-NULL
  return(res)
}
