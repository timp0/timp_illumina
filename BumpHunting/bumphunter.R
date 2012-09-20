regionMatch <- function(x, subject, verbose=TRUE) {
    require(IRanges)

    ret = matrix(NA,nrow=length(x),ncol=7)
    colnames(ret) = c("dist","matchIndex","type","amountOverlap","insideDist","size1","size2")
    sp1 = split(1:length(x), as.character(seqnames(x)))
    sp2 = split(1:length(subject), as.character(seqnames(subject)))
    for(i in names(sp1)) {
        if(verbose) cat(i," ")
        inds1 = sp1[[i]]
        if(i %in% names(sp2)) {
          inds2 = sp2[[i]]
          ret[inds1,"matchIndex"] = inds2[
               nearest(ranges(x[inds1,]),ranges(subject[inds2,]))]
        } else ret[inds1,"matchIndex"] = NA
      }
    if(verbose) cat("\n")
    y = subject[ret[,"matchIndex"],]
    ret = data.frame(ret, stringsAsFactors=FALSE)
    
    ret$type[(start(x)>start(y) & end(x)<=end(y)) |
             (start(x)>=start(y) & end(x)<end(y))] = "inside"
    ret$type[start(x)<=start(y) & end(x)>=end(y)] = "cover"
    ret$type[start(x)>end(y)] = "disjointR"
    ret$type[end(x)<start(y)] = "disjointL"
    ret$type[is.na(ret$matchIndex)] = "disjoint"
    ret$type[(start(x)>start(y) & start(x)<=end(y)) & end(x)>end(y)] = "overlapR"    
    ret$type[start(x)<start(y) & (end(x)>=start(y) & end(x)<end(y))] = "overlapL"

    ret$dist = 0
    ret$dist[ret$type=="disjoint"]  = NA
    ret$dist[ret$type=="disjointR"] = end(y)[ret$type=="disjointR"] - start(x)[ret$type=="disjointR"]
    ret$dist[ret$type=="disjointL"] = start(y)[ret$type=="disjointL"] - end(x)[ret$type=="disjointL"]
    ret$amountOverlap[ret$type=="overlapR"] = -1*(end(y)[ret$type=="overlapR"]-start(x)[ret$type=="overlapR"]+1)
    ret$amountOverlap[ret$type=="overlapL"] = end(x)[ret$type=="overlapL"]-start(y)[ret$type=="overlapL"]+1
    ret$type[ret$type%in%c("disjointR","disjointL")] = "disjoint"
    ret$type[ret$type%in%c("overlapR","overlapL")] = "overlap"

    ## insideDist column:
    insideIndex = ret$type=="inside" #no missing ret$type at this point
    tmp0 = cbind(end(x)[insideIndex]  -end(y)[insideIndex],
                 start(x)[insideIndex]-start(y)[insideIndex])
    tmp = apply(abs(tmp0),1,which.min)
    tmpinsidedist = tmp0[,1]
    tmpIndex = tmp==2
    tmpinsidedist[tmpIndex] = tmp0[tmpIndex,2]
    ret$insideDist[insideIndex] = tmpinsidedist

    ## size1 and size2 columns:
    ret$size1 = end(x) -start(x) +1
    ret$size2 = end(y)-start(y)+1
    
    ret
}

stringToLocations <- function(x){
  tmp <- sapply(strsplit(x,":"),function(y) c(y[1],y[2]))
  tmp2 <- sapply(strsplit(tmp[2,],"-"),
                 function(y) c(as.numeric(y[1]),as.numeric(y[2])))
  return(data.frame(chr=I(tmp[1,]),
                    start=tmp2[1,],
                    end=tmp2[2,]))
}

##THIS NEEDS TO BE EDITED TO USE GENOMICRANGES
## pointMatch = function(object1,object2,col1=2,col2=2,verbose=TRUE){
##   Indexes1=split(seq(along=as.character(object1$chr)),as.character(object1$chr))
##   Indexes2=split(seq(along=as.character(object2$chr)),as.character(object2$chr))
##   res=matrix(NA,nrow=nrow(object1),2)
##   for(i in seq(along=Indexes1)){
##     CHR=names(Indexes1)[i]
##     if(verbose) cat(CHR,",")
##     Index1=Indexes1[[CHR]]
##     Index2=Indexes2[[CHR]]
##     if(length(Index1)>0 & length(Index2)>0){
##       tmpmat1=object1[Index1,col1]
##       tmpmat2=object2[Index2,col2]
##       tmp=fuzzy.match(tmpmat1,tmpmat2)
##       tmp[,2] = Index2[tmp[,2]]
	  
##       res[Index1,]=tmp    
##      }
##   }
##   colnames(res)<-c("dist","matchIndex")
##   return(res)
## }

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


###COMING SOON: HARRIS' VERSION

##HARRIS IS CHANGING THIS
nearestgene <- function(object,up=5000,down=5000,species="human",
                        build="hg18",
                        datadir="/dexter/disk4/gcode/hji/cisreg/genomes/",
                        wd="./cisgenome",
                        binary="/home/bst/faculty/hji/projects/cisgenome_project/bin/refgene_getnearestgene",
                        inputFileName="tmp.cod",
                        outputFileName="tmp.txt",
                        r=0){##is the type of distance

  

  tmpscipen=.Options$scipen
  options(scipen=100)
  
  if(!file.exists(wd)) {
    cat("\nCreating directory:",wd,"\n")
    dir.create(wd)
  }
  
  inputfile=file.path(wd,inputFileName)
  outputfile=file.path(wd,outputFileName)
  write.table(data.frame(object$chr,object$start,object$end,rep("+",nrow(object))),file=inputfile,sep="\t",quote=FALSE,col.names=FALSE)

  thecall=paste(binary,"-d",file.path(datadir,species,build,"annotation/refFlat_sorted.txt"),"-dt 1 -s",species,"-i",inputfile,"-o",outputfile,"-r", r,"-up",up,"-down",down)
  system(thecall)
  res=read.delim(outputfile,as.is=TRUE,check.names=FALSE,comment.char="#",header=FALSE,col.names=c("seq_id","chr","start","end","strand","name","annotation","num","genestrand","TSS","TSE","CSS","CSE","ne","exonStarts","exonEnds"),fill=NA)
  res[res=="---"]<-NA
  res[res==""] <- NA
  options(scipen=tmpscipen)
  res
}


matchGenes <- function(object,promoterDist=2500,verbose=TRUE,species="human",
                       build="hg18",
                       r=0,up=50000000,down=50000000){
  if(verbose) cat("Matching regions to genes.\n")
  genes=nearestgene(object,up=up,down=down,species=species,build=build,r=r)
  type=rep("",nrow(genes))
  subtype=rep("",nrow(genes))
  ctype=rep("",nrow(genes))
  dist=rep(0,nrow(genes))
  insidedistance<-rep(NA,nrow(genes))
  exonnumber<-rep(NA,nrow(genes))
  nexons<-rep(NA,nrow(genes))
  
  geneL=rep(0,nrow(genes))
  codingL=rep(0,nrow(genes))
  subdist=rep(0,nrow(genes))
  if(verbose) cat("Annotating results")
  for(i in 1:nrow(genes)){
    if(verbose & i%%1000==0)  cat(".")

    if(!is.na(genes[i,10])){
    TS = genes[i,10]
    TE = genes[i,11]
    geneL[i]=TE-TS
    CS = genes[i,12]
    CE = genes[i,13]
    codingL[i]=CE-CS
    ES = as.numeric(strsplit(genes[i,15],",")[[1]])
    EE = as.numeric(strsplit(genes[i,16],",")[[1]])
    Exons= cbind(ES,EE)
    nexons[i]=nrow(Exons)
    Strand= ifelse(genes[i,9]=="+",1,-1)
    S = genes[i,3]
    E = genes[i,4]
    
    type[i]=""
    if(genes[i,3] <= TS & genes[i,4] >= TE){
      type[i]="covers"
    } else{
      if(genes[i,3] < TS & genes[i,4] < TS){
        if(Strand==1){
          type[i]="upstream" 
          dist[i]=TS-E
        } else{
          type[i]="downstream"
          dist[i]=TE-E
        }
      }
      if(genes[i,3] > TE & genes[i,4] > TE){
        if(Strand==-1){
          type[i]="upstream"
          dist[i]=S-TE
        }  else{
          type[i]="downstream"
          dist[i]=S-TS
        }
      }
      if (genes[i,3]>= TS & genes[i,4] <= TE){
        type[i]="inside"
        if(Strand==-1) dist[i]=TE-E  else dist[i]=S-TS
      }
      
      ##overlaps
      if(type[i]==""){
        if(genes[i,3] < TS & genes[i,4] <=TE){
          ##OVERLAP FRONT
          if(Strand==1) type[i]="overlaps 5'" else{
            type[i]="overlaps 3'"
            dist[i]=TE-E
          }
          S=TS
        }
        if (genes[i,3] >= TS & genes[i,4] > TE){
          ##OVERAL BACK
          if(Strand==-1) type[i]="overlaps 5'" else{
            type[i]="overlaps 3'"
            dist[i]=S-TS
          }
          E=TE
        }
      }
    }
    m1=NA;m2=NA
    if(S >= TS & E <= TE){
      
      ##INSIDE
      tmp1=fuzzy.match2(S,Exons)
      tmp2=fuzzy.match2(E,Exons)
      
      m1=tmp1[1,1]
      m2=tmp2[1,1]
      exon1=tmp1[1,2]
      exon2=tmp2[1,2]
      m1m2Index=which.min(abs(c(m1,m2)))
      
      if(exon1==exon2 & m1==0 & m2==0){
        subtype[i]="inside exon"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      
      if( (sign(m1)==sign(m2) & (m1!=0 & m2!=0) & exon1==exon2) |
         (sign(m1)==-1 & sign(m2)==1 & exon2-exon1==1) ){
        subtype[i]="inside intron"

        insidedistance[i]=c(m1,m2)[m1m2Index]
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
        
      }
      if( (exon2-exon1 > 1) |
         ( (exon2-exon1 == 1) & 
          ((sign(m1)==-1 & sign(m2)==-1) |
           (sign(m1)==1 & sign(m2)==-1) |
           (sign(m1)==1 & sign(m2)==1) |
           (sign(m1)==0 & sign(m2)==-1)|
           (sign(m1)==1 & sign(m2)==0))) |
         (exon2==exon1 & sign(m1)==1 & sign(m2)==-1)){
        subtype[i]="covers exon(s)"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      
      if( (exon2-exon1 == 1 & sign(m1)==-1 & sign(m2)==0) |
         (exon2==exon1 & sign(m1)==1 & sign(m2)==0)){
        if(Strand==1) subtype[i]="overlaps exon upstream" else subtype[i]="overlaps exon downstream"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      if( (exon2-exon1 == 1 & sign(m1)==0 & sign(m2)==1) |
         (exon2==exon1 & sign(m1)==0 & sign(m2)==-1)){
        if(Strand==-1) subtype[i]="overlaps exon upstream" else subtype[i]="overlaps exon downstream"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      if( exon2-exon1 == 1 & sign(m1)==0 & sign(m2)==0){
        subtype[i]="overlaps two exons"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }

      if(Strand!=1){
        insidedistance[i]= -insidedistance[i]
        exonnumber[i] = nrow(Exons) - exonnumber[i] + 1
      }
      
      ctype[i]="inside transcription region"
      
      if(S<CS & E<CS){
        if(Strand==1) ctype[i]="5' UTR" else ctype[i]="3'UTR"
      }
      
      if(S>CE & E>CE){
        if(Strand==-1) ctype[i]="5' UTR" else ctype[i]="3'UTR"
      }
      if(S<CS & E>CE){
        ctype[i]="covers transcription region"
      }
      if(S<CS & E>CS){
        if(Strand==1) ctype[i]="overlaps 5' UTR" else ctype[i]="overlaps 3'UTR"
      }
      if(S<CE & E>CE){
        if(Strand==-1) ctype[i]="overlaps 5' UTR" else ctype[i]="overlaps 3'UTR"
      }
    }
    if(FALSE){##graphical check
      plot(0,0,ylim=c(0,1),xlim=range(genes[i,c(3:4,10:11)]),xlab=paste("inside distance=",insidedistance[i],m1,m2,exonnumber[i]))
      polygon(c(TS,TE,TE,TS),c(0,0,0.5,0.5),density=0,col=2)
      polygon(c(CS,CE,CE,CS),c(0.1,0.1,0.4,0.4),density=0,col=3)
      abline(h=0.25,lwd=2)
      apply(Exons,1,function(x) polygon(c(x[1],x[2],x[2],x[1]),
                                        c(0.2,0.2,0.3,0.3),col=4))
      polygon(as.vector(genes[i,c(3,4,4,3)]),c(0.75,0.75,0.85,0.85),col=5)
      lines(c(TS,TS+1000),c(0.65,0.65),lwd=3)
      title(paste(i,Strand,type[i],subtype[i],ctype[i],dist[i],sep=":"))

      
    }
  }
  }
  if(verbose) cat("Done.\n")
  type[dist<=promoterDist & type=="upstream"] <- "promoter"
  type[dist<=promoterDist & type=="downstream"] <- "close to 3'"

  description=type
  tmpIndex=which(description=="inside")
  description[tmpIndex] <- subtype[tmpIndex]
  return(data.frame(name=I(genes[,6]),
                    annotation=I(genes[,7]),
                    description=factor(description,levels=c("upstream","promoter","overlaps 5'","inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons","overlaps 3'","close to 3'","downstream","covers")),
                    region=factor(type,levels=c("upstream","promoter","overlaps 5'","inside","overlaps 3'","close to 3'","downstream","covers")),
                    distance=dist,
                    subregion=factor(subtype,levels=c("inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons")),
                    insidedistance=insidedistance,
                    exonnumber=exonnumber,
                    nexons=nexons,
                    UTR=factor(ctype,levels=c("inside transcription region","5' UTR","overlaps 5' UTR","3'UTR","overlaps 3'UTR","covers transcription region")),
                    strand=genes[,9],
                    geneL=geneL,
                    codingL=codingL))
}

clusterMaker <- function(chr,pos,order.it=TRUE,maxGap=300){
  nonaIndex=which(!is.na(chr) & !is.na(pos))
  Indexes=split(nonaIndex,chr[nonaIndex])
  clusterIDs=rep(NA,length(chr))
  LAST=0
  for(i in seq(along=Indexes)){
    Index=Indexes[[i]]
    x=pos[Index]

    if(order.it){ Index=Index[order(x)];x=pos[Index] }

    y=as.numeric(diff(x)>maxGap)
    z=cumsum(c(1,y))
    clusterIDs[Index]=z+LAST
    LAST=max(z)+LAST
  }
  clusterIDs
}  


runmedByGroup <- function(y, chr=NULL, pos=NULL, group, k=5,
                          endrule="constant",
                          verbose=TRUE) {
  ##we dont use chr and pos byt
  ##Depends on limma
  ##bpSpan is in basepairs
  ##assumed y are ordered
  if(is.null(dim(y))) y=matrix(y,ncol=1) ##need to change this to be

  Indexes <- split(seq(along=group),group)
  groupL <- sapply(Indexes,length)
  spans <- rep(NA,length(Indexes))##using spans for consistency
  loessed <- rep(TRUE,nrow(y))
  for(i in seq(along=Indexes)) {
    if(verbose) if(i %% 1e4 == 0) cat(".")
    Index = Indexes[[i]]
    spans[i]<-k 
    if(groupL[i] >= k) {
      for(j in 1:ncol(y)){
        y[Index] <- runmed(y[Index],k=k,endrule=endrule)
      }
    } else{ y[Index] <- NA;smoothed[Index]=FALSE}
  }
  return(list(fitted=y,spans=spans,groupL=groupL))
}

##you can pass cutoff through the ...
regionFinder<-function(x,chr,pos,group=NULL,y=x,
                       summary=mean,ind=seq(along=x),
                       order=TRUE,oneTable=TRUE,maxGap=300,
                       ...){
  if(any(is.na(x[ind]))){
    warning("NAs found and removed. ind changed.")
    ind=intersect(which(!is.na(x)),ind)
  } 
  if(is.null(group)) group=clusterMaker(chr,pos,maxGap=maxGap)
  Indexes=getSegments(x[ind],group[ind],...)
  groupN<-table(group)[group]##very inneficient.. but we want length
  ##of original interval to know how much we had to choose from
  
  res=vector("list",2)
  for(i in 1:2){
    res[[i]]=data.frame(chr=sapply(Indexes[[i]],function(Index)chr[ind[Index[1]]]),
         start=sapply(Indexes[[i]],function(Index) min(pos[ind[Index]])),
         end=sapply(Indexes[[i]],function(Index) max(pos[ind[Index]])),
         value=sapply(Indexes[[i]],function(Index) summary(y[ind[Index]])),
         area=sapply(Indexes[[i]],function(Index) abs(sum(y[ind[Index]]))),
         group=sapply(Indexes[[i]],function(Index) group[ind[Index]][1]),
           indexStart=sapply(Indexes[[i]],function(Index) min(ind[Index])),
         indexEnd=sapply(Indexes[[i]],function(Index) max(ind[Index])))
    res[[i]]$L=res[[i]]$indexEnd-res[[i]]$indexStart+1
    res[[i]]$groupL=sapply(Indexes[[i]],function(Index) groupN[ind[Index]][1])
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


.lm <- function(mat,design,coef,B=NULL){
  ##From HCB does fast permuations
  v=design[,coef]
  A=design[,-coef,drop=FALSE]
  qa=qr(A)
  S=diag(nrow(A)) - tcrossprod(qr.Q(qa))
  
  vv <- if (is.null(B)) matrix(v,nc=1) else replicate(B, sample(v))
  (mat %*% crossprod(S,vv)) / diag(crossprod(vv,S) %*% vv)
  
  ##need to add residuals etc... to get sd for limma
}

loessByGroup <- function(y, x, group, weights=rep(1,length(y)),
                         bpSpan = 1000, 
                         minNum=7,minInSpan=5,maxSpan=1,
                         verbose=TRUE) {
  ##Depends on limma
  ##bpSpan is in basepairs
  ##assumed x are ordered
  ##if y is vector change to matrix
  if(is.null(dim(y))) y=matrix(y,ncol=1) ##need to change this to be
  ##more robust
  Indexes <- split(seq(along=group),group)
  groupL <- sapply(Indexes,length)
  spans <- rep(NA,length(Indexes))
  smoothed <- rep(TRUE,nrow(y))
  for(i in seq(along=Indexes)) {
    if(verbose) if(i %% 1e4 == 0) cat(".")
    Index = Indexes[[i]]
    if(groupL[i] >= minNum) {
      span = bpSpan/median(diff(x[Index]))/length(Index) # make a span
      if(span > maxSpan) span <- maxSpan
      spans[i]<-span
      if(span*length(Index)>minInSpan){
        ##this can be parallelized
        for(j in 1:ncol(y)){
          y[Index,j] <- limma::loessFit(y[Index,j],x[Index],span =
                                        span,weights = weights[Index])$fitted
        }
      } else{ y[Index,]<-NA;spans[i]<-NA;smoothed[Index]<-FALSE}
    } else{ y[Index,] <- NA;smoothed[Index]<-FALSE}
  }
  return(list(fitted=y,smoothed=smoothed,spans=spans,groupL=groupL))
}

###THIS IS THE WRAPPER... smoothing not permitted yet
bumphunter <- function(object, design, 
                       chr,pos,group=NULL,
                       coef=2,cutoff=NULL,cutoffQ=0.975,
                       maxGap=500,smooth=TRUE,
                       B=100,verbose=TRUE,
                       smoothFunction=loessByGroup,
                       ...){
  ##B is the number of random samples we take

  cat("Computing coefficients.\n")
  if(is.null(group)) group=clusterMaker(chr,pos,maxGap=maxGap)

  beta<-.lm(object,design,coef) 
  if(smooth){
    beta<-smoothFunction(beta,pos,group,verbose=verbose,...) ##weights come latter
    Index=which(beta$smoothed)
    beta=beta$fitted
  }
  
  if(is.null(cutoff)) cutoff=quantile(abs(beta),cutoffQ)

  cat("Finding regions.\n")
  tab=regionFinder(beta,chr,pos,group,cutoff=cutoff,ind=Index,verbose=FALSE)
  
  cat("Performing",B,"permutations.\n")
  fcs<- .lm(object,design,coef,B)
  if(smooth)
    fcs<- smoothFunction(fcs,pos,group,verbose=verbose)$fitted

  cat("Computing regions for each permutation.\n")
  L<-vector("list",B)
  V<-vector("list",B)
  A<-vector("list",B)
  for(i in 1:B){
    cat(".")
    nulltab=regionFinder(fcs[,i],chr,pos,group,cutoff=cutoff,ind=Index,
      verbose=FALSE)
    L[[i]]<-nulltab$L
    V[[i]]<-nulltab$value
    A[[i]]<-nulltab$area
  }
  tots=sapply(seq(along=V),function(i)
    apply(cbind(tab$L,abs(tab$value)),1,function(x)
          sum(L[[i]]>=x[1] & abs(V[[i]])>=x[2])))
  rate1=rowMeans(tots>0)
  rate2=rowSums(tots)/B
  tots2=sapply(seq(along=A),function(i)
    sapply(tab$area,function(x) sum(A[[i]]>=x[1])))
  rate3=rowMeans(tots2>0)
  rate4=rowSums(tots2)/B

  tab$rate=rate2
  tab$perRunRate=rate1
  tab$areaRate=rate4
  tab$areaPerRunRate=rate3

  tab=tab[order(tab$rate,tab$area),]
  return(list(table=tab,fitted=beta,null=list(value=V,length=L)))
}
  
  
