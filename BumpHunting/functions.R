as.fumeric<- function(x,...) as.numeric(as.factor(x,...))
first<- function(x,k=1) x[k]
###CHECK THIS ~pmurakam/feinberg/CharmFiles/functions.R
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

splitit <- function(x) split(seq(along=x),x)

myplclust <- function( hclust, lab=hclust$labels, lab.col=rep(1,length(hclust$labels)), hang=0.1, ... )
{
  ## modifiction of plclust for plotting hclust objects *in colour*!
  ## Copyright Eva KF Chan 2009
  ## Arguments:
  ##    hclust:    hclust object
  ##    lab:        a character vector of labels of the leaves of the tree
  ##    lab.col:    colour for the labels; NA=default device foreground colour
  ##    hang:     as in hclust & plclust
  ## Side effect:
  ##    A display of hierarchical cluster with coloured leaf labels.
  y <- rep(hclust$height,2)
  x <- as.numeric(hclust$merge)
  y <- y[which(x<0)]
  x <- x[which(x<0)]
  x <- abs(x)
  y <- y[order(x)]
  x <- x[order(x)]
  plot( hclust, labels=FALSE, hang=hang, ... )
    text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), labels=lab[hclust$order], col=lab.col[hclust$order], srt=90, adj=c(1,0.5), xpd=NA, ... )
  
}

dogo <- function(names,universe,species="human"){
  if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Hs.egREFSEQ2EG
  } else  {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA))
  
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))
  
  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = "BP", pvalueCutoff = 0.01, conditional = TRUE,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)
}




affyGO <- function(names,universe=NULL,platform="hgu133plus2",species="human",pvalueCutoff=1){
  library(paste(platform,"db",sep="."),character.only=TRUE)
  gomap=get(paste(platform,"ENTREZID",sep=""))
  symmap=get(paste(platform,"SYMBOL",sep=""))
  if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
  } else  {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]

  if(is.null(universe)) universe=ls(gomap)
  Universe=unlist(mget(universe,gomap,ifnotfound = NA))
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))
  
  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = "BP", pvalueCutoff = pvalueCutoff,
                conditional = TRUE,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs= sapply(tmp1,function(y) paste(unlist(mget(names(x)[x%in%y],symmap,ifnotfound=NA)),collapse=";"))
  return(tab)
}

myfilter2 <- function(x,filter,...){
###unlike myfilter, myfilter2 returns NAs in the edges.
  L=dim(x)[1]
  if(L>length(filter)) res=filter(x,filter,...) else{res = t(colMeans(x) %*% t(rep(1,dim(x)[1])))}

  return(res)

}


getDesc<-function(x){
 ##get gene description
  require(org.Hs.eg.db)
  genenames=sapply(mget(as.character(x),org.Hs.egREFSEQ2EG,ifnotfound = NA),function(x) x[1])
  genenames[!is.na(genenames)]= sapply(mget(genenames[!is.na(genenames)],org.Hs.egGENENAME,ifnotfound=NA),function(x) x[1])
  genenames
}


midpoints <- function(x) (x[1:(length(x)-1)]+x[2:length(x)])/2

pointsplot <- function(l,jitter=TRUE,factor=1,col=rep(1,length(l)),pch=21,...){
  if(!is.list(l)) stop("l must be a list\n")
  y=unlist(l)
  x=rep(seq(along=l),sapply(l,length))
  col=rep(col,sapply(l,length))
  if(jitter) x=jitter(x,factor)
  if(pch!=21) plot(x,y,col=col,pch=pch,...) else plot(x,y,pch=pch,bg=col,...)
}

tplot <- function(x,xlab="",ylab="",...){
  plot(as.numeric(names(x)),x,xlab=xlab,ylab=ylab,...)
}

scatterBox <- function(x,y,cuts=100,...) boxplot(split(y,cut(x,unique(quantile(x,seq(0,1,length=cuts+1))),include.lowest=TRUE)),...)


maplot <- function(x,y,...){ A=(x+y)/2;M=y-x;plot(A,M,...)}
                             
                             
