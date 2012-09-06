##This regionFinder function is identical to normal except that the methylation differences given by y (default in 450k - dm, as y=dm in regionFinder) are cubed before being added together. This gives more weight to probes with larger average methylation differences. This gives less weight to DMRs that are lots of probes with little methylation differences.

regionFinder<-function(x,regionNames,chr,position,y=x,
                       summary=mean,ind=seq(along=x),order=TRUE,oneTable=TRUE,
                       ...){
  
  Indexes=getSegments(x[ind],regionNames[ind],...)
  
  res=vector("list",2)
  ##Get up, then down from Indexes(getSegments) which finds areas of up or down comparison
  
  for(i in 1:2){
    res[[i]]=data.frame(chr=sapply(Indexes[[i]],function(Index) chr[ind[Index[1]]]),
                        start=sapply(Indexes[[i]],function(Index) min(position[ind[Index]])),
                        end=sapply(Indexes[[i]],function(Index) max(position[ind[Index]])),
                        value=sapply(Indexes[[i]],function(Index) summary(y[ind[Index]])),
                        area=sapply(Indexes[[i]],function(Index) abs(sum((y[ind[Index]])^3))),
#                        area=sapply(Indexes[[i]],function(Index) abs(sum((y[ind[Index]])))),                        
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
