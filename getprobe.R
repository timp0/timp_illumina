between <- function(num,lower,upper,strict=FALSE){
    if(strict)  return(num>lower & num<upper)
    if(!strict) return(num>=lower & num<=upper)
}

library(gdata)
x <- read.xls("~/feinberg/winston/Supplementary_Data_4.xls",sheet=2, pattern="distToCGI", header=TRUE)
load("~wtimp/cancer_dmr/From_Rafa/otherstuff.rda")
load("~wtimp/cancer_dmr/From_Rafa/sMM.rda")

j <- which(names(INDEX)=="colon:normal")
k <- which(names(INDEX)=="colon:tumor")
#(names(INDEX) are column names for sMM)
sD <- sMM[,j]-sMM[,k]
sT <- sD/sqrt(ssigmas[,j]+ssigmas[,k])

#The genomic positions in pos are only meaningful when seen together with
#pns, which tells what chromosome they're on.  They index positions in the
#chromosome, not the whole genome.
pnschr <- sapply(strsplit(pns,":"),function(v){return(v[1])})
stopifnot(identical(levels(as.factor(pns)),unique(pns)))
pbs   <- vector("list",length=nrow(x))
index <- vector("integer",length=nrow(x))
for(i in 1:nrow(x)){
    thischr <- which(pnschr==x[i,"chr"])
    pos0    <- pos[thischr]
    ins     <- which(between(pos0,x[i,2],x[i,3]))
    ans <- thischr[ins]
    pbs[[i]] <- ans
    
    pn <- unique(pns[ans])
    if(length(pn)==2 & pn[1]=="chr11:1900000-2900000") pn <- pn[2]
      #there are 2 of them and both are within the smaller subregion.
    index[i] <- which(unique(pns)==pn)
}
#pbs is the elements of pos(=>rows of sD and sT) that are between each DMR, going
#down the rows of x
#Check that index's chromosomes are right:
indexchr <- sapply(strsplit(unique(pns)[index],":"),function(v){return(v[1])})
stopifnot(all(indexchr==x[,1]))

maxesT <- sapply(pbs,function(b){ which.max(abs(sT[b])) })
maxesD <- sapply(pbs,function(b){ which.max(abs(sD[b])) })
maxprobeT <- numeric(length=length(pbs))
maxprobeD <- numeric(length=length(pbs))
for(i in 1:length(pbs)){
    maxprobeT[i] <- pbs[[i]][maxesT[i]]
    maxprobeD[i] <- pbs[[i]][maxesD[i]]
}

coordsT <- pos[maxprobeT] #the 50-mer probes start at these positions.
coordsD <- pos[maxprobeD] #the 50-mer probes start at these positions.

#A check, to make sure these positions are inside the DMRs:
for(k in 1:length(coordsT)){
    stopifnot(between(coordsT[k],x[k,2],x[k,3]))
    stopifnot(between(coordsD[k],x[k,2],x[k,3]))
}

x$index <- index
x$maxprobestartT <- coordsT
x$maxprobestartD <- coordsD
write.csv(x,file="~/feinberg/winston/SupData4.csv",row.names=FALSE,quote=FALSE)
 
#############################################################################
#Now for the plots:
#This code comes from results.R: 
