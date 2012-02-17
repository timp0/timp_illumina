##CHARM plotting functions

plotCHARM <- function(chromy,start,end,class=c(1:10)) {
  ##class is a name

  ##Set Possible CHARM Colors - don't think I will need more than these - leave them for now
  coloring=c("black","coral1","coral4","cadetblue4","cadetblue1","goldenrod1","goldenrod4","seagreen1","seagreen4","plum1","plum4")

  ##Find probes which are in my region
  in_region=( (charm_nfo$chr==chromy) & ((charm_nfo$loc>start)&(charm_nfo$loc<end)) )
  ##Which individual value classes match sMM classes being plotted
  classname=colnames(sMM)
  ##Get Correct Samples
  right=sapply(FACS,function(x) any(x==classname[class]))
  ##Assign each type a color
  typecolor=as.character(factor(FACS[right],levels=classname[class],label=coloring[class]))
  #plot sample points
  matpoints(charm_nfo$loc[in_region],M[in_region,right],pch=1,col=typecolor)
  #plot smoothed data
  matlines(charm_nfo$loc[in_region],sMM[in_region,class],col=coloring[class],lwd=2, lty=1)

  legend("topright",classname[class],col=coloring[class],lty=1,lwd=2)

}

plotfullCHARM <- function(probes,genes,avg_data,sampy,outname,charm=c(1:10),controls=TRUE) {

  ##
  ##Start Plot
  pdf(outname,width=11,height=8)
  par(mfrow=c(3,1))
  layout(rbind(1:5,rep(6,5),rep(7,5)),heights=c(0.5,0.3,0.2))
  par(mar=c(2,1,0,1),mgp=c(1.5,.5,0),oma=c(0,1,2,1)) 

  
  
  ##Loop through all regions
  for (i in unique(probes$Region)) {
    ##Find probes which are in the same region - this was previously defined by perl in the probes$Region command
    not_far<-probes[(probes$Region==i),]
    num_close<-nrow(not_far)
    ##cat("Region: ",i,"\n")
    ##Set which chromosome and coordinates we are using for the region
    chromy <- paste("chr",not_far$Chromosome[1],sep="")
    index=(min(not_far$Start_loc)-1000):(max(not_far$Finish_loc)+1000)
    final_index=length(index)
    ##Plot actual samples
    for (k in 1:5){
      ##Init sample plot
      if (k==1){
        plot(0,0,ylim=c(0,1.1),xlim=c(index[900],index[final_index-900]),ylab="Methylation",xlab="",type="n",axes=F)
        box()
        axis(2, at=seq(0,1,.2))
      } else {
        plot(0,0,ylim=c(0,1.1),xlim=c(index[900],index[final_index-900]),ylab="",xlab="",type="n",xaxt="n",yaxt="n")
      }
      if (any(sampy$Class==(k*2))) {
        ##Plot actual samples
        plotsamples(class_stuff[(k*2):(k*2+1),],probes,avg_data,sampy,not_far,i,index[1])
        ##Plot label on axis as a tick on the bottom
        axis(side=1,at=not_far$Start_loc,labels=not_far$Probe_ID, cex.axis=0.8)
        ##axis(side=1)
      }
    }
    ##Plot title to graph
    mtext(paste("ID:",i,"--",as.character(chromy),":",index[1],"-",index[final_index],sep=""),cex=2,side=3,outer=TRUE)
    
    ##Make initial plot
    plot(0,0,ylim=c(-0.5,3),xlim=c(index[1],index[final_index]),ylab="Methylation",xlab="",type="n",axes=F)
    box()
    axis(2)

    ##Plot CHARM Data
    plotCHARM(chromy,index[1],index[final_index],charm)
    ##Plot label on axis as a tick on the bottom
    axis(side=1,at=not_far$Start_loc,labels=not_far$Probe_ID, cex.axis=0.8)

    ##Plot RefSeq Gene regions
    plotgenes(genes,chromy, index[1], index[final_index])
  }
  dev.off()
}

CHARMmatch <- function (chrom, loc, charm_nfo, num=5) {


  
  rel=array(0,c(length(chrom), num))
  
  for (i in 1:length(chrom)) {
    
    ##Matching chromosome check - give a low score if matches for sorting
    chr=charm_nfo$chr!=paste("chr",chrom[i],sep="")
    ##distance to illumina probe i
    distbp=abs(charm_nfo$loc-loc[i])
    rel[i,]=order(chr,distbp)[1:num]    
    ##sort and take top num
    print(i)
  }

  return(rel)
  
}



heatCHARM <- function(namey="CHARM_heat.pdf") {

  ch_n=FACS=="colon:normal"
  ch_t=FACS=="colon:tumor"

  norm=apply(CHIL, 1, function(x) apply(M[x,ch_n], 2, median))
  tumor=apply(CHIL, 1, function(x) apply(M[x,ch_t], 2, median))
 
  divider=factor(c(rep(1, times=nrow(norm)), rep(2, times=nrow(tumor)) ))
  labels=divider
  levels(divider)=rainbow(2)
  levels(labels)=c("colon:normal", "colon:tumor")
  fully=rbind(norm, tumor)
  
  pdf(paste("Movie/", namey, sep=""), width=11, height=8)

  r=seriate(dist(fully), method="OLO", control=list(method="ward"))
  c=seriate(dist(t(fully)), method="OLO", control=list(method="ward"))
  
  col.heat <- heatmap.2(as.matrix(t(fully)), dendrogram="both", trace="none", labCol=labels, ColSideColors=as.character(divider), labRow="", col=brewer.pal(11, "RdYlBu"), hclustfun=hclust.ward,
                        Rowv=as.dendrogram(c[[1]]), Colv=as.dendrogram(r[[1]]) )

  col.heat <- heatmap.2(as.matrix(rescale_cpg(t(fully), numdev=1)), dendrogram="both", trace="none", labCol=labels,
                        ColSideColors=as.character(divider),  labRow="", col=brewer.pal(11, "RdYlBu"),
                        hclustfun=hclust.ward, Rowv=as.dendrogram(c[[1]]), Colv=as.dendrogram(r[[1]]) )

  col.heat <- heatmap.2(as.matrix(rescale_cpg(t(fully), numdev=1.5)), dendrogram="both", trace="none", labCol=labels,
                        ColSideColors=as.character(divider),  labRow="", col=brewer.pal(11, "RdYlBu"),
                        hclustfun=hclust.ward, Rowv=as.dendrogram(c[[1]]), Colv=as.dendrogram(r[[1]]), scale="row")
  
  
  dev.off()


}

plotPlatform <- function(avg, sampy,CHIL, namey="CHIL1.pdf") {
  #Ok - plot avg of 5 closest CHARM values in the area vs the avg probe value

  #First let's calculate the difference value for each probe in Illumina
  norm=(sampy$Progression==0)&(sampy$Class==3)
  tum=(sampy$Progression>2)&(sampy$Class==3)

  normvals=rowMeans(avg[,norm])
  tumvals=rowMeans(avg[,tum])
  ildiff=tumvals-normvals

  probec=ifelse(data$probes$UCSC_Dist_to_Island>0,"red", "green")
  probec[data$probes$UCSC_Dist_to_Island>2000]="blue"
  
  charmdiff=apply(CHIL,1,function(x) median(sMM[x,5]-sMM[x,4]))
  
  
  pdf(file.path("Movie", namey))
  plot(ildiff,charmdiff, col=probec, pch=19, ylab="CHARM difference", xlab="Illumina difference")
  legend("topright", c("Islands", "Shores", "Far"), col=c("red", "green", "blue"), pch=19)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
  
  dev.off() 
  
}




##CHARM probe loc load
if (!exists("charm_nfo")) {
  load("CHARM_nfo.rda")
}

##Load CHARM Data
if (!exists("M")) {
  load("CHARM_data.rda")
}

if (!exists("CHIL")) {
  load("CHIL.rda")
}

source("~/Work/Analysis/Illumina/Illumina_cluster1b.R")


##plotfullCHARM(data$probes,genes,avg_data,sampy,"Movie/colon_charm.pdf",c(4:5))
