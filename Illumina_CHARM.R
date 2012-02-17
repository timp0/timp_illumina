##CHARM plotting functions

plotCHARM <- function(chromy,start,end,class=c(1:10), pts=T) {
  ##class is a name

  ##Set Possible CHARM Colors - don't think I will need more than these - leave them for now
  ##coloring=c("black","coral1","coral4","cadetblue4","cadetblue1","goldenrod1","goldenrod4","seagreen1","seagreen4","plum1","plum4")
  coloring=c("green","purple","goldenrod4","blue","red","goldenrod1","goldenrod4","seagreen1","seagreen4","plum1","plum4")

  ##Find probes which are in my region
  in_region=( (charm_nfo$chr==chromy) & ((charm_nfo$loc>start)&(charm_nfo$loc<end)) )
  ##Which individual value classes match sMM classes being plotted
  classname=colnames(sMM)
  ##Get Correct Samples
  right=sapply(FACS,function(x) any(x==classname[class]))
  ##Assign each type a color
  typecolor=as.character(factor(FACS[right],levels=classname[class],label=coloring[class]))
  #plot sample points
  if (pts) {
    matpoints(charm_nfo$loc[in_region],M[in_region,right],pch=19,col=typecolor, cex=.4)
  }
  #plot smoothed data
  matlines(charm_nfo$loc[in_region],sMM[in_region,class],col=coloring[class],lwd=2, lty=1)

  legend("topright",classname[class],col=coloring[class],lty=1,lwd=2)

  return(in_region)
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



heatCHARMcancer <- function(which,namey="CHARM_heat_cancer.pdf") {

  ch_n=FACS=="colon:normal"
  ch_t=FACS=="colon:tumor"

  norm=apply(which, 1, function(x) apply(M[x,ch_n], 2, median))
  tumor=apply(which, 1, function(x) apply(M[x,ch_t], 2, median))
 
  divider=factor(c(rep(1, times=nrow(norm)), rep(2, times=nrow(tumor)) ))
  labels=divider
  levels(divider)=c("blue", "red")
  levels(labels)=c("colon:normal", "colon:tumor")
  fully=rbind(norm, tumor)
  
  pdf(paste("Movie/", namey, sep=""), width=11, height=8)

  r=seriate(dist(fully), method="OLO", control=list(method="ward"))
  c=seriate(dist(t(fully)), method="OLO", control=list(method="ward"))
  
  col.heat <- heatmap.2(as.matrix(t(fully)), dendrogram="both", trace="none", labCol=labels, ColSideColors=as.character(divider),
                        labRow="", col=brewer.pal(11, "RdYlBu"), hclustfun=hclust.ward,Rowv=as.dendrogram(c[[1]]),
                        Colv=as.dendrogram(r[[1]]))

  col.heat <- heatmap.2(as.matrix(rescale_cpg(t(fully), numdev=1)), dendrogram="both", trace="none", labCol=labels,
                        ColSideColors=as.character(divider),  labRow="", col=brewer.pal(11, "RdYlBu"),
                        hclustfun=hclust.ward, Rowv=as.dendrogram(c[[1]]), Colv=as.dendrogram(r[[1]]) )

  col.heat <- heatmap.2(as.matrix(rescale_cpg(t(fully), numdev=1.5)), dendrogram="both", trace="none", labCol=labels,
                        ColSideColors=as.character(divider),  labRow="", col=brewer.pal(11, "RdYlBu"),
                        hclustfun=hclust.ward, Rowv=as.dendrogram(c[[1]]), Colv=as.dendrogram(r[[1]]), scale="row")
  
  
  dev.off()


}


heatCHARMtissue <- function(which,namey="CHARM_heat_tissue.pdf") {

  ch_c=FACS=="colon:normal"
  ch_l=FACS=="liver:normal"
  ch_s=FACS=="spleen:normal"
  ch_b=FACS=="brain:normal"
  
  colon=apply(which, 1, function(x) apply(M[x,ch_c], 2, median))
  liver=apply(which, 1, function(x) apply(M[x,ch_l], 2, median))
  spleen=apply(which, 1, function(x) apply(M[x,ch_s], 2, median))
  brain=apply(which, 1, function(x) apply(M[x,ch_b], 2, median))

  divider=factor(c(rep(1, times=nrow(colon)), rep(2, times=nrow(liver)), rep(3, times=nrow(spleen)), rep(4, times=nrow(brain)) ))
  labels=divider
  levels(divider)=c("green","purple","goldenrod4","blue")
  levels(labels)=c("colon:normal", "liver:normal", "spleen:normal", "brain:normal")
  fully=rbind(colon,liver, spleen, brain)
  
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


