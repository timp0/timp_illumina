##This code is used to setup clustering for Illumina samples

hclust.ward <- function(x) {
  hc = hclust(x, method="ward")
}

fclust.ward <- function(x) {
  hc= flashClust(x, method="ward")
}

sigmoid_color<- function(long) {
  ##This function makes a sigmoid colormap - better for heatmaps than the std for R linear colormaps
  x=seq(-2*pi, 2*pi, length.out=long)
  y=tanh(x)
  r=(y<0)*-y
  g=(y>0)*y

  newy=rgb(r,g,0)

  newy
}



rescale_cpg<- function(data,numdev=3) {
##Rescale the data along rows akin to MATLAB
  rescaled=data
  for (i in 1:nrow(data)) {
    mid=mean(data[i,])
    dev=sd(data[i,])
    rescaled[i,]=rescaled[i,]-mid;
    rescaled[i,(rescaled[i,]>(numdev*dev))]=numdev*dev
    rescaled[i,(rescaled[i,]<(-numdev*dev))]=-numdev*dev
  }

  return(rescaled)
}

simplecluster <- function(avg, sampy, which_tissue=3, which_type=NA, which_grade=NA, which_progression=NA, divider=NA, namey="a.pdf")  {
 
  ##Initialize with if not set to all possibles
  if (is.na(which_type[1])) {
    which_type=unique(sampy$Type)
  }
  if (is.na(which_grade[1])) {
    which_grade=unique(sampy$Grade)
  }
  if (is.na(which_progression[1])) {
    which_progression=unique(sampy$Progression)
  }
  if (is.na(divider[1])) {
    divider=sampy$Progression>0
  }

        


  cgoody=is.element(sampy$Class,which_tissue)&is.element(sampy$Type,which_type) & is.element(sampy$Grade, which_grade) & is.element(sampy$Progression, which_progression)
  ngoody=is.element(sampy$Class,which_tissue)&is.element(sampy$Type,which_type) & is.element(sampy$Grade, which_grade) & is.element(sampy$Progression, 0)
  noody=avg[,ngoody]
  rowo=apply(noody, 1, median)
  colony=avg[order(rowo),cgoody]



  newy=sigmoid_color(63)
  
  divider=factor(divider[cgoody])
  levels(divider)=rainbow(nlevels(divider))
  
  pdf(paste("Movie/", namey, sep=""), width=11, height=8)
  #heatmap.2(as.matrix(colony),col=colorpanel(255,low="black",high="green"),trace="none",hclustfun=hclust.ward, ColSideColors=as.character(divider),key=T)
  #col.heat <- heatmap.2(as.matrix(colony), dendrogram="column", trace="none", scale="none", labCol=sampy$Sample_ID[cgoody], hclustfun=hclust.ward,
  ##ColSideColors=as.character(divider), col=brewer.pal(11,"RdYlBu"))

  r=seriate(dist(colony), method="OLO", control=list(method="ward"))
  c=seriate(dist(t(colony)), method="OLO", control=list(method="ward"))
  ##col.heat <- heatmap.2(as.matrix(colony), dendrogram="both", trace="none", labCol=sampy$Sample_ID[cgoody],ColSideColors=as.character(divider),  labRow="", col=brewer.pal(11, "RdYlBu"), hclustfun=hclust.ward)
  
  col.heat <- heatmap.2(as.matrix(colony), dendrogram="column", trace="none", labCol=sampy$Sample_ID[cgoody],ColSideColors=as.character(divider),  labRow="", col=brewer.pal(11, "RdYlBu"), hclustfun=hclust.ward,
                        Rowv=as.dendrogram(r[[1]]), Colv=as.dendrogram(c[[1]]) )


  
  
  dev.off()


}


library(gplots)
  
library(RColorBrewer)

library(seriation)
##Main function - generate colon tn plot
