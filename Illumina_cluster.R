##This code is used to setup clustering for Illumina samples

hclust.ward <- function(x) {
  hclust(x, method="ward")
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



rescale_cpg<- function(data) {
##Rescale the data along rows akin to MATLAB
  rescaled=data
  for (i in 1:nrow(data)) {
    mid=mean(data[i,])
    dev=sd(data[i,])
    rescaled[i,]=rescaled[i,]-mid;
    rescaled[i,(rescaled[i,]>(3*dev))]=1*dev
    rescaled[i,(rescaled[i,]<(-3*dev))]=-1*dev
  }

  rescaled
}

library(gplots)
##Main function - generate colon tn plot

cgoody=( (fusesampy$Class==3)&( (fusesampy$Other_Note==-1) | (fusesampy$Other_Note==1) | (fusesampy$Other_Note==5) ) )
colony=fuseavg[,cgoody]
newy=sigmoid_color(255)
tn=factor(fusesampy$Other_Note[cgoody])
levels(tn)=rainbow(nlevels(tn))

pdf("Movie/colonRclust_unscale.pdf", width=11, height=8)
heatmap.2(as.matrix(colony),col=colorpanel(255,low="black",high="green"),trace="none",hclustfun=hclust.ward, ColSideColors=as.character(tn),key=T)
dev.off()
