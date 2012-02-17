##Ok - make some plots specific for webpage
##Adding mds and dot plot

setwd("~/Data/Illumina_Bead/Analysis/")

load("predata.rda")

library(RColorBrewer)

library(plotrix)

hector=read.csv("~/Data/Illumina_Bead/Analysis/ordered_cpgs.csv", stringsAsFactors=F)

top25=hector$X[1:25]

##PCA

p=prcomp(t(data$fqbeta[top25,]))

##MDS

z=cmdscale( dist( t( data$fqbeta[top25,]) ), k=2)

norms=data$fsamp$Progression==0
carcy=data$fsamp$Progression>2
ady=data$fsamp$Progression>0&data$fsamp$Progressoin<3

##MDS colors
colly=c("green", "orange", "blue", "red","brown")

type_col=factor(data$fsamp$Class)

types=levels(type_col)

levels(type_col)=colly


##Get coordinates and ranges of ellipse defining normal range


nprof=list(xcoor=rep(0,5), ycoor=rep(0,5), xrad=rep(0,5), yrad=rep(0,5), col=rep(as.character("black"),5))

for (i in 1:5) {
  type_norms=data$fsamp$Class==types[i]&data$fsamp$Progression==0
  ztn=z[type_norms,]
  nprof$xcoor[i]=median(ztn[,1])
  nprof$ycoor[i]=median(ztn[,2])
  nprof$xrad[i]=mad(ztn[,1])*3
  nprof$yrad[i]=mad(ztn[,2])*3
  nprof$col[i]=colly[i]
}

pca_range_x=range(p$x[,1])
pca_range_y=range(p$x[,2])

mds_range_x=range(z[,1])
mds_range_y=range(z[,2])

pdf("Movie/norms_mds.pdf")

plot(z[norms,1], z[norms,2], bg=as.character(type_col[norms]), pch=21, xlim=mds_range_x, ylim=mds_range_y)

draw.ellipse(x=nprof$xcoor, y=nprof$ycoor, a=nprof$xrad, b=nprof$yrad, lty=2, lwd=2,border=nprof$col)

legend("topright", c("Breast", "Colon", "Lung", "Wilms", "Thyroid"), col=as.character(levels(type_col)), pch=16)

dev.off()

pdf("Movie/carc_mds.pdf")

plot(z[carcy,1], z[carcy,2], bg=as.character(type_col[norms]), pch=21, xlim=mds_range_x, ylim=mds_range_y)

draw.ellipse(x=nprof$xcoor, y=nprof$ycoor, a=nprof$xrad, b=nprof$yrad, lty=2, lwd=2,border=nprof$col)

legend("topright", c("Breast", "Colon", "Lung", "Wilms", "Thyroid"), col=as.character(levels(type_col)), pch=16)

dev.off()

## ok - for dot plot, need to set y values to the different tissue types, and add jitter 
tissue.y=as.numeric(type_col)*2-as.numeric(carcy)
tissue.y=jitter(tissue.y)

pdf("Movie/dotplot_mds.pdf")

plot(z[(norms|carcy),1], tissue.y[(norms|carcy)], bg=as.character(type_col[(norms|carcy)]), pch=21, xlim=mds_range_x, ylim=c(0,11))

##draw.ellipse(x=nprof$xcoor, y=nprof$ycoor, a=nprof$xrad, b=nprof$yrad, lty=2, lwd=2,border=nprof$col)

legend("topright", c("Breast", "Colon", "Lung", "Wilms", "Thyroid"), col=as.character(levels(type_col)), pch=16)

dev.off()


##Just Breast for Rafa

pdf("Movie/just_breast_norm.pdf")

plot(z[(norms&(data$fsamp$Class==2)),1],z[(norms&(data$fsamp$Class==2)),2], pch=21, xlim=c(-.75, .25), ylim=c(-.6, .1))

draw.ellipse(x=nprof$xcoor[1], y=nprof$ycoor[1], a=nprof$xrad[1], b=nprof$yrad[1], lty=2, lwd=2,border=nprof$col[1])

dev.off()


##Find norms and cancers
norms=data$fsamp$Progression==0
carcy=data$fsamp$Progression>2

##Far is 2, shore is 1, island is 0
islstatus=(data$probes$UCSC_Dist_to_Island>0)+(data$probes$UCSC_Dist_to_Island>2000)


cgi.colors <- brewer.pal(8,"Dark2")[c(4:5,8)]

normvar=apply(data$fqbeta[,norms],1,mad)
cancvar=apply(data$fqbeta[,carcy],1,mad)

pdf("Movie/allvar1.pdf")

rangy=max(max(normvar), max(cancvar))

plot(normvar,cancvar, xlab="Normal", ylab="Cancer",
     xlim=c(0, rangy), ylim=c(0,rangy),
       bg=cgi.colors[islstatus+1], pch=21)


##Signifcance lines
cc <- qf(.99, sum(carcy)-1, sum(norms)-1)
abline(0,sqrt(cc), lty=2)
abline(0,1)

legend("bottomright", c("Island", "Shore", "Far"), pch=21, pt.bg=cgi.colors)

dev.off()

pdf("Movie/tisvar1.pdf")


for (i in c(2, 3, 4, 6, 7)) {
  norms=(data$fsamp$Progression==0)&(data$fsamp$Class==i)
  carcy=(data$fsamp$Progression>2)&(data$fsamp$Class==i)
  

  cgi.colors <- brewer.pal(8,"Dark2")[c(4:5,8)]
  
  normvar=apply(data$fqbeta[,norms],1,mad)
  cancvar=apply(data$fqbeta[,carcy],1,mad)
  

  
  rangy=max(max(normvar), max(cancvar))
  
  plot(normvar,cancvar, xlab="Normal", ylab="Cancer",
       xlim=c(0, rangy), ylim=c(0,rangy),
       bg=cgi.colors[islstatus+1], pch=21,
       main=i)
  
  
  ##Signifcance lines
  cc <- qf(.99, sum(carcy)-1, sum(norms)-1)
  abline(0,sqrt(cc), lty=2)
  abline(0,1)
  
  legend("bottomright", c("Island", "Shore", "Far"),
         pch=21, pt.bg=cgi.colors)  
  
}

dev.off()
