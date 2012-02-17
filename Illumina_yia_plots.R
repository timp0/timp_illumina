##This code is used to make thyroid plots I want

##dev.off()
library(plotrix)

## Maha Dist from Normal

#dat_med=apply(data$fqbeta[,norms],1,median)
#dat_mad=apply(data$fqbeta[,norms],1,mad)
#dat_mad=apply(data$fqbeta[,carcy],1,mad)

#comp_mad=dat_mad-dat_mad

#dat_mad=apply(data$fqbeta, 1, mad)

#top25=order(comp_mad, decreasing="T")[1:25]

hector=read.csv('/Users/timp/Big_Data/Illumina_Bead/Analysis/ordered_cpgs.csv', stringsAsFactors=F)

top25=hector$X[1:25]

##PCA

p=prcomp(t(data$fqbeta[top25,]))

##MDS

z=cmdscale( dist( t( data$fqbeta[top25,]) ), k=2)

norms=data$fsamp$Progression==0
carcy=data$fsamp$Progression>2
ady=data$fsamp$Progression>0&data$fsamp$Progressoin<3

##Now thyroid progression plot
colly=c("green", "orange", "blue", "red","brown")

type_col=factor(data$fsamp$Class)

types=levels(type_col)

levels(type_col)=colly

pnorms=p$x[norms,]
pn_med=apply(pnorms,2,median)
pn_mad=apply(pnorms,2,mad)

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

plot(z[carcy,1], z[carcy,2], bg=as.character(type_c?mol[norms]), pch=21, xlim=mds_range_x, ylim=mds_range_y)

draw.ellipse(x=nprof$xcoor, y=nprof$ycoor, a=nprof$xrad, b=nprof$yrad, lty=2, lwd=2,border=nprof$col)

legend("topright", c("Breast", "Colon", "Lung", "Wilms", "Thyroid"), col=as.character(levels(type_col)), pch=16)

dev.off()

#Colon adenoma/tumor sample stuff
colons=data$fsamp$Class==3

colony=list(fqbeta=data$fqbeta[,colons], fsamp=data$fsamp[colons,])

##MDS

z=cmdscale( dist( t( colony$fqbeta[top25,]) ), k=2)

norms=colony$fsamp$Progression==0
carcy=colony$fsamp$Progression>2
ady=colony$fsamp$Progression>0&colony$fsamp$Progression<3

##Now progression plot
pcolly=c("blue", "green", "red")

prog_col=factor(ifelse(colony$fsamp$Progression>0,1,0)+ifelse(colony$fsamp$Progression>2,1,0))

levels(prog_col)=pcolly

pnorms=p$x[norms,]
pn_med=apply(pnorms,2,median)
pn_mad=apply(pnorms,2,mad)

ztn=z[norms,]
c_x=median(ztn[,1])
c_y=median(ztn[,2])
c_a=mad(ztn[,1])*3
c_b=mad(ztn[,2])*3


mds_range_x=range(z[,1])
mds_range_y=range(z[,2])

pdf("Movie/colons_prog_mds.pdf")

plot(z[,1], z[,2], bg=as.character(prog_col), pch=21, xlim=mds_range_x, ylim=mds_range_y)

draw.ellipse(x=rep(c_x,2), y=rep(c_y,2), a=rep(c_a,2), b=rep(c_b,2), lty=2, lwd=2, border="green")

legend("topleft", c("Normal", "Adenoma", "Carcinoma"), col=as.character(levels(prog_col)), pch=16)

dev.off()






#anti_norm=sweep( abs(sweep(thy_data$fqbeta,1,thyn_med)), 1, thyn_mad, FUN="/")

#anti_norm_score=apply(anti_norm[top25,],2,sum)


#prog_x=factor(thy_data$fsamp$Progression)

#levels(prog_x)=c(1,2,3,4)

#pdf("Movie/thy_prog_maha_scat.pdf")

#plot(jitter(as.numeric(prog_x)), anti_norm_score, bg=as.character(prog_col), pch=21, xaxt="n", xlab="")
#legend("topleft", c("Normal", "Nodule", "Adenoma", "Carcinoma"), col=levels(prog_col), pch=16)

#dev.off()

#source("~/Work/Analysis/Illumina/Illumina_cluster1b.R")


##simplecluster <- function(avg, sampy, which_tissue=3, which_type=NA, which_grade=NA, which_progression=NA, divider=NA, namey="a.pdf")  {
#simplecluster(thy_data$fqbeta, thy_data$fsamp, which_tissue=7, divider=thy_data$fsamp$Progression, namey="thy_cluster_prog.pdf")
#simplecluster(thy_data$fqbeta, thy_data$fsamp, which_tissue=7, divider=abs(thy_data$fsamp$Other_Note), namey="thy_cluster_type.pdf")

