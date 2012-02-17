setwd("~/Data/Illumina_Bead/Analysis")
source("data_run1a.R")

source("~/Work/Analysis/Illumina/Illumina_CHARM3a.R")
source("~/Work/Analysis/Illumina/Illumina_annotation1.R")

##Load libraries in
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)
library(plotrix)

##Hector's list of ranked CpGs - just as something to use
horder<-read.csv("ordered_cpgs.csv",stringsAsFactors=FALSE)  

top25=horder$X[1:25]

##SupData4 from 2009 Nat Gen paper
##natgen2009<-read.csv("cDMRs_NatGen09.csv", stringsAsFactors=FALSE)
##SupData2 from 2009 Nat Gen paper
##natgen2009<-read.csv("fullcDMRs_NatGen09.csv", stringsAsFactors=FALSE)

##natgen2009$chromy=substr(natgen2009$chr, 4, nchar(natgen2009$chr))

##Ok - find the right CHARM probes for the DMRs
##natprobes=CHARMmatch(natgen2009$chromy, (natgen2009$start+natgen2009$end)/2,charm_nfo)



load("natgenprobes.rda")


##Load UCSC Isl file
ucsc_isl<-read.delim("ucsc_cpgisl.txt",stringsAsFactors=FALSE)


  
heatCHARMcancer(CHIL)
##heatCHARMtissue(natprobes)


pdf("Movie/CHARM_tissue_allregion1.pdf", version="1.4")

##layout(rbind(1,2),heights=c(0.8,0.2))
##par(mar=c(2,1,0,1),mgp=c(1.5,.5,0),oma=c(0,1,2,1)) 



for (i in 1:290)
{
  ##Define Region I want to plot the CHARM results in:
  ##Let's do DUSP5 - Illumina probe DUSP5-2 is one of the best
  ##It's Chromosome 10, 112249240
  chromy=paste("chr", horder$Chromosome[i], sep="")
  ##Take ~+/- 1000 bp
  index=c(horder$Start_loc[i]-1000, horder$Start_loc[i]+1000)    
  ##Make initial plot
  plot(0,0,ylim=c(-0.25,1.5),xlim=c(index[1],index[2]),ylab="Methylation",xlab="",type="n",axes=F)
  box()
  axis(2)
  axis(1)
  title(paste("Chromosome ", horder$Chromosome[i], sep=""))
  ##Plot CHARM Data
  plotCHARM(chromy,index[1],index[2], c(1:5), pts=F)
  ##Plot label on axis as a tick on the bottom
  ##axis(side=1,at=not_far$Start_loc,labels=not_far$Probe_ID, cex.axis=0.8)
  ##Plot UCSC Island as orangle translucent rectangle
  plotucsc(ucsc_isl,chromy,index[1],index[2])

  ##Plot selected Cpgs/probes as a dotted vertical line
  abline(v=horder$Start_loc[i], lty=2)
  ##Plot cpgs as a rug
  cpgsrug(chromy,index[1],index[2])


}
dev.off()


##Get Table of samples
stable=data.frame(Normal=numeric(5), Adenoma=numeric(5), Adenocarcinoma=numeric(5))
sclasses=c(2,3,4,6,7)
rownames(stable)=c("Breast", "Colon", "Lung", "Wilms", "Thyroid")


for (i in 1:5) {
  stable$Normal[i]=sum((data$fsamp$Class==sclasses[i])&(data$fsamp$Progression==0))
  stable$Adenoma[i]=sum((data$fsamp$Class==sclasses[i])&((data$fsamp$Progression>0)&(data$fsamp$Progression<3)))
  stable$Adenocarcinoma[i]=sum((data$fsamp$Class==sclasses[i])&(data$fsamp$Progression>2))

}


##Thyroid progression plot - needed - coming from illumina_yia_plots1
##Plot normal v adenoma v cancer

thyroids=data$fsamp$Class==7

thyry=list(fqbeta=data$fqbeta[,thyroids], fsamp=data$fsamp[thyroids,])

##MDS

z=cmdscale( dist( t( thyry$fqbeta[top25,]) ), k=2)

norms=thyry$fsamp$Progression==0
carcy=thyry$fsamp$Progression>2
ady=thyry$fsamp$Progression>0&thyry$fsamp$Progression<3

##Now progression plot
pthy=c("blue", "green", "red")

prog_thy=factor(ifelse(thyry$fsamp$Progression>0,1,0)+ifelse(thyry$fsamp$Progression>2,1,0))

levels(prog_thy)=pthy

#pnorms=p$x[norms,]
#pn_med=apply(pnorms,2,median)
#pn_mad=apply(pnorms,2,mad)

ztn=z[norms,]
c_x=median(ztn[,1])
c_y=median(ztn[,2])
c_a=mad(ztn[,1])*3
c_b=mad(ztn[,2])*3


mds_range_x=range(z[,1])
mds_range_y=range(z[,2])

pdf("Movie/thyroid_prog_mds.pdf")

plot(z[,1], z[,2], bg=as.character(prog_thy), pch=21, xlim=mds_range_x, ylim=mds_range_y)

draw.ellipse(x=rep(c_x,2), y=rep(c_y,2), a=rep(c_a,2), b=rep(c_b,2), lty=2, lwd=2, border="green")

legend("topleft", c("Normal", "Adenoma", "Carcinoma"), col=as.character(levels(prog_thy)), pch=16)

dev.off()




