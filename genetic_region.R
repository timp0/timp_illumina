##Run using R CMD BATCH findseq1.R
##Finding chromosomal locations of sequences

plotucsc <-function(ucsc_isl,chromy,start,end){
  ##Check for UCSC Island
  in_ucscisl=ucsc_isl[( (ucsc_isl$X.chrom==chromy)& ( ((ucsc_isl$chromStart>start)&(ucsc_isl$chromStart<end))|((ucsc_isl$chromEnd>start)&(ucsc_isl$chromEnd<end)))),]
  
  ##Plot UCSC Island
  if (nrow(in_ucscisl)>0) {
    for(j in 1:nrow(in_ucscisl)){
      ##cat("Island in ", i, "\n")
      isl_start=in_ucscisl$chromStart[j]
      isl_end=in_ucscisl$chromEnd[j]
      
      polygon(c(isl_start,isl_end,isl_end,isl_start),c(0.0125,0.0125,0.05,0.05),density=10,col="blue",angle=45)
    }
  }
}

plothmm <-function(hmm_isl, chromy, start, end){
  ##Check for HMM Island
  in_hmmisl=hmm_isl[( (hmm_isl$chr==chromy)& ( ((hmm_isl$start>start)&(hmm_isl$start<end))|((hmm_isl$end>start)&(hmm_isl$end<end)))),]
  
  ##plot HMM Island
  if (nrow(in_hmmisl)>0) {
    for(j in 1:nrow(in_hmmisl)){
      ##cat("Island in ", i, "\n")
      isl_start=in_hmmisl$start[j]
      isl_end=in_hmmisl$end[j]
      
      polygon(c(isl_start,isl_end,isl_end,isl_start),c(0,0,0.0375,0.0375),density=10,col="red",angle=-45)
    }
  }
}

plotcpgdens <- function(chromy,start,end){
  ##Get the chromosome sequence
  seq <- Hsapiens[[chromy]]
  ##cat("Plot CPGDens ",chromy, start, end)
  ##Subset the sequence
  subseq <- subseq(seq,start=start,end=end)
  ##Find CpGs
  cpgs=start(matchPattern("CG",subseq))+start-1
  ##Create the CpG bins
  cuts=seq(start,end,8)
  ##Bin the CpGs
  scpgs=cut(cpgs,cuts,include.lowest=TRUE)
  
  x = (cuts[1:(length(cuts)-1)]+cuts[2:(length(cuts))])/2
  y = table(scpgs)/8
  SPAN=400/diff(range(x))
  d = loess(y~x,span=SPAN,degree=1)
  
  plot(x,d$fitted,type="l",ylim=c(0,0.2),xlab="",
       ylab="CpG density",xlim=c(start,end),xaxt="n")
  
  rug(cpgs)
}

plotcontrols <- function(avg_data,sampy,not_far,region,start,end) {

  ##Plot actual probe locations
  ##Control_samples

  cntrl_data<-avg_data[(probes$Region==region),(sampy$Class==1)]
  cntrl_info<-sampy[sampy$Class==1,]
  
  plot(0,0,ylim=c(0,1),xlim=c(start,end),"axes"=FALSE, type="n",xlab="",ylab="")
  axis(side=4)

  for (j in 1:nrow(not_far)) {
    ##Plot control points
    loc=not_far$Start_loc[j]-50+cntrl_info$Other_Note-start
    points(loc+start,cntrl_data[j,],col="green")
    ##Linear Fit
    d=lm((as.numeric(cntrl_data[j,]))~loc)
    ##plot on the same scale
    x=unique(loc)
    y=d$coefficients[[1]]+d$coefficients[[2]]*x
    lines(x+start,y,col="green")
    ##Plot smoothed curve
    lines(lowess(loc+start,as.numeric(cntrl_data[j,])),col="red",lty=2)
  } 
  ##Plot label on axis as a tick on the top
  axis(side=1,at=not_far$Start_loc,labels=not_far$Probe_ID, cex.axis=0.8)
}

plotsamples <- function(class_stuff,probes,avg_data,sampy,not_far,region,start) {
  for (j in 1:length(class_stuff$nums)) {
    data<-avg_data[(probes$Region==region),(sampy$Class==class_stuff$nums[j])]
    samp_info<-sampy[sampy$Class==class_stuff$nums[j],]
    ##cat(class_stuff$nums[j],"\n")
    num_samp=ncol(data)
    for (k in 1:nrow(not_far)) {
      loc=jitter(rep(not_far$Start_loc[k]-start,num_samp),amount=50)
      points(loc+start,data[k,],col=class_stuff$coloring[j])
      ##cat("X: ",loc+start," Y: ",data[k,1],"\n")
    }
  }
  
  legend("topright",class_stuff$nam,col=class_stuff$coloring,pch=1)
  
}

    
plotgenes <- function(genes,chromy, start, end){
  ##Find Genes in Region
  in_genes=genes[( (genes$chrom==chromy)& ( ((genes$txStart>start)&(genes$txStart<end))|((genes$txEnd>start)&(genes$txEnd<end)))),]    
  ##Plot Genes
  plot(0,0,ylim=c(-1.5,1.5),xlim=c(start,end),yaxt="n",ylab="Genes",
       xlab="")
  ##Label Y-axis
  axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1)
  ##Add dotted line
  abline(h=0,lty=3)

  #Plot each gene
  if(nrow(in_genes)>0){
    for(j in 1:nrow(in_genes)){
      TS = in_genes$txStart[j]
      TE = in_genes$txEnd[j]
      CS = in_genes$cdsStart[j]
      CE = in_genes$cdsEnd[j]
      ES = as.numeric(strsplit(in_genes$exonStarts[j],",")[[1]])
      EE = as.numeric(strsplit(in_genes$exonEnds[j],",")[[1]])
      Exons=cbind(ES,EE)
      Strand=ifelse(in_genes$strand[j]=="+",1,-1)
      polygon(c(TS,TE,TE,TS),Strand/2+c(-0.4,-0.4,0.4,0.4),density=0,col=j)
      polygon(c(CS,CE,CE,CS),Strand/2+c(-0.25,-0.25,0.25,0.25),density=0,
              col=j) 
      apply(Exons,1,function(x) polygon(c(x[1],x[2],x[2],x[1]),
                                        Strand/2+c(-0.2,-0.2,0.2,0.2),col=j))     
      text((max(start,TS)+min(end,TE))/2,Strand*0.9,in_genes$X.geneName[j],cex=1,
           pos=Strand+2)
      ##Strand +2 works cuase -1 turns into 1 (below) and 1 
      ##turns into 3 (above)
    }
  }
}

plotClass <- function(probes,ucsc_isl,hmm_isl,genes,avg_data,sampy,class_stuff,which,outname,controls=TRUE) {
  ##
  ##Start Plot
  pdf(outname,width=11,height=8)
  par(mfrow=c(3,1))
  layout(c(1:3),heights=c(0.5,0.3,0.2))
  par(mar=c(2,1,0,1),mgp=c(1.5,.5,0),oma=c(0,1,2,1)) 
  
  
  ##Loop through all regions
  for (i in unique(probes$Region)) {
    ##Find probes which are in the same region - this was previously defined by perl in the probes$Region command
    not_far<-probes[(probes$Region==i),]
    num_close<-nrow(not_far)
    
    ##Set which chromosome and coordinates we are using for the region
    chromy <- paste("chr",not_far$Chromosome[1],sep="")
    index=(min(not_far$Start_loc)-1000):(max(not_far$Finish_loc)+1000)
    final_index=length(index)
    plot(0,0,ylim=c(0,1.1),xlim=c(index[1],index[final_index]),ylab="Methylation",xlab="",type="n",axes=F)
    box()
    axis(2, at=seq(0,1,.2))
    ##Plot actual samples
    plotsamples(class_stuff[which,],probes,avg_data,sampy,not_far,i,index[1])
    ##Plot label on axis as a tick on the bottom
    axis(side=1,at=not_far$Start_loc,labels=not_far$Probe_ID, cex.axis=0.8)
    
    ##Plot title to graph
    mtext(paste("ID:",i,"--",as.character(chromy),":",index[1],"-",index[final_index],sep=""),cex=2,side=3,outer=TRUE)

    ##Plot CpGDensity
    plotcpgdens(chromy,index[1],index[final_index])

    ##Plot Controls
    if (controls) {
      ##Stay on same plot
      par(new=T)
      plotcontrols(avg_data,sampy,not_far,i,index[1],index[final_index])
    }
    ##Plot Island locations
    plotucsc(ucsc_isl,chromy, index[1],index[final_index]) 
    plothmm(hmm_isl,chromy,index[1],index[final_index])
    legend("topleft",c("UCSC Islands", "HMM Islands"),col=c("blue", "red"),lty=1,lwd=2)
    
    ##Plot RefSeq Gene regions
    plotgenes(genes,chromy, index[1], index[final_index])
  }
  dev.off()
}



plotfullTN <- function(probes,ucsc_isl,hmm_isl,genes,avg_data,sampy,class_stuff,outname,controls=TRUE) {
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
    
    ##Set which chromosome and coordinates we are using for the region
    chromy <- paste("chr",not_far$Chromosome[1],sep="")
    ##cat("Region ", i, "\n")
    index=(min(not_far$Start_loc)-1000):(max(not_far$Finish_loc)+1000)
    final_index=length(index)

    
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

    ##Plot CpGDensity
    plotcpgdens(chromy,index[1],index[final_index])

    ##Plot Island locations
    plotucsc(ucsc_isl,chromy, index[1],index[final_index]) 
    plothmm(hmm_isl,chromy,index[1],index[final_index])
    legend("topleft",c("UCSC Islands", "HMM Islands"),col=c("blue", "red"),lty=1,lwd=2)

    ##Plot Controls
    if (controls) {
      ##Stay on same plot
      par(new=T)
      plotcontrols(avg_data,sampy,not_far,i,index[1],index[final_index])
    }
    
    ##Plot RefSeq Gene regions
    plotgenes(genes,chromy, index[1], index[final_index])
  }
  dev.off()
}


##Load libraries in
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)

##Load CSV file
probes<-read.csv("New_norm_mat_targer_ucscisl_hmmisl.csv",stringsAsFactors=FALSE)
old_probes<-read.csv("Old_Illumina_pri_mat_targer_ucscisl_hmmisl.csv",stringsAsFactors=FALSE)
##Load UCSC Isl file
ucsc_isl<-read.delim("ucsc_cpgisl.txt",stringsAsFactors=FALSE)
##Load HMM Isl file
hmm_isl<-read.csv("hmm_isl1.txt",stringsAsFactors=FALSE)  
##Load RefSeq Gene file
genes<-read.delim("ref_genes.txt",stringsAsFactors=FALSE)
##Load Data file
raw<-read.csv("New_norm_mat_data.csv",stringsAsFactors=FALSE)
raw_old<-read.csv("Old_Illumina_pri_mat_data.csv",stringsAsFactors=FALSE)
##Get Avg Data
avg_data<-raw[,seq(2,ncol(raw)-1,by=5)]
old_avg_data<-raw_old[,seq(2,ncol(raw_old)-1,by=5)]
##Load Sample file
sampy<-read.csv("New_norm_mat_sample.csv",stringsAsFactors=FALSE)
old_sampy<-read.csv("Old_Illumina_pri_mat_sample.csv",stringsAsFactors=FALSE)



##Define Classes/Colors
class_stuff=data.frame(nums=c(1,2,3,4,5,6,7,8,9,10,11),
  nam=c("Controls","Breast Tumor","Breast Normal","Colon Tumor", "Colon Normal","Lung Tumor","Lung Normal","Ovary Tumor","Ovary Normal","Wilms Tumor","Wilms Normal"),
  coloring=c("black","coral1","coral4","cadetblue1","cadetblue4","goldenrod1","goldenrod4","seagreen1","seagreen4","plum1","plum4"),stringsAsFactors=FALSE)


##Plot new data all sets
plotfullTN(probes,ucsc_isl,hmm_isl,genes,avg_data,sampy,class_stuff,"Movie/full_bad.pdf")
##Plot old data wilms/colon(all I have right now)
plotfullTN(old_probes,ucsc_isl,hmm_isl,genes,old_avg_data,old_sampy,class_stuff,"Movie/old_full.pdf",FALSE)
##Plot new normals against each other
##plotClass(probes,ucsc_isl,hmm_isl,genes,avg_data,sampy,class_stuff,c(3,5,7,9,11),"Movie/normals_R.pdf")
#plotTN(probes,ucsc_isl,hmm_isl,genes,avg_data,sampy,class_stuff,c(6:7),"lung_R.pdf")
#plotTN(probes,ucsc_isl,hmm_isl,genes,avg_data,sampy,class_stuff,c(8:9),"ovary_R.pdf")
#plotTN(probes,ucsc_isl,hmm_isl,genes,avg_data,sampy,class_stuff,c(10:11),"wilms_R.pdf")
