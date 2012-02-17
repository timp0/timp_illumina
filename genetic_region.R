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
      
      polygon(c(isl_start,isl_end,isl_end,isl_start),c(0.05,0.05,0.1,0.1),density=10,col="blue",angle=45)
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
      
      polygon(c(isl_start,isl_end,isl_end,isl_start),c(0,0,0.05,0.05),density=10,col="red",angle=-45)
    }
  }
}

plotcpgdens <- function(chromy,start,end){
  ##Get the chromosome sequence
  seq <- Hsapiens[[chromy]]
  
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
  
  plot(x,d$fitted,type="l",ylim=c(0,0.2),xlab="Location",
       ylab="CpG density",xlim=c(start,end))
  
  rug(cpgs)
}

plotgenes <- function(genes,chromy, start, end){
  ##Find Genes in Region
  in_genes=genes[( (genes$chrom==chromy)& ( ((genes$txStart>start)&(genes$txStart<end))|((genes$txEnd>start)&(genes$txEnd<end)))),]    
  ##Plot Genes
  plot(0,0,ylim=c(-1.5,1.5),xlim=c(start,end),yaxt="n",ylab="Genes",
       xlab="Location")
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

  
##Load libraries in
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)

##Load CSV file
probes<-read.csv("New_norm_mat_targer_ucscisl_hmmisl.csv",stringsAsFactors=FALSE)
##Load UCSC Isl file
ucsc_isl<-read.delim("ucsc_cpgisl.txt",stringsAsFactors=FALSE)
##Load HMM Isl file
hmm_isl<-read.csv("hmm_isl1.txt",stringsAsFactors=FALSE)  
##Load RefSeq Gene file
genes<-read.delim("ref_genes.txt",stringsAsFactors=FALSE)
##Sort the probe info based on chromosome
sorted<-probes[order(probes$Chromosome,probes$Start_loc) ,]



##Start Plot




pdf("b.pdf",width=11,height=8)
par(mfrow=c(2,1))
par(mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0)) 


i <- 1
while (i <= nrow(sorted)) {
  ##Find probes which have the same start name(they will be the same)
  ##cat(i);
  ##cat(" running.\n")
  not_far<-sorted[((sorted$Chromosome[i]==sorted$Chromosome)&(abs(sorted$Start_loc[i]-sorted$Start_loc)<10000)),]
  num_close<-nrow(not_far)
  
  chromy <- paste("chr",not_far$Chromosome[1],sep="")
  index=(not_far$Start_loc[1]-1000):(not_far$Finish_loc[num_close]+1000)
  final_index=length(index)
  
  ##Plot CpGDensity
  plotcpgdens(chromy,index[1],index[final_index])

  ##Plot actual probe locations
  points(not_far$Start_loc, rep(0.07, num_close))
  text(not_far$Start_loc,rep(0.07,num_close),labels=not_far$Probe_ID, cex=0.5, pos=1)
  

  ##Plot Island locations
  plotucsc(ucsc_isl,chromy, index[1],index[final_index]) 
  plothmm(hmm_isl,chromy,index[1],index[final_index])
  legend("topright",c("UCSC Islands", "HMM Islands"),col=c("blue", "red"),lty=1,lwd=2)
  
  mtext(paste("ID:",i,"--",as.character(chromy),":",index[1],"-",index[final_index],sep=""),cex=2)

  plotgenes(genes,chromy, index[1], index[final_index])
  
  
  
  
  i=i+num_close;
}
dev.off()
