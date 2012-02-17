#Run using R CMD BATCH findseq1.R
#Finding chromosomal locations of sequences


#Load libraries in
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)

#Load CSV file
probes<-read.csv("New_norm_mat_targer_ucscisl_hmmisl.csv",stringsAsFactors=FALSE)
genes<-read.delim("ref_genes.txt",stringsAsFactors=FALSE)
sorted<-probes[order(probes$Chromosome,probes$Start_loc) ,]

#Start Plot




pdf("b.pdf",width=11,height=8)
i <- 1
while (i <= nrow(sorted)) {
  #Find probes which have the same start name(they will be the same)
  par(mfrow=c(2,1))
  cat(i);
  cat(" running.\n")
  not_far<-sorted[((sorted$Chromosome[i]==sorted$Chromosome)&(abs(sorted$Start_loc[i]-sorted$Start_loc)<10000)),]
  num_close<-nrow(not_far)
  
  chromy <- paste("chr",not_far$Chromosome[1],sep="")
  #Get the chromosome sequence
  seq <- Hsapiens[[chromy]]
  
  index=(not_far$Start_loc[1]-1000):(not_far$Finish_loc[num_close]+1000)
  final_index=length(index)

  #Subset the sequence
  subseq <- subseq(seq,start=index[1],end=index[final_index])
  #Find CpGs
  cpgs=start(matchPattern("CG",subseq))+index[1]-1
  #Create the CpG bins
  cuts=seq(index[1],index[final_index],8)
  #Bin the CpGs
  scpgs=cut(cpgs,cuts,include.lowest=TRUE)
  
  x = (cuts[1:(length(cuts)-1)]+cuts[2:(length(cuts))])/2
  y = table(scpgs)/8
  SPAN=400/diff(range(x))
  d = loess(y~x,span=SPAN,degree=1)

  plot(x,d$fitted,type="l",ylim=c(0,0.15),xlab="Location",
           ylab="CpG density",xlim=c(index[1],index[final_index]))
  for(j in c(1:num_close)) {
    points(not_far$Start_loc[j], .07)
  }
  rug(cpgs)


  #Find Genes in Region
  in_genes=genes[( (genes$chrom==chromy)& ( ((genes$txStart>index[1])&(genes$txStart<index[final_index]))|((genes$txEnd>index[1])&(genes$txEnd<index[final_index])))),]
  
  
  #Plot Genes
  plot(0,0,ylim=c(-1.5,1.5),xlim=c(index[1],index[final_index]),yaxt="n",ylab="Genes",
             xlab="Location")
  #Label Y-axis
  axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1)
  #Add dotted line
  abline(h=0,lty=3)
  
  
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
      
      text((max(index[1],TS)+min(index[final_index],TE))/2,Strand*0.9,in_genes$X.geneName[j],cex=1,
           pos=Strand+2)
      ##Strand +2 works cuase -1 turns into 1 (below) and 1 
      ##turns into 3 (above)
    }
  }
  mtext(paste("ID:",i,"--",as.character(chromy),":",index[1],"-",index[final_index],sep=""),
        side=3,cex=1.5,outer=TRUE)



  i=i+num_close;
}
dev.off()
