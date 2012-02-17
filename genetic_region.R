#Run using R CMD BATCH findseq1.R
#Finding chromosomal locations of sequences


#Load libraries in
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)

#Load CSV file
probes<-read.csv("New_norm_mat_targer_ucscisl_hmmisl.csv")

sorted<-probes[order(probes$Chromosome,probes$Start <- loc) ,]

#Start Plot
pdf("a.pdf",width=11,height=8)

for (i in 1:nrow(probes)) {

  chromy <- paste("chr",probes$Chromosome[i],sep="")
  #Get the chromosome sequence
  seq <- Hsapiens[[chromy]]
  
  index=(probes$Start_loc[i]-1000):(probes$Finish_loc[i]+1000)
  #Subset the sequence
  subseq <- subseq(seq,start=index[1],end=index[length(index)])
  #Find CpGs
  cpgs=start(matchPattern("CG",subseq))+index[1]-1
  #Create the CpG bins
  cuts=seq(index[1],index[length(index)],8)
  #Bin the CpGs
  scpgs=cut(cpgs,cuts,include.lowest=TRUE)
  
  x = (cuts[1:(length(cuts)-1)]+cuts[2:(length(cuts))])/2
  y = table(scpgs)/8
  SPAN=400/diff(range(x))
  d = loess(y~x,span=SPAN,degree=1)

  plot(x,d$fitted,type="l",ylim=c(0,0.15),xlab="",xaxt="n",
           ylab="CpG density",xlim=c(index[1],index[length(index)]))
  rug(cpgs)
}
dev.off()
