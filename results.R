#This is just results.R, but I've customized it for winston's purpose.
#So this is not a file for general use.

nummy <- 20

#May have to use ~ririzarr/bin/R-2.5.0 if BSgenome species genome fails.
#If so, need to use the contemporaneous versions of all the packages:
.libPaths("/home/bst/student/bcarvalh/bin/R-2.5.0/library")

#Rlibspath <- "~/Rlibs" #location of BSgenome packages, if you 
#                       #installed them yourself.  If not, say NULL.
inpath  <- "~wtimp/cancer_dmr/From_Rafa"
  #path to where your sMM.rda, M.rda, otherstuff.rda, 
  #and object.rda charm objects are, which must exist.
outpath <- "~/thumper/Work/Notes/Genomics/Cancer DMRs/Regions"
  #path to where you want your output files to go.
numplots <- 300 #number of plots you want. e.g., if numplots=100, then the 
                #top 100 regions will be plotted in the plot.pdf file.
spec <- "human" #cis-genome.R currently only implements "human" or "mouse"
#Obj.i <- c(1,2,3) #which data frames in object.rda you want to care about.  
        #i.e., which rows in COMPS do you care about.
all.lines <- FALSE #Only relevant if there are >2 exposure groups.
    #Do you want to plot smooth lines for not just the
    #two groups being compared, but for all groups?  If TRUE, note 
    #that the regions in each plot may have >2 lines, but those regions
    #were identified only as DMRs between the first 2 groups listed in 
    #the legend, not the others.

#########Load in your charm objects:###########################
setwd(inpath)
load("otherstuff.rda") #need CHRLEVELS, pos, and pns
if(!exists("sMM")) load("sMM.rda") #need FACS, INDEX, and sMM
if(!exists("M")){
  load("M.rda")
  library(genefilter)
  tmpmed=rowMedians(t(M))
  M=sweep(M,2,tmpmed,"/")
  M=M*median(tmpmed)
}
#load("object.rda") #need NAMES and object
obj <- read.csv("~/thumper/Work/Notes/Genomics/Cancer DMRs/Regions/Top_try.csv")

###############################################################
####Load in Rafa's objects (hopefully he won't move them):#####
setwd("~ririzarr/projects/cegs/where-is-methylation/forgrant")
mypar <-function(a=1,b=1,brewer.n=8,brewer.name="Dark2",...){
  require(RColorBrewer)
  par(mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0)) 
  par(mfrow=c(a,b),...)
  palette(brewer.pal(brewer.n,brewer.name)) 
}
#source("../functions.R") #bad path to fm.R inside, and matchGenes()
#needs to modified to handle different species, and functions2mouse.R
#needs to have every "mm#" be mm8, so replace with:
if(spec=="human") source("~pmurakam/feinberg/CharmFiles/functions.R") 

#source("../cis-genome.R") #handles both human and mouse
source("~pmurakam/feinberg/CharmFiles/cis-genome.R") #changed paths 
  #inside for human so need not be on node 36, and changed mouse 
  #version mm9 to mm8, since arrays were designed with mm8.
  
if(exists("Rlibspath")) Rlibspath <- c(Rlibspath,.libPaths())
#library(BSgenome,lib.loc=Rlibspath)
    #newer R needed for newer versions of BSgenome.
library(BSgenome) #if using Rafa's R-2.5.0
if(spec=="human"){
    load("../../islands/current-cpg.rda")
#    library(BSgenome.Hsapiens.UCSC.hg18,lib.loc=Rlibspath)
    library(BSgenome.Hsapiens.UCSC.hg18)

#    CTCF=read.delim("~ririzarr/data/cegs/where-is-methylation/CTCFHG18.txt",header=FALSE,as.is=TRUE)
#    names(CTCF)=c("chr","start","end")
#    load("~ririzarr/projects/cegs/where-is-methylation/rdas/bivalent.rda")
}

ocpgi = data.frame(chr=I(cpg.cur[,1]),start=as.numeric(cpg.cur[,2]),
                   end=as.numeric(cpg.cur[,3]))
ocpgi=ocpgi[ocpgi$chr%in%CHRLEVELS,]

library("genefilter")
setwd(outpath) #some of matchGenes()' intermediate files will go here too.

NAMES=names(INDEX)  #INDEX is from sMM.rda. 
  #(We already have NAMES from object.rda actually, so this is unnecessary.)
Indexes=split(seq(along=pns),pns) #pns is from otherstuff.rda
                                  #pns=probeNames(Data).
   
#################Make output table and pdf file of plots##############
#for(object.i in Obj.i){
    #################Make output table:###############################
#    obj=object[[object.i]]
#    absArea <- abs(obj$diff)*obj$npr
#    obj$absArea <- absArea
#    obj <- obj[order(-absArea),]
#    
#    tmp   <- matchGenes(obj,species=spec)
#    obj$gene           <- tmp$name
#    obj$dist2gene      <- tmp$dist
#    obj$relation2gene  <- tmp$description
#    islas <- regionMatch(obj,ocpgi)
#    obj$dist2island     <- islas$dist
#    obj$relation2island <- islas$type
#    if(spec=="human"){
#        biva  <- regionMatch(obj,bivalent)
#        obj$dist2bivalent     <- biva$dist
#        obj$relation2bivalent <- biva$type
#        ctcf  <- regionMatch(obj,CTCF)
#        obj$dist2ctcf     <- ctcf$dist
#        obj$relation2ctcf <- ctcf$type
#    }

#   write.csv(obj,file=paste("table",object.i,".csv",sep=""), row.names=FALSE,quote=FALSE)
    
    ################Make pdf file of plots###########################
    pdf("topprobes1.pdf",width=11,height=8,pointsize=10)
    mypar()
    palette(rev(brewer.pal(8,"Dark2")[c(2,1,3,4,8)]))
    ADD1=1500
    ###colors for plots
    
    #comps <- c(match(colnames(obj)[8],NAMES),match(colnames(obj)[9],NAMES)) 
        #I added this line.  use for col inds for sMM below.
    comps <- c(4,5)
    tcols=as.numeric(factor(FACS,levels=NAMES[comps])) #used only in matplot
    if(all.lines==TRUE) comps <- c(comps,setdiff(1:length(NAMES),comps))

    for(i in 1:min(nrow(obj),numplots)){
      YLIM=c(-0.5,3)
      thechr=obj$chr[i]
      ADD=max(ADD1*2,obj$end[i]-obj$start[i]+132)
      start=floor((obj$start[i]+obj$end[i])/2)-ADD
      end=start+2*ADD
    
      Index=Indexes[[obj$index[i]]] #index is what probeset in pns its in
      x=pos[Index] #pos is from otherstuff.rda.  pos=pmPosition(Data).
      start=max(start,min(x))
      end=min(end,max(x))
      
      Index= Index[x>=start & x <= end]
      Index=Index[order(pos[Index])]
      
      ##PLOT TISSUES
      par(mfrow=c(3,1))
      layout(matrix(1:3,ncol=1),heights=c(0.6,0.2,0.2))
      par(mar=c(0,2.5,0.25,1.1),oma=c(0,0,2,0))
      matplot(pos[Index],M[Index,],pch=1,col=tcols,cex=0.6,ylim=YLIM,xlab="",
              xaxt="n",ylab="M",xlim=c(start,end),lty=1)
      
      matlines(pos[Index],sMM[Index,comps],lwd=2,lty=1) #colors cycle
      
      abline(h=0,lwd=2,lty=3)
      abline(h=1,lwd=2,lty=3)
      
      legend("topright",NAMES[comps],col=1:length(comps),lty=1,lwd=2)
      #I add this to just this file:
      #This shows the location of the top probes

      #From Rafa book on R(Gentleman et al, R and Bioconductor) p. 35
      #Set a colormap

      rainbow <- brewer.pal(nummy, "Spectral")
      
      for(f in 0:nummy-1){
        probey=obj[i,f*2+12]
        toppy=1+2*(f%%2)
        if(!identical(probey, NA)) rug(c(probey,probey+50),col=rainbow[(f%%10)+1],lwd=3,side=toppy)
      }
      
      
      ##PLOT CPG ISLANDS
      if(spec=="human") seq<-Hsapiens[[as.character(thechr) ]]
      if(spec=="mouse") seq<-Mmusculus[[as.character(thechr) ]]
      Index=max(start,1):min(end,length(seq))
      subseq<-seq[Index]  #using older version of Biostrings. For newer:
    #  subseq<-subseq(seq,start=Index[1],end=Index[length(Index)])
      cpgs=start(matchPattern("CG",subseq))+Index[1]-1
      ##mcrbc=sort(c(start(matchPattern("ACG",subseq))+Index[1]-1,
      ##  start(matchPattern("GCG",subseq))+Index[1]-1))
      cuts=seq(min(Index),max(Index),8)
      scpgs=cut(cpgs,cuts,include.lowest=TRUE)
      x = (cuts[1:(length(cuts)-1)]+cuts[2:(length(cuts))])/2
      y = table(scpgs)/8
      SPAN=400/diff(range(x))
      d = loess(y~x,span=SPAN,degree=1)
      plot(x,d$fitted,type="l",ylim=c(0,0.15),xlab="",xaxt="n",
           ylab="CpG density",xlim=c(start,end))
      rug(cpgs)
    ##  rug(mcrbc,side=3)
      Index1 = which(ocpgi[,1]==as.character(thechr) &
                     ((ocpgi[,2] > start & ocpgi[,2]< end) |
                      (ocpgi[,3] > start & ocpgi[,3]< end)))
      if(length(Index1)>0)
          sapply(Index1,function(j) rug(unlist(ocpgi[j,2:3]),
                 col=5,lwd=3,side=1))
          #WHY WARNINGS? cpg.cur NOT CORRECT? Rafa said it happens.
      
      ##PLOT GENES
      par(mar=c(2.5,2.5,0.25,1.1))
      print(i)
      if(TRUE){
        tmp=neighborgenes(data.frame(chr=thechr,start=start,end=end),
          g=ADD,up=10,down=10,species=spec)
        tmp=tmp[!duplicated(tmp[,6]),]
        tmp=tmp[,c(2,14,15)]
        names(tmp)<- c("chr","start","end")
        tmp[is.na(tmp[,2]),2] <- "---"
        tmp=tmp[tmp[,2]!="---",]
        plot(0,0,ylim=c(-1.5,1.5),xlim=c(start,end),yaxt="n",ylab="Genes",
             xlab="Location")
        axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1)
        abline(h=0,lty=3)
        if(length(tmp[,1])>0){
          genes=nearestgene(tmp,species=spec)
          for(l in 1:nrow(tmp)){
            TS = genes[l,10]
            TE = genes[l,11]
            CS = genes[l,12]
            CE = genes[l,13]
            ES = as.numeric(strsplit(genes[l,15],",")[[1]])
            EE = as.numeric(strsplit(genes[l,16],",")[[1]])
            Exons= cbind(ES,EE)
            Strand= ifelse(genes[l,9]=="+",1,-1)
            polygon(c(TS,TE,TE,TS),Strand/2+c(-0.4,-0.4,0.4,0.4),density=0,col=l)
            polygon(c(CS,CE,CE,CS),Strand/2+c(-0.25,-0.25,0.25,0.25),density=0,
                    col=l) 
            apply(Exons,1,function(x) polygon(c(x[1],x[2],x[2],x[1]),
                                              Strand/2+c(-0.2,-0.2,0.2,0.2),col=l))
            
            text((max(start,TS)+min(end,TE))/2,Strand*0.9,genes[l,6],cex=1,
                  pos=Strand+2)
            ##Strand +2 works cuase -1 turns into 1 (below) and 1 
            ##turns into 3 (above)
          }
        }
        mtext(paste("ID:",i,"--",as.character(thechr),":",start,"-",end,sep=""),
              side=3,cex=1.5,outer=TRUE)
      }
    }
    dev.off()

#}

cat("\n results finished.\n")
