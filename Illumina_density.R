plotdsamples <- function(class_stuff,probes,avg_data,sampy,not_far,region,highlight) {
  ##Determine how many probes per region
  numprobes=nrow(not_far)
  ##How many sample types
  numsamp=4

  ##Setup that many plots
  
  layout(rbind(matrix(1:(numsamp*numprobes),nrow=numprobes,byrow=T), rep((numprobes*numsamp+1),numsamp), rep((numprobes*numsamp+2),numsamp)),heights=c(rep((0.7/numprobes),numprobes),0.15,0.15))
  ##layout(c(1:3),heights=c(0.5,0.3,0.2))
  par(mar=c(0,0,0,0),mgp=c(1.5,.5,0),oma=c(1.5,3,2,3))
  
  inorder=not_far[order(not_far$Start_loc),]

  just=c(1,2,3,5)
  
  for (m in 1:numprobes) {
    probenum=as.numeric(row.names(inorder)[m])

    for (n in just) {
      plot(0,0,ylim=c(0,7),xlim=c(-0,1),xlab="",type="n",axes="F")
      box(lwd=2)
      if (n==1) {
        axis(side=2,lwd=2,at=c(0,2,4,6))
        if (highlight[probenum]) {
          mtext(text=probes$Probe_ID[probenum],side=2,line=1.5,col="red")
        } else {
          mtext(text=probes$Probe_ID[probenum],side=2,line=1.5,col="black")
        }
      }

      if (n==tail(just,1)) {
        mtext(text=probes$Illumina_ID[probenum],side=4,line=1.5)
      }
      
      if (m==numprobes) {
        axis(side=1,lwd=2,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1"))
        mtext(text=class_stuff$nam[n*2],side=1,line=1.5)
      }
      
      
      rawnum<-as.numeric(avg_data[probenum,(sampy$Class==2*n)])
      d_data<-density(rawnum,from=0, to=1)
      lines(d_data,col=class_stuff$coloring[2*n],lwd=2)
      rug(rawnum,side=1,col=class_stuff$coloring[2*n])

      rawnum<-as.numeric(avg_data[probenum,(sampy$Class==(2*n+1))])
      d_data<-density(rawnum,from=0, to=1)
      lines(d_data,col=class_stuff$coloring[2*n+1],lwd=2)
      rug(rawnum,side=3,col=class_stuff$coloring[2*n+1])

        
    }

  }
  par(mar=c(0,1,4,1))
}
  
plotsubdTN <- function(probes,genes,avg_data,sampy,outname,controls=TRUE,interest=logical(1)) {
  ##
  ##Start Plot
  pdf(outname,width=11,height=8)

  regions=numeric()
  
  ##Identify good regions
  if (!any(interest)) {
    interest=logical(length(probes$Region))
    interest[]=T
  }

  
  ##Loop through all regions
  for (i in unique(probes$Region[interest])) {

  
  ##Find probes which are in the same region - this was previously defined by perl in the probes$Region command
    not_far<-probes[(probes$Region==i),]
    num_close<-nrow(not_far)
    
    ##Set which chromosome and coordinates we are using for the region
    chromy <- paste("chr",not_far$Chromosome[1],sep="")
    ##cat("Region ", i, "\n")
    index=(min(not_far$Start_loc)-1000):(max(not_far$Finish_loc)+1000)
    final_index=length(index)

    plotdsamples(class_stuff,probes, avg_data,sampy,not_far,i,interest)

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
    par(mar=c(0,1,2,1))    
    ##Plot RefSeq Gene regions
    plotgenes(genes,chromy, index[1], index[final_index])
  }
  dev.off()
}

