##Old Illumina sample plotting code

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
plotsampreg <- function(class_stuff,probes,avg_data,sampy,not_far,region) {
  
  for (j in 1:length(class_stuff$nums)) {
    data<-avg_data[(probes$Region==region),(sampy$Class==class_stuff$nums[j])]
    
    samp_info<-sampy[sampy$Class==class_stuff$nums[j],]
    ##cat(class_stuff$nums[j],"\n")
    num_samp=ncol(data)
    num_pts=nrow(not_far)
    for (k in 1:num_samp) {
      lines(not_far$Start_loc,data[,k],col=class_stuff$coloring[j])
      
    }
  }
  legend("topright",class_stuff$nam,col=class_stuff$coloring,pch=1)
  
}


plotsampregdelta <- function(class_stuff,probes,avg_data,sampy,not_far,region) {
  
  for (j in 1:length(class_stuff$nums)) {
    data<-avg_data[(probes$Region==region),(sampy$Class==class_stuff$nums[j])]
    samp_info<-sampy[sampy$Class==class_stuff$nums[j],]
    ##cat(class_stuff$nums[j],"\n")
    num_samp=ncol(data)
    num_pts=nrow(not_far)
    if (num_pts>1) {
      delta=apply(data,2,diff)
      locs=not_far$Start_loc[1:num_pts-1]+diff(not_far$Start_loc)/2
      for (k in 1:num_samp) {
        if (length(locs)==1) {
          points(locs,delta[k],col=class_stuff$coloring[j])
        } else {
          lines(locs,delta[,k],col=class_stuff$coloring[j])
        }
      }
    }
  }
  
  legend("topright",class_stuff$nam,col=class_stuff$coloring,pch=1)
  
}
  


plotClass <- function(probes,genes,avg_data,sampy,which,outname,controls=TRUE,interest=logical(1)) {
  ##
  ##Start Plot
  pdf(outname,width=11,height=8)
  par(mfrow=c(3,1))
  layout(c(1:3),heights=c(0.5,0.3,0.2))
  par(mar=c(2,1,0,1),mgp=c(1.5,.5,0),oma=c(0,1,2,1)) 
  


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
    index=(min(not_far$Start_loc)-1000):(max(not_far$Finish_loc)+1000)
    final_index=length(index)
    plot(0,0,ylim=c(-1,1.1),xlim=c(index[1]+900,index[final_index]-900),ylab="Methylation",xlab="",type="n",axes=F)
    box()
    axis(2, at=seq(0,1,.2))

    ##Plot actual samples
    plotsampreg(class_stuff[which,],probes,avg_data,sampy,not_far,i)
    ##Plot label on axis as a tick on the bottom
    axis(side=1,at=not_far$Start_loc,labels=not_far$Probe_ID, cex.axis=0.8)
    
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
