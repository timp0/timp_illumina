tissue_density <- function(avg, sampy, probes, interest=1:384,namey="a.pdf")  {

  divider=ifelse(sampy$Progression>0,1,4)
  divider[((sampy$Class==3)&(sampy$Progression>0)&(sampy$Progression<3))]=2
  divider[((sampy$Class==7)&(sampy$Progression>1)&(sampy$Progression<3))]=2

  
  ##Start Plot
  pdf(paste("Movie/", namey, sep=""), width=11, height=8)

  regions=numeric()
  
  ##Loop through all regions
  for (i in interest) {

    ##Find probes which are in the same region - this was previously defined by perl in the probes$Region command
    not_far<-probes[(probes$Region==probes$Region[i]),]
    num_close<-nrow(not_far)
           
    ##Determine how many probes per region
    numprobes=nrow(not_far)
    ##How many different tissue types plotted
    numcomp=5

    
   
    ##Setup that many plots

    bottom_mar=((.5+(7/3)*(3-numprobes)))
    bottom_mar=ifelse(bottom_mar<0, .5, bottom_mar)
    par(mfrow=c(numprobes,numcomp),mar=c(0,0,0,0),mgp=c(1.5,.5,0),omi=c(bottom_mar,.5,.5,.5) )

    
    inorder=not_far[order(not_far$Start_loc),]
    for (m in 1:numprobes) {
      probenum=as.numeric(row.names(inorder)[m])

      for (n in c(2,3,4,6,7)) {    
        plot(0,0,ylim=c(0,7),xlim=c(-0,1),xlab="",type="n",axes="F")
        if (m==1) {
          mtext(text=classy[n],side=3,line=1.5,col="red")
        }
        
        box(lwd=2)
        ##Bottom Text
        if (m==numprobes) {
          axis(side=1, lwd=2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
        }
        ##Left Text
        if (n==2) {
          axis(side=2,lwd=2,at=c(0,2,4,6))
          if (probenum==i) {
            mtext(text=probes$Probe_ID[probenum],side=2,line=1.5,col="red")
          } else {
            mtext(text=probes$Probe_ID[probenum],side=2,line=1.5,col="black")
          }
        }
        ##Right Text
        if (n==7) {
          ##mtext(text=probes$Illumina_ID[probenum],side=4,line=1.5)
          if (data$probes$UCSC_Dist_to_Island[probenum] == 0) {
            mtext(text="Island", side=4,line=1.5)
          } else {
            if (data$probes$UCSC_Dist_to_Island[probenum] < 2000) {
              mtext(text="Shores", side=4, line=1.5)
            } else {
              mtext(text="Far", side=4, line=1.5)
            }
          }
        }
        
        
        
        for (o in c(1,2,4)) {
          samples=(o==divider)&(sampy$Class==n)
          if (any(samples)) {
            ##Get all data points for this line
            rawnum<-avg[probenum,samples]
            ##Assign Color
            if (o==2) {
              firstcolor="orange"
            } else {
              firstcolor=paste(coloringbook[n],as.character(o))
            }
            ##
            d_data<-density(rawnum,from=0, to=1)
            lines(d_data,col=firstcolor,lwd=2)
            rug(rawnum,side=1,col=firstcolor)
          }
        }
        
      }
      
    }
  
  }
  dev.off()
}


col_thy_density <- function(avg, sampy, probes, interest=1:384,namey="a.pdf")  {

  divider=ifelse(sampy$Progression>0,1,4)
  divider[((sampy$Class==3)&(sampy$Progression>0)&(sampy$Progression<3))]=2
  divider[((sampy$Class==7)&(sampy$Progression>1)&(sampy$Progression<3))]=2

  
  ##Start Plot
  pdf(paste("Movie/", namey, sep=""), width=8, height=8)

  par(mfrow=c(2,2) )
  
  
  for (m in interest){
    
    for (n in c(3,7)) {    

      #Get points
      normy=avg[m,((1==divider)&(sampy$Class==n))]
      normal=density(normy, from=0, to=1)
      ady=avg[m,((2==divider)&(sampy$Class==n))]
      adenoma=density(ady, from=0, to=1)
      carcy=avg[m,((4==divider)&(sampy$Class==n))]
      carcinoma=density(carcy, from=0, to=1)

      #Do plots
      plot(normal, type="l", col="green", lwd=2, xlim=c(0, 1), ylim=c(0,7))
      rug(normy, side=1,col="green")
      lines(adenoma, col="red", lwd=2)
      rug(ady, side=1, col="red")
      lines(carcinoma, col="blue", lwd=2)
      rug(carcy, side=1,col="blue")


      
      mtext(text=classy[n],side=3,line=1.5,col="red")
      
      mtext(text=probes$Probe_ID[m],side=2,line=1.5,col="red")

      if (data$probes$UCSC_Dist_to_Island[m] == 0) {  
        mtext(text="Island", side=4,line=1.5)
      } else {
        if (data$probes$UCSC_Dist_to_Island[m] < 2000) {
          mtext(text="Shores", side=4, line=1.5)
        } else {
          mtext(text="Far", side=4, line=1.5)
        }
      }

    }
  }

  dev.off()
}




coloringbook=c("black","coral","cadetblue","seagreen", "black","mediumorchid","royalblue")
classy=c("Cross_Controls", "Breast", "Colon", "Lung", "Ovary", "Wilms", "Thyroid", "Pancreas")
