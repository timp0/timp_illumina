tissue_density <- function(avg, sampy, probes,which_tissue=c(2,3,4,6,7), divider=NA, interest=F,namey="a.pdf")  {
 
  ##Initialize with if not set to all possibles
  if (is.na(divider[1])) {
    divider=ifelse(sampy$Progression>0,0,1)
  }
  diffdiv=unique(divider)

  colory=factor(diffdiv)
  levels(colory)=coloringbook[1:nlevels(colory)]
  
  ##Start Plot
  pdf(paste("Movie/", namey, sep=""), width=11, height=8)

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
           
    ##Determine how many probes per region
    numprobes=nrow(not_far)
    ##How many different tissue types plotted
    numcomp=length(which_tissue)

    
   
    ##Setup that many plots
    ##layout(matrix(1:(numcomp*numprobes),nrow=numprobes,byrow=T),heights=c(rep(.2,numprobes)))
    bottom_mar=((.5+(7/3)*(3-numprobes)))
    bottom_mar=ifelse(bottom_mar<0, .5, bottom_mar)
    par(mfrow=c(numprobes,numcomp),mar=c(0,0,0,0),mgp=c(1.5,.5,0),omi=c(bottom_mar,.5,.5,.5) )

    
    inorder=not_far[order(not_far$Start_loc),]
    for (m in 1:numprobes) {
      probenum=as.numeric(row.names(inorder)[m])

      for (n in which_tissue) {    
        plot(0,0,ylim=c(0,7),xlim=c(-0,1),xlab="",type="n",axes="F")
        if (m==1) {
          mtext(text=classy[n],side=3,line=1.5,col="red")
        }
        
        box(lwd=2)
        ##Left Text
        if (n==head(which_tissue,1)) {
          axis(side=2,lwd=2,at=c(0,2,4,6))
          if (interest[probenum]) {
            mtext(text=probes$Probe_ID[probenum],side=2,line=1.5,col="red")
          } else {
            mtext(text=probes$Probe_ID[probenum],side=2,line=1.5,col="black")
          }
        }
        ##Bottom Text
        if (n==tail(which_tissue,1)) {
          mtext(text=probes$Illumina_ID[probenum],side=4,line=1.5)
        }
        
        for (o in 1:length(diffdiv)) {
          samples=(diffdiv[o]==divider)&(sampy$Class==n)
          ##Get all data points for this line
          rawnum<-avg[probenum,samples]
          ##Assign Color
          firstcolor=as.character(colory[o])
          ##
          d_data<-density(rawnum,from=0, to=1)
          lines(d_data,col=firstcolor,lwd=2)
          rug(rawnum,side=1,col=firstcolor)
        }
        
      }
      
    }
  }
  
  dev.off()
}



coloringbook=c("black","red","green","blue", "purple","orange","brown")
classy=c("Cross_Controls", "Breast", "Colon", "Lung", "Ovary", "Wilms", "Thyroid", "Pancreas")
