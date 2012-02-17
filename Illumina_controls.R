plotcontrols <- function(avg_data,sampy,not_far,region,start,end) {

  ##Plot actual probe locations
  ##Control_samples

  cntrl_data<-avg_data[(probes$Region==region),(sampy$Class==1)]
  cntrl_info<-sampy[sampy$Class==1,]
  
  plot(0,0,ylim=c(0,1),xlim=c(start,end),"axes"=FALSE, type="n",xlab="",ylab="")
  axis(side=4,at=c(0,0.5,1))

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
  axis(side=3,at=not_far$Start_loc,labels=not_far$Probe_ID, cex.axis=0.8)
  axis(side=1)
}
