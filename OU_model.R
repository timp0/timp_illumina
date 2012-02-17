##Run single ornstein-uhlenbeck simulation
ou_sim1 <- function (m0=1,mu=0.5, theta=1,sigma=.1,tfinal=10,deltat=0.05) {
  ##Parameters - m0 is initial methylation, mu is the equilibrium spring ppoint
  ##theta is the magnitude of the restoring force, sigma is the magnitude of the noise
  ##tfinal is the end time of the simulation, deltat is the time step size

  ##Create all time points
  t=seq(from=0,to=tfinal,by=deltat)
  ##number of time points
  n=length(t)
  ##Initialize methylation vector
  y=numeric(n);

  ##Set initial point
  y[1]=m0;

  ##Iterate through time
  for (i in 2:length(y)) {
    ##Euler-Maruyama method
    y[i]=y[i-1]+deltat*(theta*(mu-y[i-1] ) )+sigma*sqrt(deltat)*rnorm(1);
    ##Can't be less than 0 methylation
    y[i]=max(c(0,y[i]))
    ##Can't be greater than 1 methylation
    y[i]=min(c(1,y[i]))
  }

  ##Set up a data frame with time and methylation
  m=data.frame(t=t, y=y)
  return(m)
 
}

ou_sim2 <- function (nsimul=10, m0=0.5, mu=0.5, on_theta=1, off_theta=0,on_sigma=0.05, off_sigma=0.05,tfinal=20, deltat=0.05, namey="ou_time1", draw=F) {
  ##Parameters - m0 is initial methylation, mu is the equilibrium spring ppoint
  ##theta is the magnitude of the restoring force, sigma is the magnitude of the noise
  ##tfinal is the end time of the simulation, deltat is the time step size
  ##nsimul is the number of times to run the simulation
  ##on_theta is the force magnitude of the spring when the spring is working (before cancer)
  ##off_theta is the force magnitude of the spring when the spring is broken (after cancer)
  ##on_sigma and off_sigma are the noise mangtidues before and after cancer genesis event
  ##namey is the name of the output file
  ##draw is a boolean asking whether to draw the results to a pdf or png

  ##Make mirrored time point, from -tfinal to positive tfinal, with 0 being the location of the cancer genesis event
  t=seq(from=-tfinal,to=tfinal,by=deltat)
  n=length(t)

  ##initialize array of methylation results
  runs=array(numeric(),dim=c(n,nsimul) )

  ##if draw, then make to a pdf, if not to a png(png is smaller)
  if (draw==T) {
    pdf(paste("Simul/", namey, ".pdf", sep=""), width=11, height=8)
  } else {
    png(paste("Simul/", namey, ".png", sep=""))
  }

  ##Initialize plot axes
  plot(t,t,type="n", ylim=c(0,1))

  ##loop through the number of simulations to run
  for (j in 1:nsimul) {
    ##initialize methylation result
    y=numeric(n);
    ##initilize initila methylation value(given)
    y[1]=m0;

    ##Iterate through values
    for (i in 2:length(y)) {
      ##determine if we are before or after cancer genesis event
      theta=ifelse(t[i]>0,off_theta,on_theta)
      sigma=ifelse(t[i]>0,off_sigma,on_sigma)
      ##Euler-Maruyama
      y[i]=y[i-1]+deltat*(theta*(mu-y[i-1] ) )+sigma*sqrt(deltat)*rnorm(1);
      ##Methylation can't be greater than 1 or less than 0 - this may be a naive way to do it
      y[i]=max(c(0,y[i]))
      y[i]=min(c(1,y[i]))
    }
    
    runs[,j]=y
    ##Plot
    lines(t, y, type="l", lwd=1, col="red")
  }

  abline(v=0, lty=2, col="green", lwd=1)
  abline(v=5, lty=2, col="blue", lwd=1)
  
  ##legend("topright",c("Methylated", "Unmethylated"),col=c("green", "red"),lty=1,lwd=1) 
  dev.off()
  data=list(t=t, simul=runs)
  return(data)

}



sense_test <- function() {
  
  sensetheta=seq(from=0, to=1, by=0.05)
  sensesigma=seq(from=0, to=1, by=0.05)
  

  for (i in sensetheta) {
    for (j in sensesigma) {

      
      z=mult_sim1(theta=i, sigma=j, mu=1,nsimul=1e2,namey=paste(i, j, "trace.png", sep="_"))
      d=density(z[201,], from=0, to=1)
      png(paste("Simul/",paste(i,j,"dense.png",sep="_"), sep=""))
      plot(d)
      dev.off()
      
          
      
    }
  }
}

ou_bar<-function(norms, nsimul=10, off_theta=0, namey="progbar", on_sigma=0.05, off_sigma=0.05, on_theta=1) {

  numprobes=length(norms)
  
  normy=array(numeric, dim=c(nsimul,numprobes))
  ady=array(numeric, dim=c(nsimul,numprobes))
  carcy=array(numeric, dim=c(nsimul,numprobes))


  for (i in 1:length(norms)) {
    res=ou_sim2(mu=norms[i], m0=norms[i], nsimul=nsimul, off_theta=off_theta, on_sigma=on_sigma, tfinal=30, off_sigma=off_sigma, on_theta=on_theta, namey="p")

    normy[,i]=res$simul[400,]
    ady[,i]=res$simul[res$t==3,]
    carcy[,i]=res$simul[res$t==30,]
    
  }

  pdf(paste("Simul/", namey, ".pdf", sep=""), width=8, height=8)

  box
}
  
ou_density<- function(mu=0.5, off_theta=0, namey="progress", on_sigma=0.05, off_sigma=0.05, on_theta=1) {

  
  res=ou_sim2(mu=mu, off_theta=off_theta, nsimul=10, m0=mu, on_sigma=on_sigma, off_sigma=off_sigma, on_theta=on_theta,namey=paste(namey, "trace",sep="_"), draw=T)
  
  res=ou_sim2(mu=mu, m0=mu, nsimul=100, off_theta=off_theta, on_sigma=on_sigma, tfinal=30, off_sigma=off_sigma, on_theta=on_theta,namey="p")

  ##Normal values
  normy=res$simul[400,]
  normal=density(normy, from=0, to=1)

  ##Adenoma values
  ady=res$simul[res$t==3,]
  adenoma=density(ady,from=0, to=1)

  ##Carcinoma values
  carcy=res$simul[res$t==30,]
  carcinoma=density(carcy, from=0, to=1)

  pdf(paste("Simul/", namey, ".pdf" , sep=""), width=8, height=8)
  par(mfrow=c(2,2))
  plot(normal, type="l", col="green", lwd=2, xlim=c(0, 1))
  rug(normy, side=1,col="green")
  lines(adenoma, col="red", lwd=2)
  rug(ady, side=1, col="red")
  lines(carcinoma, col="blue", lwd=2)
  rug(carcy, side=1,col="blue")
  dev.off()

}



