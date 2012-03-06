##Run single ornstein-uhlenbeck simulation
ou.sim.single <- function (m0=1,mu=0.5, theta=1,sigma=.1,tfinal=10,deltat=0.05) {
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

ou.sim.monte <- function (nsimul=10, m0=0.5, mu=0.5, on.theta=1, off.theta=0,on.sigma=0.05, off.sigma=0.05,tfinal=20,
                          deltat=0.05, namey="ou.time1", draw=F) {

  ##Parameters - m0 is initial methylation, mu is the equilibrium spring ppoint
  ##theta is the magnitude of the restoring force, sigma is the magnitude of the noise
  ##tfinal is the end time of the simulation, deltat is the time step size
  ##nsimul is the number of times to run the simulation
  ##on.theta is the force magnitude of the spring when the spring is working (before cancer)
  ##off.theta is the force magnitude of the spring when the spring is broken (after cancer)
  ##on.sigma and off.sigma are the noise mangtidues before and after cancer genesis event
  ##namey is the name of the output file
  ##draw is a boolean asking whether to draw the results to a pdf or png

  ##Make mirrored time point, from -tfinal to positive tfinal, with 0 being the location of the cancer genesis event
  t=seq(from=-tfinal,to=tfinal,by=deltat)
  n=length(t)

  ##initialize array of methylation results
  runs=array(numeric(),dim=c(n,nsimul) )

  ##if draw, then make to a pdf, if not skip
  if (draw) {
    pdf(paste("Simul/", namey, ".pdf", sep=""), width=11, height=8)
 
    ##Initialize plot axes
    plot(t,t,type="n", ylim=c(0,1))
  }

    
  ##loop through the number of simulations to run
  for (j in 1:nsimul) {
    ##initialize methylation result
    y=numeric(n);
    ##initilize initila methylation value(given)
    y[1]=m0;

    ##Iterate through values
    for (i in 2:length(y)) {
      ##determine if we are before or after cancer genesis event
      theta=ifelse(t[i]>0,off.theta,on.theta)
      sigma=ifelse(t[i]>0,off.sigma,on.sigma)
      ##Euler-Maruyama
      y[i]=y[i-1]+deltat*(theta*(mu-y[i-1] ) )+sigma*sqrt(deltat)*rnorm(1);
      ##Methylation can't be greater than 1 or less than 0 - this may be a naive way to do it
      y[i]=max(c(0,y[i]))
      y[i]=min(c(1,y[i]))
    }
    ##Assign the results of iteration to one of the array
    runs[,j]=y
    
    if (draw) {
      ##Plot
      lines(t, y, type="l", lwd=1, col="red")
    }
  }

  if (draw) {
    ##Mark vertical lines at cancer genesis event and potential location of adenoma time point?
    abline(v=0, lty=2, col="green", lwd=1)
    abline(v=5, lty=2, col="blue", lwd=1)
  }
  
  ##legend("topright",c("Methylated", "Unmethylated"),col=c("green", "red"),lty=1,lwd=1) 

  data=list(t=t, simul=runs)
  return(data)

}

ou.bar<-function(norms, nsimul=10, off.theta=0, namey="progbar", on.sigma=0.05, off.sigma=0.05, on.theta=1) {
#Make boxplot of results at different time points - for real data comparison
  numprobes=length(norms)
  
  normy=array(numeric, dim=c(nsimul,numprobes))
  ady=array(numeric, dim=c(nsimul,numprobes))
  carcy=array(numeric, dim=c(nsimul,numprobes))


  for (i in 1:length(norms)) {
    res=ou.sim.monte(mu=norms[i], m0=norms[i], nsimul=nsimul, off.theta=off.theta, on.sigma=on.sigma, tfinal=30,
      off.sigma=off.sigma, on.theta=on.theta, namey="p")

    normy[,i]=res$simul[400,]
    ady[,i]=res$simul[res$t==3,]
    carcy[,i]=res$simul[res$t==30,]
    
  }

  pdf(paste("Simul/", namey, ".pdf", sep=""), width=8, height=8)
  ##Incomplete hereVVVV
  box
}
  
ou.density <- function(mu=0.5, off.theta=0, off.sigma=0.05, on.sigma=0.05,
                       on.theta=1, nsimul=1e3,
                       panel=F,loc=c(0.25, 0.25, 0.4, 0.4)) {
  
  ##Make density plot of results at different time points

  ##Do 100 simulations to make the density plots
  res=ou.sim.monte(mu=mu, m0=mu, nsimul=nsimul,
    off.theta=off.theta, on.sigma=on.sigma, on.theta=on.theta,
    tfinal=30)

  ##Normal values
  normy=res$simul[400,]
  normal=density(normy, from=0, to=1)

  ##Adenoma values
  ady=res$simul[res$t==3,]
  adenoma=density(ady,from=0, to=1)

  ##Carcinoma values
  carcy=res$simul[res$t==30,]
  carcinoma=density(carcy, from=0, to=1)

  if (panel) {
    vp=viewport(x = loc[1], y=loc[2], height = loc[3], width = loc[4],
      xscale=c(-0.05, 1.05), yscale=extendrange(c(carcinoma$y, adenoma$y, normal$y)))
    pushViewport(vp)

    grid.lines(normal$x, normal$y, gp=gpar(col="purple", lwd=2), default.units="native")
    panel.rug(normy, col="purple")
    grid.lines(adenoma$x, adenoma$y, gp=gpar(col="green", lwd=2), default.units="native")
    panel.rug(ady, col="green")
    grid.lines(carcinoma$x, carcinoma$y, gp=gpar(col="orange", lwd=2), default.units="native")
    panel.rug(normy, col="orange")


    grid.rect()
    grid.yaxis()
    grid.xaxis()

    popViewport()
               
  } else {
    plot(normal, type="l", col="green", lwd=2, xlim=c(0, 1))
    rug(normy, side=1,col="green")
    lines(adenoma, col="red", lwd=2)
    rug(ady, side=1, col="red")
    lines(carcinoma, col="blue", lwd=2)
    rug(carcy, side=1,col="blue")
  }
}


ou.hex <- function (mu=0.5, off.theta=0, off.sigma=0.05, on.sigma=0.05,
                    on.theta=1, namey="monte", nsimul=1e2, red=F,
                    loc=c(0.25, 0.5, 0.4, 0.4)) {
  require(hexbin)
  require(reshape)
  
  ##this function is designed to do ~1e3 runs of the oe simulation, and spit out a hexbin
  ##density plot of the result

  ##Do 100 simulations to make the density plots
  res=ou.sim.monte(mu=mu, m0=mu, nsimul=nsimul,
    off.theta=off.theta, on.sigma=on.sigma, on.theta=on.theta,
    tfinal=30)

  
  sim.vals=res$simul
  rownames(sim.vals)=res$t
  xy=melt(sim.vals)

  xbnds=extendrange(xy$X1)
  ybnds=c(-0.05, 1.05)
  
  meth.hex=hexbin(xy$X1, xy$value, xbnds=xbnds, ybnds=ybnds)

  pushViewport(vp=viewport(x = loc[1], y=loc[2],
                 height = loc[3], width = loc[4],
                 xscale=xbnds, yscale=ybnds))
  

  grid.hexagons(meth.hex)

  ##Red line traces
  if (red) {
    red.res=ou.sim.monte(mu=mu, m0=mu, nsimul=10,
      off.theta=off.theta, on.sigma=on.sigma, on.theta=on.theta,
      tfinal=30)    
    for (i in 1:10) {
      grid.lines(red.res$t, red.res$simul[,i], gp=gpar(col="#FF0000BF", lwd=0.4),
                 default.units="native")
    }
  }
  
  ##Mark vertical lines at cancer genesis event and potential location of adenoma time point?
  panel.abline(v=-10, lty=2, col="purple", lwd=2)
  panel.abline(v=3, lty=2, col="green", lwd=2)
  panel.abline(v=30, lty=2, col="orange", lwd=2)
  
  grid.rect()
  grid.yaxis()
  ##plot(meth.hex, newpage=F)
  popViewport()


  
}




