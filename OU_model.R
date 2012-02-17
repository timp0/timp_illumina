ou_sim1 <- function (m0=1,mu=0.5, theta=1,sigma=.1,tfinal=10,deltat=0.05) {

  t=seq(from=0,to=tfinal,by=deltat)
  n=length(t)
    
  y=numeric(n);

  y[1]=m0;
  
  for (i in 2:length(y)) {
    y[i]=y[i-1]+deltat*(theta*(mu-y[i-1] ) )+sigma*sqrt(deltat)*rnorm(1);
    y[i]=max(c(0,y[i]))
    y[i]=min(c(1,y[i]))
  }

  m=data.frame(t=t, y=y)
  return(m)
 
}

ou_sim2 <- function (nsimul=10, m0=0.5, mu=0.5, on_theta=1, off_theta=0,sigma=0.05, tfinal=20, deltat=0.05, namey="ou_time1.png") {

  t=seq(from=-tfinal,to=tfinal,by=deltat)
  n=length(t)    
  runs=array(numeric(),dim=c(n,nsimul) )

  png(paste("Simul/", namey, sep=""))
  plot(t,t,type="n", ylim=c(0,1))
  
  for (j in 1:nsimul) {
    y=numeric(n);
    
    y[1]=m0;

    for (i in 2:length(y)) {
      theta=ifelse(t[i]>0,off_theta,on_theta)
      y[i]=y[i-1]+deltat*(theta*(mu-y[i-1] ) )+sigma*sqrt(deltat)*rnorm(1);
      y[i]=max(c(0,y[i]))
      y[i]=min(c(1,y[i]))
    }
    
    runs[,j]=y
    ##Plot
    lines(t, y, type="l", lwd=1, col="red")
  }
  
  ##legend("topright",c("Methylated", "Unmethylated"),col=c("green", "red"),lty=1,lwd=1) 
  dev.off()
  return(runs)

}




mult_sim1 <- function (nsimul=10,m0=1,mu=0.5, theta=1,sigma=.1,tfinal=10,deltat=0.05, namey="simul1.png") {

  runs=array(numeric(),dim=c(tfinal/deltat+1,nsimul) )

  t=seq(from=0,to=tfinal,by=deltat)

  png(paste("Simul/", namey, sep=""))

  plot(t,t,type="n", ylim=c(0,1))
  
  for (i in 1:nsimul) {
    z=ou_sim1(m0=m0,mu=mu,theta=theta,sigma=sigma, tfinal=tfinal, deltat=deltat) 
    ##Get approximations
    runs[,i]=z$y
    ##Plot
    lines(z$t, z$y, type="l", lwd=1, col="red")
  }
  
  ##legend("topright",c("Methylated", "Unmethylated"),col=c("green", "red"),lty=1,lwd=1) 
  dev.off()
  return(runs)
  
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

