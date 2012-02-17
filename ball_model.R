ball_plot1 <- function (m=0,mu=0.5) {
  
  r=0.05
  
  x=seq(from=0, to=1, by=0.001)
  
  y=(x-mu)^2

  a=1/max(y)
  
  theta=atan(2*a*(x-mu))-pi/2
  
  x2=x+r*cos(theta)
  y2=a*(x-mu)^2+r*sin(theta)
  l=length(x2)
  #z=ou_sim2(deltat=0.25)
  
  #library(plotrix)
  
  plot(x2,y2,type="l",lwd=2, yaxt="n", ylab="",asp=1, yaxt="n", frame.plot="F", xlim=c(-.5,1.05))

  abline(v=-0.05, lty=2,lwd=2,col="red")
  abline(v=1.05, lty=2,lwd=2,col="red")

  draw.circle(m, a*(m-mu)^2, r, col="red")
       
}

spring_plot1 <- function(m=0,mu=0.5) {  

  #Spring attachment point is -(1-mu)
  e=-1.2+mu
  s=0.05
  h=0.05
  
  plot(c(e,e+s), c(0,0), type="l", lwd=2, xlim=c(-.5, 1.05), asp=1, axes="F", frame.plot="F")

  lines(c(e, e), c(-1e3, 1e3), lwd=2)

  spring_x=seq(from=e+s, to=m-s, length=20)
  spring_y=c(0,h,-h,h,-h,h,-h,h,-h,h,-h,h,-h,h,-h,h,-h,h,-h,0)
  
  lines(spring_x, spring_y, lwd=2)
  lines(c(m-s, m), c(0,0), lwd=2)
  draw.circle(m,0,h, col='red')

  abline(v=-0.05, lty=2,lwd=2,col="red")
  abline(v=1.05, lty=2,lwd=2,col="red")
  

}



spring_para_plot1 <- function(x=0,mu=0.5) {

  layout(c(1,2), heights=c(0.2,0.8))
  par(mar=c(1,0.2,1,0.2))

  spring_plot1(m=x,mu=mu)
  ball_plot1(m=x,mu=mu)
}


under_harmonic1 <- function (gamma=0.1) {


  t=seq(from=0, to=3e1, length=2e2)

  x=-0.5*exp(-gamma*t)*cos( (1-gamma^2)^(1/2)*t)+0.5

  for (i in 1:length(t)) {
    png( paste("Ball/u", i, ".png", sep=""))

    spring_para_plot1(x[i])
    dev.off()
  }

}

over_harmonic1 <- function () {

  gamma=1.1
  c1=-(gamma - (gamma^2 -1)^(1/2))/(4*(gamma^2-1)^(1/2))
  c2=-(c1+.5)

  t=seq(from=0, to=3e1, length=2e2)

  x=exp(-gamma*t)*( c1 * exp(t) + c2 * exp(-t) ) +0.5

  for (i in 1:length(t)) {
    png( paste("Ball/o", i, ".png", sep=""))

    spring_para_plot1(x[i])
    dev.off()
  }

}

brown_harmonic1 <- function () {
  
  z=ou_sim2(tfinal=3e1, deltat=3e1/2e2,m0=0,on_theta=1/1.1, off_theta=1/1.1)  

  for (i in 1:length(z$t)) {
    png( paste("Ball/b", i, ".png", sep=""))

    spring_para_plot1(z$simul[i,1])
    dev.off()
  }
}

brown_harmonic2 <- function () {
  
  z=ou_sim2(tfinal=3e1, deltat=3e1/2e2,m0=0,on_theta=1/1.1, off_theta=1/1.1, mu=1,draw=T)  

  for (i in 1:length(z$t)) {
    png( paste("Ball/b", i, ".png", sep=""))

    spring_para_plot1(z$simul[i,1], mu=1)
    dev.off()
  }
}


library(plotrix)
