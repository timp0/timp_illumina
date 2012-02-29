ball.plot <- function (m=0,mu=0.5, a=4) {
  ##Make ball in parabolic potential well plot, using m as location, and mu=0.5 as the bottom of potential well
  ##radius of circle is 0.05
  r=0.05

  ##X coordinate of potential
  x=seq(from=0, to=1, by=0.001)
  ##Y coordinate of potential
  y=a*(x-mu)^2

  ##Max value of potential
  ##a=1/max(y)

  ##Get angle for drawing the edge of potential - this is complicated, but its needed to
  ##draw a line so that the ball(falling the parabola), touches in one and only one point
  theta=atan(2*a*(x-mu))-pi/2
  
  x2=x+r*cos(theta)
  y2=a*(x-mu)^2+r*sin(theta)
  l=length(x2)
  #z=ou.sim.monte(deltat=0.25)
  
  #library(plotrix)

  ##Plot the parabola
  plot(x2,y2,type="l",lwd=2, yaxt="n", ylab="",asp=1, yaxt="n", frame.plot="F",
       xlim=c(-.5,1.05), ylim=c(0, 1), axes="F")
  axis(1, at=c(0, 0.5, 1.0))
  ##Draw dotted lines showing the edge of the field - methylation can't be bigger than
  ##1 or less than 0
  ##abline(v=0-r, lty=2,lwd=2,col="red")
  ##abline(v=1+r, lty=2,lwd=2,col="red")

  ##Draw a red circle
  draw.circle(m, a*(m-mu)^2, r, col="red")       
}

spring.plot <- function(m=0,mu=0.5) {  
  ##Draw a spring with a ball attached to the end of it - equilibrium of spring is 
  
  ##Spring attachment point is -(1.2-mu)
  e=-1+mu
  ##Spring width
  s=0.05
  ##Spring height 
  h=0.05

  ##Setup axes
  plot(c(e,e+s), c(0,0), type="l", lwd=2, xlim=c(-.5, 1.05), asp=1, axes="F", frame.plot="F")

  ##Make wall that spring is connected to
  lines(c(e, e), c(-1e3, 1e3), lwd=2)

  ##Make spring, it has 20 segments, 
  spring.x=seq(from=e+s, to=m-s, length=20)
  spring.y=c(0,h,-h,h,-h,h,-h,h,-h,h,-h,h,-h,h,-h,h,-h,h,-h,0)

  ##Draw spring
  lines(spring.x, spring.y, lwd=2)
  ##Draw line that connects to ball
  lines(c(m-s, m), c(0,0), lwd=2)
  ##Draw ball
  draw.circle(m,0,h, col='red')

  ##Draw dotted lines showing edge of field
  abline(v=0-h, lty=2,lwd=2,col="red")
  abline(v=1+h, lty=2,lwd=2,col="red")
  

}



broken.spring.plot <- function(m=0,mu=0.5) {  
  ##Draw a spring with a ball attached to the end of it - equilibrium of spring is 
  
  ##Spring attachment point is -(1.2-mu)
  e=-1+mu
  ##Spring width
  s=0.05
  ##Spring height 
  h=0.05

  ##Setup axes & plot wall-spring connector
  plot(c(e,e+s), c(0,0), type="l", lwd=2, xlim=c(-.5, 1.05), asp=1, axes="F", frame.plot="F")

  ##Make wall that spring is connected to
  lines(c(e, e), c(-1e3, 1e3), lwd=2)

  ##Make halves of broken spring
  wall.half.x=seq(from=(e+s), to=-.2, length=10)
  wall.half.y=c(0,h,-h,h,-h,h,-h,h,-h,h)

  ball.half.x=seq(from=m-s-.2, to=m-s, length=10)
  ball.half.y=c(-h,h,-h,h,-h,h,-h,h,-h,0)
  
  ##Draw spring
  lines(wall.half.x, wall.half.y, lwd=2)
  lines(ball.half.x, ball.half.y, lwd=2)
  ##Draw line that connects to ball
  lines(c(m-s, m), c(0,0), lwd=2)
  ##Draw ball
  draw.circle(m,0,h, col='red')

  ##Draw dotted lines showing edge of field
  abline(v=0-h, lty=2,lwd=2,col="red")
  abline(v=1+h, lty=2,lwd=2,col="red")
  

}



spring.para.plot <- function(x=0,mu=0.5) {
##To plot parabolic ball and ball and spring
  layout(c(1,2), heights=c(0.2,0.8))
  par(mar=c(1,0.2,1,0.2))
  
  ##Max value of potential  
  a=1/max(c(mu^2, (1-mu)^2))

  spring.plot(m=x,mu=mu)
  ball.plot(m=x,mu=mu, a=a)
}

breaking.spring.plot <- function(x, mu=0.5, to.a=0.1, num=20) {
  from.a=1/max(c(mu^2, (1-mu)^2))
  a.vals=seq(from=from.a, to=to.a, length=num)

  for (i in 1:length(a.vals)) {
    png(paste("Ball/e", formatC(i, digits=0, width=4, flag=0), ".png", sep=""), type="cairo")

    ##To plot parabolic ball and ball and spring
    layout(c(1,2), heights=c(0.2,0.8))
    par(mar=c(1,0.2,1,0.2))

    broken.spring.plot(m=x, mu=mu)
    ball.plot(m=x,mu=mu, a=a.vals[i])
    dev.off()
  }

}

broken.spring.para.plot <- function(x=0, mu=0.5, a=.1) {
##To plot parabolic ball and ball and spring
  layout(c(1,2), heights=c(0.2,0.8))
  par(mar=c(1,0.2,1,0.2))
  
  broken.spring.plot(m=x,mu=mu)
  ball.plot(m=x,mu=mu, a=a)

}

  
under.harmonic <- function (gamma=0.1) {
##Simulating underdamped harmonic oscillator

  t=seq(from=0, to=3e1, length=2e2)

  x=-0.5*exp(-gamma*t)*cos( (1-gamma^2)^(1/2)*t)+0.5

  for (i in 1:length(t)) {
    png(paste("Ball/u", i, ".png", sep=""), type="cairo")

    spring.para.plot(x[i])
    dev.off()
  }

}

over.harmonic <- function () {
##Simulating overdamped harmonic oscillator
  gamma=1.1
  c1=-(gamma - (gamma^2 -1)^(1/2))/(4*(gamma^2-1)^(1/2))
  c2=-(c1+.5)

  t=seq(from=0, to=3e1, length=2e2)

  x=exp(-gamma*t)*( c1 * exp(t) + c2 * exp(-t) ) +0.5

  for (i in 1:length(t)) {
    png( paste("Ball/o", i, ".png", sep=""))

    spring.para.plot1(x[i])
    dev.off()
  }

}

brown.harmonic <- function () {
  ##Simulating brownian harmonic oscillator
  z=ou.sim.monte(tfinal=3e1, deltat=3e1/2e2,m0=0,on.theta=1/1.1, off.theta=1/1.1)  

  for (i in 1:length(z$t)) {
    png( paste("Ball/b", formatC(i, digits=0, width=4, flag=0), ".png", sep=""),
        type="cairo")

    spring.para.plot(z$simul[i,1])
    dev.off()
  }
}

brown.harmonic.wall <- function () {
  ##Simulating brownian harmonic oscillator - using O-E
  z=ou.sim.monte(tfinal=3e1, deltat=3e1/2e2,m0=0,on.theta=1/1.1, off.theta=1/1.1, mu=1,draw=T)  

  for (i in 1:length(z$t)) {
    png( paste("Ball/b", i, ".png", sep=""), type="cairo")

    spring.para.plot(z$simul[i,1], mu=1)
    dev.off()
  }
}

cancer.harmonic <- function () {
  
  ##Simulating brownian harmonic oscillator - using O-E
  z=ou.sim.monte(tfinal=3e1, deltat=3e1/2e2,m0=0,on.theta=1/1.1, off.theta=0,
    mu=0.5,draw=T)  

  
  for (i in 1:200) {
    png( paste("Ball/cb", formatC(i, digits=0, width=4, flag=0), ".png", sep=""),
        type="cairo")    
    spring.para.plot(z$simul[i,1], mu=0.5)
    dev.off()
  }
  breaking.spring.plot(z$simul[200,1], num=20)

  for (i in 201:length(z$t)) {
    png( paste("Ball/ca", formatC((i-200), digits=0, width=4, flag=0, format='d')
               , ".png",sep=""), type="cairo")
    broken.spring.para.plot(x=z$simul[i,1])
    dev.off()
  }
}

cancer.noise <- function () {
  ##Simulating brownian harmonic oscillator with higher noise level  
  z=ou.sim.monte(tfinal=3e1, deltat=3e1/2e2, m0=0, on.theta=1/1.1, off.theta=1/1.1,
    mu=0.5, draw=T, off.sigma=.2)
  
  ##Plot pngs for movie
  for (i in 1:400) {
    ##Cairo allows us to access X11 backend w/o starting server
    png( paste("Ball/noi", formatC(i, digits=0, width=4, flag=0), ".png", sep=""),
        type="cairo")
    spring.para.plot(z$simul[i,1], mu=0.5)
    dev.off()                
  }

}
  
  
  




library(plotrix)
source("~/Code/timp_illumina/OU_model.R")
