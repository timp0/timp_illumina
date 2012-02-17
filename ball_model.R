ball_plot1 <- function () {
  
  a=4
  mu=0.5
  r=0.05
  
  x=seq(from=0, to=1, by=0.001)
  y=a*(x-mu)^2
  theta=atan(2*a*(x-mu))-pi/2
  
  x2=x+r*cos(theta)
  y2=a*(x-mu)^2+r*sin(theta)
  
  z=ou_sim2(deltat=0.25)
  
  library(plotrix)
  
  for (i in 1:length(z$t)) {
    jpeg(paste("Ball/a", i, ".jpg", sep=""))
    
    plot(x2,y2,type="l",lwd=2)
    v=z$simul[i,1]
    
    draw.circle(v, a*(v-mu)^2, r, col="red")
    dev.off()
  }
       
}

spring_plot1 <- function(x=0) {  
  s=0.05
  h=0.05
  
  plot(c(-1,-1+s), c(0,0), type="l", lwd=2, xlim=c(-1, 1.05), asp=1)

  lines(c(-1, -1), c(-1e3, 1e3), lwd=2)

  spring_x=seq(from=-1+s, to=x-s, length=10)
  spring_y=c(0,h,-h,h,-h,h,-h,h,-h,0)
  
  lines(spring_x, spring_y, lwd=2)
  lines(c(x-s, x), c(0,0), lwd=2)
  draw.circle(x,0,h, col='red')

}



spring_para_plot1 <- function(x=0) {

  layout(c(1,2), heights=c(0.3,0.7))
  par(mar=c(1,2.2,1,0.2))



  
