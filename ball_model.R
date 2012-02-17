a=4
mu=0.5

x=seq(from=0, to=1, by=0.001)
y=a*(x-mu)^2
theta=atan(2*a*(x-mu))-pi/2

x2=x+0.05*cos(theta)
y2=a*(x-mu)^2+0.05*sin(theta)

library(plotrix)

plot(x2,y2,type="l")
lines(x,y,col="red")

draw.circle(x[300], y[300], 0.05)

b=1
z=b*(x-1)^2
tz=atan(2*b*(x-1))-pi/2

x3=x+0.05*cos(tz)
y3=b*(x-1)^2+0.05*sin(tz)

plot(x3,y3,type="l",xlim=c(-0.05,1.05),ylim=c(-0.05,1.05))
lines(c(1,1.05,1.05),c(-0.05,-0.05,1e5))
lines(x,z,col="red")



