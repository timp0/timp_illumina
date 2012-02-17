methyl_sim <- function(M_0=100, U_0=0, Em=.95, Ed=.05,nsimul=10,namey="simul1.png")  {


##Parameters
parms=c(Em=Em, Ed=Ed)

##Intial State(Methylated)
x0=c(M=M_0, U=U_0)

##Equations - methy maintain, methyl maintain fail, denovo methyl, denovo methyl fail
a=c("M*Em","M*(1-Em)", "U*Ed", "U*(1-Ed)")

##Direction that the states go in if the event(in columns) occurs
nu=matrix(c(0,-1,1,0, 0,1,-1,0),nrow=2, byrow=T)

final=50

t=0:final

runs=list(meth=numeric(),unmeth=numeric())

png(paste("Simul/", namey, sep=""))


plot(x=t, y=t, type="n", ylim=c(0,max(c(M_0, U_0))))

for (i in 1:nsimul) {
  z=ssa(x0,a,nu,parms,tf=final,simName="Methyl1")
  ##Get rid of redundant times
  redun=unique(z$data,margin=1)
  ##Get approximations
  m=approx(x=redun[,1], y=redun[,2], xout=t)
  u=approx(x=redun[,1], y=redun[,3], xout=t)
  runs$meth=cbind(runs$meth, m$y)
  runs$unmeth=cbind(runs$unmeth,u$y)

  ##Plot
  lines(z$data[,1], z$data[,2], type="l", lwd=1, col="green")
  lines(z$data[,1], z$data[,3], type="l", lwd=1, col="red")
}

legend("topright",c("Methylated", "Unmethylated"),col=c("green", "red"),lty=1,lwd=1)

dev.off()
return(runs)
}

sense_test <- function() {
  
  senseEm=seq(from=0, to=1, by=0.05)
  senseEd=seq(from=0, to=1, by=0.05)
  
  for (i in senseEm) {
    for (j in senseEd) {
      
      methyl_sim(Em=i, Ed=j, namey=paste(i, j, "test.png", sep="_"))
    }
  }
}


library(GillespieSSA)
