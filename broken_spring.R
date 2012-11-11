##Get path for code on this machine
if (!exists("codedir")) {
  codedir=getwd()
}


##Put us in the right dir for output
plotdir="~/Data/Genetics/Infinium/111012_model/"


source(paste(codedir, "ball_model.R", sep="/"))

brown.harmonic(plotdir=plotdir)

pdf(file.path(plotdir, "spring_para.pdf"))
broken.spring.para.plot(x=.5)
spring.para.plot(x=.5)
dev.off()


##cancer.harmonic()

##cancer.noise()
