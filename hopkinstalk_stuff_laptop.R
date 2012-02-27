##Get path for code on this machine
if (!exists("codedir")) {
  codedir=getwd()
}


##Put us in the right dir for output
setwd("~/Data/Infinium/011412_analysis/")


source(paste(codedir, "ball_model.R", sep="/"))

brown_harmonic1()

cancer.harmonic()

cancer.noise()
