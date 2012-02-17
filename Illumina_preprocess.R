##Ok - taken from preprocessing.R from Hector Bravo

###Lets subset out the controls
#keepIndex=!dat$ipd$control
#R=dat$R[,keepIndex]
#G=dat$G[,keepIndex]
#dat$ipd2 <- dat$ipd[keepIndex,]

###normalize data
normG=preprocessCore::normalize.quantiles(data$green)
normR=preprocessCore::normalize.quantiles(data$red)
data$qbeta=normG/(normR+normG)

## remove bad arrays
library(matrixStats)
cc <- colMads(qlogis(data$qbeta))
mG <- colMedians(log(normG))
mR <- colMedians(log(normR))
a <- 0.5*(mG+mR)
keep <- cc>=1.9 & a>7

## Get rid of bad samples in either way
data$samp$keep=keep&(data$samp$Exclude==0)

data$fqbeta=data$qbeta[,data$samp$keep]
data$fsamp=data$samp[data$samp$keep,]


##pdf("Movie/qlog_qc1.pdf")
#hist(cc, nc=100)
#abline(v=1.7, lty=1)
#dev.off()

save(data, data_norm, file="predata.rda")

##


