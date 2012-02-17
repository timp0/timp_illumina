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
rm(normG)
rm(normR)

## remove bad arrays
library(matrixStats)
cc <- colMads(qlogis(data$qbeta))
pdf("Movie/qlog_qc1.pdf")
hist(cc, nc=100)
abline(v=1.7, lty=1)
dev.off()

