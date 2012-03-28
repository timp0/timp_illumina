source("~/Code/timp_illumina/450k_cancer_loadin.R")

setwd("~/rgp")

rgp.data=preprocessIllumina(RGset[,RGset$Source=="RafaGP"])


rgp.beta=getBeta(rgp.data)

rgp.pData=pData(rgp.data)

write.csv(rgp.beta, file="rgp.beta.csv")

write.csv(rgp.pData, file="rgp.anno.csv")

