##Ok - make some plots specific for webpage


setwd("~/Data/Illumina_Bead/Analysis/")

load("predata.rda")

library(RColorBrewer)


##Find norms and cancers
norms=data$fsamp$Progression==0
carcy=data$fsamp$Progression>2

##Far is 2, shore is 1, island is 0
islstatus=(data$probes$UCSC_Dist_to_Island>0)+(data$probes$UCSC_Dist_to_Island>2000)


cgi.colors <- brewer.pal(8,"Dark2")[c(4:5,8)]

normvar=apply(data$fqbeta[,norms],1,mad)
cancvar=apply(data$fqbeta[,carcy],1,mad)

pdf("Movie/allvar1.pdf")

rangy=max(max(normvar), max(cancvar))

plot(normvar,cancvar, xlab="Normal", ylab="Cancer",
     xlim=c(0, rangy), ylim=c(0,rangy),
       bg=cgi.colors[islstatus+1], pch=21)


##Signifcance lines
cc <- qf(.99, sum(carcy)-1, sum(norms)-1)
abline(0,sqrt(cc), lty=2)
abline(0,1)

legend("bottomright", c("Island", "Shore", "Far"), pch=21, pt.bg=cgi.colors)

dev.off()
