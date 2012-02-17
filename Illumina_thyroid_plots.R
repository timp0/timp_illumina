##This code is used to make thyroid plots I want

##dev.off()
library(plotrix)

thy=data$fsamp$Class==7

thy_data=list(fqbeta=data$fqbeta[,thy], fsamp=data$fsamp[thy,])

p=prcomp(t(thy_data$fqbeta))

norms=thy_data$fsamp$Progression==0
tums=thy_data$fsamp$Progression>0





##First Thyroid type plot - set colors - open circles are normal, closed circles are tumor
type_col=factor(abs(thy_data$fsamp$Other_Note))

levels(type_col)=c("red", "green", "blue", "black", "darkorchid", "orange", "aquamarine")

pdf("Movie/thy_types_pca.pdf")

plot(p$x[tums,1], p$x[tums,2], col=type_col[tums], pch=16,xlim=range(p$x[,1]), ylim=range(p$x[,2]))
points(p$x[norms,1], p$x[norms,2], col=type_col[norms], pch=1)

legend("topright", c("ADN", "FA", "FC", "FVPTC", "HA", "HC", "PTC"),col=levels(type_col),pch=16)

dev.off()


##Now thyroid progression plot
prog_col=factor(thy_data$fsamp$Progression)

levels(prog_col)=c("green", "orange", "blue", "red")

pnorms=p$x[norms,]
pn_med=apply(pnorms,2,median)
pn_mad=apply(pnorms,2,mad)

pdf("Movie/thy_prog_pca.pdf")

plot(p$x[,1], p$x[,2], bg=as.character(prog_col), pch=21)
draw.ellipse(x=rep(pn_med[1],2), y=rep(pn_med[2],2), a=rep(pn_mad[1]*5,2), b=rep(pn_mad[2]*5,2), lty=2, border="red" )
legend("topright", c("Normal", "Nodule", "Adenoma", "Carcinoma"), col=levels(prog_col), pch=16)

dev.off()

##Thyroid Maha Dist from Normal

thyn_med=apply(thy_data$fqbeta[,norms],1,median)
thyn_mad=apply(thy_data$fqbeta[,norms],1,mad)
thyt_mad=apply(thy_data$fqbeta[,tums],1,mad)

comp_mad=thyt_mad-thyn_mad

thya_mad=apply(thy_data$fqbeta, 1, mad)

top25=order(comp_mad, decreasing="T")[1:25]

anti_norm=sweep( abs(sweep(thy_data$fqbeta,1,thyn_med)), 1, thyn_mad, FUN="/")

anti_norm_score=apply(anti_norm[top25,],2,sum)


prog_x=factor(thy_data$fsamp$Progression)

levels(prog_x)=c(1,2,3,4)

pdf("Movie/thy_prog_maha_scat.pdf")

plot(jitter(as.numeric(prog_x)), anti_norm_score, bg=as.character(prog_col), pch=21, xaxt="n", xlab="")
legend("topleft", c("Normal", "Nodule", "Adenoma", "Carcinoma"), col=levels(prog_col), pch=16)

dev.off()

pdf("Movie/thy_type_maha_scat.pdf")

plot(jitter(as.numeric(prog_x)), anti_norm_score, bg=as.character(type_col), pch=21, xaxt="n", xlab="")
legend("topleft", c("ADN", "FA", "FC", "FVPTC", "HA", "HC", "PTC"), col=levels(type_col), pch=16)

dev.off()

source("~/Work/Analysis/Illumina/Illumina_cluster1b.R")

##simplecluster <- function(avg, sampy, which_tissue=3, which_type=NA, which_grade=NA, which_progression=NA, divider=NA, namey="a.pdf")  {
simplecluster(thy_data$fqbeta, thy_data$fsamp, which_tissue=7, divider=thy_data$fsamp$Progression, namey="thy_cluster_prog.pdf")
simplecluster(thy_data$fqbeta, thy_data$fsamp, which_tissue=7, divider=abs(thy_data$fsamp$Other_Note), namey="thy_cluster_type.pdf")

#Zeiger requested - tumor only
simplecluster(thy_data$fqbeta[,tums], thy_data$fsamp[tums,], which_tissue=7, divider=thy_data$fsamp$Progression[tums], namey="thy_cluster_tum_prog.pdf")

##Now thyroid progression plot
prog_col=factor(thy_data$fsamp$Progression[tums])

levels(prog_col)=c("orange", "blue", "red")

#pnorms=p$x[norms,]
#pn_med=apply(pnorms,2,median)
#pn_mad=apply(pnorms,2,mad)
tonly=prcomp(t(thy_data$fqbeta[,tums]))

pdf("Movie/thy_prog_tums_pca.pdf")

plot(tonly$x[,1], tonly$x[,2], bg=as.character(prog_col), pch=21)
#draw.ellipse(x=rep(pn_med[1],2), y=rep(pn_med[2],2), a=rep(pn_mad[1]*5,2), b=rep(pn_mad[2]*5,2), lty=2, border="red" )
legend("topright", c("Nodule", "Adenoma", "Carcinoma"), col=levels(prog_col), pch=16)

dev.off()

#Figure out best way to plot thyroid probes - find ones which distinguish adn from carc best?
