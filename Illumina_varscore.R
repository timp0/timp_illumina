##Based on Hector's runit.R

##srcdir <- "~hcorrada/barcode/src/succs"
source(file.path(srcdir,"succs2.R"))

                                        
beta <- beta[,ipd$tissue != "Ovary"]
ipd <- ipd[ipd$tissue != "Ovary",]

y <- ifelse(ipd$type=="Normal", 0, 1)
succsig <- succit(y, beta, sigsize=25, mult=5, filtersize=100, unidir=TRUE)

pred <- predict.succs(succsig, beta)
cnts <- attr(pred, "counts")
o <- order(cnts)

pdf("score.pdf")
mypar(1,1)
plot(cnts[o], col=y[o]+1, type="h", main="Score",
     xlab="Sample", ylab="vMeth Score")
dev.off()
     
require(matrixStats)
MED <- rowMedians(beta[,ipd$type=="Normal"])
ann$mstat <- "Meth"
ann$mstat[MED<0.5] <- "Unmeth"

succ.ann <- ann[succsig$filtered,]

xx <- succify(beta[succsig$filtered,], succsig$stats)
succ.ann$SON <- rowSums(xx)
tmp <- succ.ann[succsig$indices,]
succ.ann <- rbind(tmp, succ.ann[-succsig$indices,])

tmp <- ann[-succsig$filtered,]
tmp$SON <- NA
succ.ann <- rbind(succ.ann, tmp)

save(succsig, succ.ann, file="succsig_meth.rda")

loo.res <- vector("list", length=length(y))
for (i in 1:length(y)) {
  train.e <- beta[,-i]
  train.y <- y[-i]

  succ.res.mult <- succit(train.y, train.e, sigsize=20, filtersize=100, mult=5)
  pred.mult <- predict.succs(succ.res.mult, beta[,i,drop=FALSE])

  succ.res.auto <- succit(train.y, train.e, sigsize=20, filtersize=100, mult=-1)
  pred.auto <- predict.succs(succ.res.auto, beta[,i,drop=FALSE])
  loo.res[[i]] <- list(sig.mult=succ.res.mult, sug.auto=succ.res.auto,
                       pred.mult=pred.mult, pred.auto=pred.auto)
}

tissues <- unique(ipd$tissue)
loco.res <- vector("list", length=length(tissues)) 
for (i in 1:length(tissues)) {
  tt <- tissues[i]
  test.inds <- which(ipd$tissue == tt)

  train.e <- beta[,-test.inds]
  train.y <- y[-test.inds]

  test.e <- beta[,test.inds]
  test.y <- y[test.inds]
  
  succ.res.mult <- succit(train.y, train.e, sigsize=40, filtersize=120, mult=5)
  pred.mult <- predict.succs(succ.res.mult, test.e)

  succ.res.auto <- succit(train.y, train.e, sigsize=40, filtersize=120, mult=-1)
  pred.auto <- predict.succs(succ.res.auto, test.e)
  
  loco.res[[i]] <- list(sig.mult=succ.res.mult, sig.auto=succ.res.auto,
                        pred.mult=pred.mult, y=test.y, pred.auto=pred.auto)
}
names(loco.res) <- tissues
save(loo.res, loco.res, file="succmeth.rda")

