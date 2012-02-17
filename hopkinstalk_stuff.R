##Put us in the right dir for output
setwd("~/Data/Infinium/011412_analysis/")

##Experimental file location
expdatapath="/thumper2/feinbergLab/core/arrays/illumina/"

##Current 450k package name
##Kasper says ignore warnings about "contains no R code"
##This is the local(expanded) version of minfi
library(minfi)
library(minfiLocal)

##This will list commands
##ls("package:minfi")


##Read in plates
plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))

##Read in data
RGset=read.450k.exp(base=expdatapath, targets=plates)

load("~/Data/Infinium/121311_analysis/probe_obj_final.rda")

##Ok - regions
##DMR0 ncbi36:11, 2125904-2126160 - Ito et al
##DMR0 hg19 chr11:2169328-2169584
dmr0=GRanges(seqnames="chr11", ranges=IRanges(start=2169328,
                                 end=2169584), id="dmr0")
##ICR is -2kb to -4kb of h19
h19.start=2019065
icr=GRanges(seqnames="chr11", ranges=IRanges(start=2019065+2000,
                                end=2019065+4000), id="icr")


load("~/Data/Genetics/072111_blocks/gene_island.rda")
##Get IGF2 gene location
igf2.gene=reduce(refseq.genes[values(refseq.genes)$gene.name=="IGF2"])

igf2.reg=igf2.gene
values(igf2.reg)=NULL
values(igf2.reg)$id="IGF2"

loi.reg=c(dmr0, icr, igf2.reg)




##Get out just colon samples and mets
tis.samp=(pData(RGset)$Tissue %in% c("colon", "lung", "breast", "thyroid",
  "kidney", "pancreas")) & (pData(RGset)$Phenotype %in% c("normal", "cancer"))
##Alternatively pData(RGset) gets out the plates variable again

tis.data=RGset[,(tis.samp)]

tis.MSet=preprocessMinfi(tis.data)



##Add something indexing gprobes vs actual data order
values(gprobes)$minfi.idx=match(values(gprobes)$name,
                 rownames(getM(tis.MSet[,1])))

##Obtain those probes
icr.probes=values(gprobes)$minfi.idx[(gprobes %in% loi.reg[2])]


##Keep just probes with no SNPs
good.probes=values(gprobes)$minfi.idx[ (!values(gprobes)$sbe.snp.boo)&
  (!values(gprobes)$boo.snps)]

## ok - for dot plot, need to set y values to the different tissue types, and add jitter 
tissue.y=as.numeric(factor(pData(tis.data)$Tissue))*2-
  (as.numeric(factor(pData(tis.data)$Status))-1)

tis.y.lev=levels(factor(pData(tis.data)$Tissue))


##Coloring for dots tumor normal
colly=c("red", "green")

state.col=factor(pData(tis.data)$Status)

levels(state.col)=colly

ybig=max(tissue.y)+1

tissue.y=jitter(tissue.y)

tis.beta=getBeta(tis.data)

pdf("Plots/icr1.pdf")


##Get methyl and unmethyl signal
for (i in icr.probes) {

  plot(tis.beta[i,], tissue.y,
       bg=as.character(state.col),
       pch=21, xlim=c(0,1), ylim=c(0,ybig),
       yaxt="n", ylab="")
  axis(2, at=as.numeric(factor(tis.y.lev))*2-0.5,
       labels=tis.y.lev)
  
}
dev.off()

##Get Block/LOCK/LAD info
load("~/Data/Genetics/072111_blocks/lg_regions2.rda")

smdmr.type=factor(values(ccancer.seq.dmrs)$type)
table(smdmr.type)

tapply(values(ccancer.seq.dmrs)$numCpGs, smdmr.type, median)
tapply(width(ccancer.seq.dmrs), smdmr.type, median)

hyper.blocks=values(ccancer.blocks)$delta.M>0
hypo.blocks=values(ccancer.blocks)$delta.M<0

sum(hyper.blocks)
range(values(ccancer.blocks)$numCpGs[hyper.blocks])
median(values(ccancer.blocks)$numCpGs[hyper.blocks])
sum(values(ccancer.blocks)$numCpGs[hyper.blocks])
range(width(ccancer.blocks[hyper.blocks]))
median(width(ccancer.blocks[hyper.blocks]))
sum(width(ccancer.blocks[hyper.blocks]))


sum(hypo.blocks)
range(values(ccancer.blocks)$numCpGs[hypo.blocks])
median(values(ccancer.blocks)$numCpGs[hypo.blocks])
sum(values(ccancer.blocks)$numCpGs[hypo.blocks])
range(width(ccancer.blocks[hypo.blocks]))
median(width(ccancer.blocks[hypo.blocks]))
sum(width(ccancer.blocks[hypo.blocks]))

##All blocks
length(ccancer.blocks)
range(values(ccancer.blocks)$numCpGs)
median(values(ccancer.blocks)$numCpGs)
sum(values(ccancer.blocks)$numCpGs)
range(width(ccancer.blocks))
median(width(ccancer.blocks))
sum(width(ccancer.blocks))


##Find gene/block overlap
z=countOverlaps(refseq.genes, ccancer.blocks)
sum(z==0)
sum(z!=0)
length(z)

y=countOverlaps(ccancer.blocks, refseq.genes)
sum(y==0)
sum(y!=0)



load("~/Data/Genetics/092011_Capture/dmr_cap.rda")

z=countOverlaps(tc.hg19, ccancer.blocks)
sum(z==0)
sum(z!=0)
length(z)

y=countOverlaps(tc.hg19, ccancer.seq.dmrs)
sum(y==0)
sum(y!=0)



x=countOverlaps(tc.hg19, tdmrs.hg19)
sum(x==0)
sum(x!=0)

w=countOverlaps(tc.hg19, cdmrs.hg19)
sum(w==0)
sum(w!=0)

v=countOverlaps(tc.hg19, rdmrs.hg19)
sum(v==0)
sum(v!=0)


source("~/Work/Analysis/Illumina/ball_model4.R")

under_harmonic1()
