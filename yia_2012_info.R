##This code is to get the info for the YIA 2012 essay

#Load libraries in
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

##Get chromosome names
seqnames <- seqnames(Hsapiens)


library(GenomicRanges)


##refseq.genes=updateObject(refseq.genes)
##ucsc.isl=updateObject(ucsc.isl)


##Get DMR info
load("~/Data/Genetics/092011_Capture/dmr_cap.rda")

##Get Block/LOCK/LAD info
load("~/Data/Genetics/072111_blocks/lg_regions2.rda")


hyper.blocks=values(ccancer.blocks)$delta.M>0
hypo.blocks=values(ccancer.blocks)$delta.M<0

##average hypo in blocks
mean(values(ccancer.blocks)$delta.M[hypo.blocks])
mean(rep(values(ccancer.blocks)$delta.M[hypo.blocks],
         values(ccancer.blocks)$numCpGs[hypo.blocks]))

mean(values(ccancer.blocks)$numCpGs[hypo.blocks])

mean(values(ccancer.blocks)$delta.M[hyper.blocks])
mean(rep(values(ccancer.blocks)$delta.M[hyper.blocks],
         values(ccancer.blocks)$numCpGs[hyper.blocks]))

mean(values(ccancer.blocks)$numCpGs[hyper.blocks])
