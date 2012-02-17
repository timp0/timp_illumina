##This code takes aligned probes, find the sbe, and figures out snp status
##As well as blocks, etc.


#Load libraries in
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

##Get chromosome names
seqnames <- seqnames(Hsapiens)


library(GenomicRanges)


setwd("~/Data/Infinium/121311_analysis")


load(file="probe_obj_final.rda")

##Need to calc dist to CpG Island, dist to genes
##Load islands and genes object
load("~/Data/Genetics/072111_blocks/gene_island.rda")

##Find nearest UCSC CpG islands to each probe (including overlap)
close=nearest(gprobes, hmm.isl)

##UCSC island index, and dist to ucsc isl
values(gprobes)$nearest.hmmisland.index=close
values(gprobes)$dist.hmmisland=width(pgap(ranges(gprobes),
                 ranges(hmm.isl[close])))

remap.index=order(match(as.character(seqnames(remap.probe)), seqnames),
  start(remap.probe))

gprobes.index=order(match(as.character(seqnames(gprobes)), seqnames),
  start(gprobes))

values(gprobes)$map.idx=remap.index[values(gprobes)$map.idx]

values(remap.probe)$ori.index=gprobes.index[values(remap.probe)$ori.index]
values(sbe)$ori.index=gprobes.index[values(sbe)$ori.index]

gprobes=gprobes[gprobes.index]
sbe=sbe[remap.index]
remap.probe=remap.probe[remap.index]

save(file="probe_obj_rafa.rda", list=c("gprobes", "sbe", "remap.probe"))

