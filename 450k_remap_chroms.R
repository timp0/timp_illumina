
##This code is deisgned to test an entire chromosome, leaving out the broken parts

##This pulls in the command line stuff
args=(commandArgs(TRUE))
 
##Check if there are any variables, and if not, assign default
if(length(args)==0) {
  print("No arguments.")
  noise=1
}else{
  for(i in 1:length(args)) {
    print(args[[i]])
    eval(parse(text=args[[i]]))
  }
}

##Ok - let's try to align to the full human genome - why not, shouldn't take that long - and estimate the accuracy, etc.

##cnum=24
##fsize=5e4


##Finding chromosomal locations of sequences - part II(split by chromosome


#Load libraries in
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

##Get chromosome names
seqnames <- seqnames(Hsapiens)


library(GenomicRanges)


setwd("~/Data/Infinium/121311_analysis")

##Load file from Part 1
load("tmp.rda")


##Set chromosome number from read in value
i=cnum
##Get chromosome name
seqname <- seqnames[i]
##Load this chromosome
chromo <- Hsapiens[[seqname]]
##Make BS strand
chromo.v=as(chromo, "XStringViews")
##Or masks(chromo) <- NULL??

p.bs.chromo=chartr("C", "T", chromo.v)
##Rather than reverse complement(which screws up the numbering) - just switch
##G to A, and you have the reverseComplment of the bs minus strand
m.rc.bs.chromo=chartr("G", "A", chromo.v)

  
##For Illumina, the probes will either be the same as the
##reverse of the minus strand (m.rc.bs.chromo) or the reverse complement
##of the plus strand - alternatively, the reverse complment of the probe
##will be equal to the plus strand


cat(">>> Finding all hits in chromosome", seqname, "...\n")

in.chrom=which((as.character(seqnames(gprobes))==seqname) & !no.info)
##Source seq not empty, and chromosome matches one we are probing

chrom.probes=GRanges()
num.chrom=length(in.chrom)
ill=ill[in.chrom,]
numplus=numeric(num.chrom)
numminus=numeric(num.chrom)
full.strand=character(num.chrom)
full.r=IRanges(start=numeric(num.chrom), end=numeric(num.chrom))

for (j in 1:num.chrom) {   

  probe=DNAString(ill$pattern[j])
  rc.probe=reverseComplement(probe)
  cat(">> Finding probe", j, "of", num.chrom, "...\n")
  ##Search for matches
  plus.match=matchPattern(rc.probe,p.bs.chromo, fixed=F)
  minus.match=matchPattern(probe, m.rc.bs.chromo, fixed=F)
  ##Count num matches
  numplus[j]=length(plus.match)
  numminus[j]=length(minus.match)
  ##If matches, make a GRange of hits
  if (numplus[j]>0) {
    plus.g=GRanges(seqname=seqname, ranges=ranges(plus.match),
      strand="+", name=rep(ill$Name[j], numplus[j]),
      ori.index=rep(in.chrom[j], numplus[j]))
  } else {
    plus.g=GRanges()
  }

  if (numminus[j]>0) {
    minus.g=GRanges(seqname=seqname, ranges=ranges(minus.match),
      strand="-", name=rep(ill$Name[j], numminus[j]),
      ori.index=rep(in.chrom[j], numminus[j]))
  } else {
    minus.g=GRanges()
  }

  ##Concat all hits each loop - this is stupid and slow, but should work 
  chrom.probes=c(chrom.probes, plus.g, minus.g)  
}


cat("!!! Plus Matches = ", sum(numplus), ">1 plus match = ",
    sum(numplus>1), "...\n")

cat("!!! Minus Matches = ", sum(numminus), ">1 minus match = ",
    sum(numminus>1), "...\n")

mult.match=(numplus+numminus)>1

cat("!!! Plus and Minus Matches same probe = ", sum(mult.match), "...\n")

save(file=paste("chr_remap", cnum, ".rda", sep=""),
     list="chrom.probes")





