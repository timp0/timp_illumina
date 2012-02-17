
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

##cnum=22
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
p.bs.chromo=chartr("C", "T", chromo.v)
##Rather than reverse complement(which screws up the numbering) - just switch
##G to A, and you have the reverseComplment of the bs minus strand
m.rc.bs.chromo=chartr("G", "A", chromo.v)

  
##For Illumina, the probes will either be the same as the
##reverse of the minus strand (m.rc.bs.chromo) or the reverse complement
##of the plus strand - alternatively, the reverse complment of the probe
##will be equal to the plus strand


cat(">>> Finding all hits in chromosome", seqname, "...\n")

in.chrom=(as.character(seqnames(gprobes))==seqname) & !no.info
##Source seq not empty, and chromosome matches one we are probing

chrom.probes=gprobes[in.chrom]
num.chrom=length(chrom.probes)
ill=ill[in.chrom,]
numplus=numeric(num.chrom)
numminus=numeric(num.chrom)
full.strand=character(num.chrom)
full.start=numeric(num.chrom)
full.end=numeric(num.chrom)

for (j in 1:num.chrom) {   

  probe=DNAString(ill$pattern[j])
  rc.probe=reverseComplement(probe)
  cat(">> Finding probe", j, "of", num.chrom, "...\n")
  ##Plus first
  plus.match=matchPattern(rc.probe,p.bs.chromo, fixed=F)
  minus.match=matchPattern(probe, m.rc.bs.chromo, fixed=F)
  
  numplus[j]=length(plus.match)
  numminus[j]=length(minus.match)
  
  if (numplus[j] > 0) {
    full.strand[j]="+"
    full.start=start(plus.match[1])
    full.end=end(plus.match[1])
  } else {
    if (numminus[j] > 0) {
      full.strand[j]="-"
      full.start=start(minus.match[1])
      full.end=end(minus.match[1])
    }
  }
  
}


cat("!!! Plus Matches = ", sum(numplus), " >1 plus match = ",
    sum(numplus>1), "...\n")

cat("!!! Minus Matches = ", sum(numminus), " >1 minus match = ",
    sum(numminus>1), "...\n")

start(chrom.probes)=full.start
end(chrom.probes)=full.end
strand(chrom.probes)=full.strand

        
save(file=paste("chr_remap", cnum, ".rda", sep=""),
     list="chrom.probes")





