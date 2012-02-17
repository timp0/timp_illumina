
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

##Align bisulfite probes to bisulfite genome

##We can't trust that the probes don't align multiple locations.  So.  I am taking
##each probe against entire genome.  Divided probes into bunches of 25000

##pgsize=25000
pgsize=250
##pgn is group number - ranges from 0 to 19

##pgn=19



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

##Init range of probes
prange=numeric(2)
prange[1]=1+pgn*pgsize
prange[2]=min(c(pgn+1)*pgsize, dim(ill)[1])
pindex=prange[1]:prange[2]

##Take just those probes
ill=ill[pindex,]



probe.map=GRanges()
num.probe=dim(ill)[1]
numplus=numeric(num.probe)
numminus=numeric(num.probe)

probe.set=DNAStringSet(ill$pattern)
rc.probe.set=reverseComplement(probe.set)

##Cycle chromosomes - this is computational intensive(BS conv), so do minimum times
##It only takes a second for X chromosome(example), but still.
for (i in 1:24) {
  ##Get chromosome name
  seqname=seqnames[i]
  ##Load this chromosome
  chromo=Hsapiens[[seqname]]
  ##Strip masks
  masks(chromo) <- NULL
  ##OR chromo.v=as(chromo, "XStringViews")

  ##Make BS strand
  p.bs.chromo=chartr("C", "T", chromo)
  ##Rather than reverse complement(which screws up the numbering) - just switch
  ##G to A, and you have the reverseComplment of the bs minus strand
  m.rc.bs.chromo=chartr("G", "A", chromo)
 
  ##For Illumina, the probes will either be the same as the
  ##reverse of the minus strand (m.rc.bs.chromo) or the reverse complement
  ##of the plus strand - alternatively, the reverse complment of the probe
  ##will be equal to the plus strand
  




  for (j in 1:num.probe) {   
    cat(">> Finding probe", j, "of", num.probe, " in chromosome", seqname, "...\n")
    ##Search for matches
    plus.match=matchPattern(rc.probe.set[[j]],p.bs.chromo, fixed=c(pattern=F, subject=T))
    minus.match=matchPattern(probe.set[[j]], m.rc.bs.chromo, fixed=c(pattern=F, subject=T))
    ##Count num matches
    numplus[j]=length(plus.match)
    numminus[j]=length(minus.match)

    ##If matches, make a GRange of hits
    if (numplus[j]>0) {
      plus.g=GRanges(seqname=seqname, ranges=ranges(plus.match),
        strand="+", name=rep(ill$Name[j], numplus[j]),
        ori.index=rep(pindex[j], numplus[j]))
    } else {
      plus.g=GRanges()
    }
    
    if (numminus[j]>0) {
      minus.g=GRanges(seqname=seqname, ranges=ranges(minus.match),
        strand="-", name=rep(ill$Name[j], numminus[j]),
        ori.index=rep(pindex[j], numminus[j]))
    } else {
      minus.g=GRanges()
    }
    
    ##Concat all hits each loop - this is stupid and slow, but should work 
    probe.map=c(probe.map, plus.g, minus.g)  
  }
}


cat("!!! Plus Matches = ", sum(numplus), ">1 plus match = ",
    sum(numplus>1), "...\n")

cat("!!! Minus Matches = ", sum(numminus), ">1 minus match = ",
    sum(numminus>1), "...\n")

mult.match=(numplus+numminus)>1

cat("!!! Plus and Minus Matches same probe = ", sum(mult.match), "...\n")

save(file=paste("remap", pga, ".rda", sep=""),
     list="probe.map")





