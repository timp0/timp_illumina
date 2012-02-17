
##This code is designed to align the probes - using matchPdict is superior, so I've
##Switched to that from matchPattern, removing ambiguous bases

##Align bisulfite probes to bisulfite genome

##We can't trust that the probes don't align multiple locations.  So.  I am taking
##each probe against entire genome. 


##Finding chromosomal locations of sequences - part II

st.time=proc.time()

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

probe.set=DNAStringSet(ill$pattern)

##Bases in probe.set
probe.bases=alphabetFrequency(probe.set)

##Probes with odd bases - G(which shouldn't exist post-bisulfite,
##Y(Which is from the SNP probes, deal with these seperately
strange.probes=(probe.bases[,3]!=0)|(probe.bases[,9]!=0)

##Subset strange probes
stg.ill=ill[strange.probes,]

##Take just those probes which aren't strange
ill=ill[!strange.probes,]
probe.set=probe.set[!strange.probes]


remap.probe=GRanges()
num.probe=dim(ill)[1]
numplus=numeric(num.probe)
numminus=numeric(num.probe)

##Since we are aligning to completely unmethylated genome(all Cs to Ts),
##We can assume R(G or A) is only A
probe.set=chartr("R", "A", probe.set)
rc.probe.set=reverseComplement(probe.set)
probe.map=PDict(probe.set)
rc.probe.map=PDict(rc.probe.set)


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
  

  cat(">> Finding probes"," in chromosome", seqname, "...\n")
  ##Search for matches
  plus.match=matchPDict(rc.probe.map,p.bs.chromo)
  minus.match=matchPDict(probe.map, m.rc.bs.chromo)

  ##Find hits against this chromosome
  numplus=countIndex(plus.match)
  numminus=countIndex(minus.match)

  ##Get name/indexes

  
  ##Init GRanges
  plus.g=GRanges(seqname=seqname, ranges=unlist(plus.match),
    strand="+", name=rep(ill$Name, numplus), ori.index=rep(ill$pidx, numplus))
  minus.g=GRanges(seqname=seqname, ranges=unlist(minus.match),
    strand="-", name=rep(ill$Name, numminus), ori.index=rep(ill$pidx, numminus))
  
  ##Concat all hits each loop - this is stupid and slow, but should work 
  remap.probe=c(remap.probe, plus.g, minus.g)  
}

fi.time=proc.time()
fi.time-st.time
##25 min to run


save(file="remap.rda", list="remap.probe")

##Map strange probes
probe.stg=DNAStringSet(stg.ill$pattern)
##Move Gs(don't understand why they are there - to Ns
probe.stg=chartr("G", "N", probe.stg)
rc.probe.stg=reverseComplement(probe.stg)

remap.stg=GRanges()
num.probe=length(probe.stg)
numplus=numeric(num.probe)
numminus=numeric(num.probe)

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
    plus.match=matchPattern(rc.probe.stg[[j]],p.bs.chromo, fixed=c(pattern=F, subject=T))
    minus.match=matchPattern(probe.stg[[j]], m.rc.bs.chromo, fixed=c(pattern=F, subject=T))
    ##Count num matches
    numplus[j]=length(plus.match)
    numminus[j]=length(minus.match)

    ##If matches, make a GRange of hits
    if (numplus[j]>0) {
      plus.g=GRanges(seqname=seqname, ranges=ranges(plus.match),
        strand="+", name=rep(stg.ill$Name[j], numplus[j]),
        ori.index=rep(stg.ill$pidx[j], numplus[j]))
    } else {
      plus.g=GRanges()
    }
    
    if (numminus[j]>0) {
      minus.g=GRanges(seqname=seqname, ranges=ranges(minus.match),
        strand="-", name=rep(stg.ill$Name[j], numminus[j]),
        ori.index=rep(stg.ill$pidx[j], numminus[j]))
    } else {
      minus.g=GRanges()
    }
    
    ##Concat all hits each loop - this is stupid and slow, but should work 
    remap.stg=c(remap.stg, plus.g, minus.g)  
  }
}


cat("!!! Plus Matches = ", sum(numplus), ">1 plus match = ",
    sum(numplus>1), "...\n")

cat("!!! Minus Matches = ", sum(numminus), ">1 minus match = ",
    sum(numminus>1), "...\n")

mult.match=(numplus+numminus)>1

cat("!!! Plus and Minus Matches same probe = ", sum(mult.match), "...\n")

remap.probe=c(remap.probe, remap.stg)

save(file="remap_f.rda", list=c("remap.probe","strange.probes") )

##Optional - save out the multiple matches
##Load back ill input file
##load("tmp.rda")

##two.match=ill$Name[as.numeric(names(y[y==2]))]
##three.match=ill$Name[as.numeric(names(y[y==3]))]
##six.match=ill$Name[as.numeric(names(y[y==6]))]
##seventeen.match=ill$Name[as.numeric(names(y[y==17]))]

##num.mult=c(rep(2,length(two.match)), rep(3, length(three.match)),
##  rep(6, length(six.match)), rep(17, length(seventeen.match)))

##multmatch=data.frame(num.match=num.mult, name=character(length(num.mult)),
##  stringsAsFactors=F)
##multmatch$name=c(two.match, three.match, six.match, seventeen.match)

##save(file="multi.rda", list="multmatch")


