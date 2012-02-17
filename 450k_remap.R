#Finding chromosomal locations of sequences


#Load libraries in
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

##Get chromosome names
seqnames <- seqnames(Hsapiens)


library(GenomicRanges)


setwd("~/Data/Infinium/121311_analysis")

#Load the Probe Manifest File
ill <- read.csv("HumanMethylation450_15017482_v.1.1.csv",stringsAsFactors=FALSE, skip=7)

cntrl.probes=ill[485578:486428,]
ill=ill[1:485577,]
 

##Change strand to "+" or "-"
ill$Strand[ill$Strand=="F"]="+"
ill$Strand[ill$Strand=="R"]="-"
ill$Strand[ill$Strand==""]="*"

##Annotate the 65 rs probes with no info
no.info=ill$CHR==""
## I went online and queried these rs##
dbsnp.query=read.csv(file="dbsnp_query1.csv", stringsAsFactors=F)
## Keep just the hg19 
dbsnp.query=dbsnp.query[dbsnp.query$a=="GRCh37.p5",]
##Make numbers into rsnumber for names to match
dbsnp.query$names=paste("rs", dbsnp.query$rs., sep="")
##SNP rs13369115 was merged into rs10155413
dbsnp.query$names[57]="rs13369115"
##Match to locs in probe table
q=match(dbsnp.query$names, ill$Name)

##ill$Genome_Build[q]=37
ill$CHR[q]=dbsnp.query$chr
ill$MAPINFO[q]=dbsnp.query$chrpos
##ill$avHet[q]=dbsnp.query$avghet
##ill$avHetSE[q]=dbsnp.query$s.e.het


##Make all the patterns for matching
##Probe II
ptype=ill$Infinium_Design_Type=="II"
ill$pattern[ptype]=ill$AlleleA_ProbeSeq[ptype]

##Probe I - much more complicated
ptype=which(ill$Infinium_Design_Type=="I")

probe.a=DNAStringSet(ill$AlleleA_ProbeSeq[ptype])
probe.b=DNAStringSet(ill$AlleleB_ProbeSeq[ptype])
probe.com=compareStrings(probe.a, probe.b)

for (i in 1:length(ptype)) {
  probe.mis=gregexpr("\\?", probe.com[i])
  for (j in as.numeric(probe.mis[[1]])) {
    substr(probe.com[i],start=j, stop=j)=
      mergeIUPACLetters(as.character(c(probe.a[[i]][j], probe.b[[i]][j])))
  }
}

ill$pattern[ptype]=probe.com

#Initialize hit table

##Make GRanges of probes
gprobes=GRanges(seqname=paste("chr", ill$CHR, sep=""),
  strand=ill$Strand,ranges=IRanges(start=ill$MAPINFO, end=ill$MAPINFO))

values(gprobes)$name=ill$Name


##Booleans from probe design origins
values(gprobes)$random.loci=!is.na(ill$Random_Loci)
values(gprobes)$M27k.array=!is.na(ill$Methyl27_Loci)
values(gprobes)$fantom=!(ill$Phantom=="")
values(gprobes)$fantom.nfo=ill$Phantom
values(gprobes)$enhancer=!is.na(ill$Enhancer)
values(gprobes)$dnase.hypersensitivity=!is.na(ill$DHS)

##Type of probe info
values(gprobes)$probe.type=ill$Infinium_Design_Type
values(gprobes)$probe.nxt.base=ill$Next_Base
values(gprobes)$probe.color.chan=ill$Color_Channel

save(file="tmp.rda", list=c("gprobes","ill", "no.info")) 

numplus=numeric(length(gprobes))
numminus=numeric(length(gprobes))
full.strand=character(length(gprobes))
full.start=numeric(length(gprobes))
full.end=numeric(length(gprobes))
full.probes=gprobes


#Just the first 24 entries - these are the useful chromosomes
for (i in 1:24) {
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

  in.chrom=which((as.character(seqnames(gprobes))==seqname) & !no.info)
  ##Source seq not empty, and chromosome matches one we are probing

  for (j in in.chrom) {   

    probe=DNAString(ill$pattern[j])
    rc.probe=reverseComplement(probe)
    cat(">> Finding probe", j, "...\n")
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
}









                        
