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

##Assign Index num to probes
ill$pidx=1:dim(ill)[1]

#Initialize hit table

##Make GRanges of probes
gprobes=GRanges(seqname=paste("chr", ill$CHR, sep=""),
  strand=ill$Strand,ranges=IRanges(start=ill$MAPINFO, end=ill$MAPINFO))

values(gprobes)$name=ill$Name
values(gprobes)$pidx=ill$pidx
values(gprobes)$type=ill$Infinium_Design_Type

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


                        
