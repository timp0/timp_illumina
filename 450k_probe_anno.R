##Put us in the right dir for output
setwd("~/Data/Infinium/121311_analysis/")

##Current 450k package name
##Kasper says ignore warnings about "contains no R code"
##This is the local(expanded) version of minfi
library(minfi)
library(minfiLocal)

##This will list commands
##ls("package:minfi")

library(GenomicRanges)
library(Biostrings)

##Needs this data object apparently?  Let's test!!
##data(IlluminaHumanMethylation450kanno)


##From Chris -
chris.anno=read.csv("~/Data/Infinium/121311_analysis/450kannotationTableMEGA.csv", header=T, stringsAsFactors=F)

##Annotate the 65 rs probes with no info
no.info=chris.anno$Name[is.na(chris.anno$pos)]
## I went online and queried these rs##
dbsnp.query=read.csv(file="dbsnp_query1.csv", stringsAsFactors=F)
## Keep just the hg19 
dbsnp.query=dbsnp.query[dbsnp.query$a=="GRCh37.p5",]
##Make numbers into rsnumber for names to match
dbsnp.query$names=paste("rs", dbsnp.query$rs., sep="")
##SNP rs13369115 was merged into rs10155413
dbsnp.query$names[57]="rs13369115"
##Match to locs in probe table
q=match(dbsnp.query$names, chris.anno$Name)

chris.anno$Genome_Build[q]=37
chris.anno$chr[q]=dbsnp.query$chr
chris.anno$pos[q]=dbsnp.query$chrpos
chris.anno$avHet[q]=dbsnp.query$avghet
chris.anno$avHetSE[q]=dbsnp.query$s.e.het

##Change strand to "+" or "-"
chris.anno$Strand[chris.anno$Strand=="F"]="+"
chris.anno$Strand[chris.anno$Strand=="R"]="-"
chris.anno$Strand[chris.anno$Strand==""]="*"

##Need to calc dist to CpG Island, dist to genes
##Load islands and genes object
load("~/Data/Genetics/072111_blocks/gene_island.rda")

##Make GRanges of probes
gprobes=GRanges(seqname=paste("chr", chris.anno$chr, sep=""),
  strand=chris.anno$Strand,
  ranges=IRanges(start=chris.anno$pos, end=chris.anno$pos))

values(gprobes)$name=chris.anno$Name

##This is needed so nearest will work, seems to be an issue with the version
##Of GRanges which made the object??
refseq.genes=updateObject(refseq.genes)
ucsc.isl=updateObject(ucsc.isl)

##Find nearest genes to each probe(including overlap)
z=nearest(gprobes, refseq.genes)

##Put gene name and dist to gene
values(gprobes)$nearest.gene=values(refseq.genes[z])$gene.name
values(gprobes)$dist.gene=width(pgap(ranges(gprobes), ranges(refseq.genes[z])))



##Find nearest UCSC CpG islands to each probe (including overlap)
z=nearest(gprobes, ucsc.isl)

##UCSC island index, and dist to ucsc isl
values(gprobes)$nearest.island.index=z
values(gprobes)$dist.island=width(pgap(ranges(gprobes), ranges(ucsc.isl[z])))

##Get DMR info

##Get Block/LOCK/LAD info

##Pass in SNP info
##Load snp granges
##load("~/Big_Data/Genetics/121511_dbSNP/snp135.rda")


values(gprobes)$il.probe.snps=chris.anno$Probe_SNPs
values(gprobes)$il.probe.snps.10=chris.anno$Probe_SNPs_10

##Booleans from probe design origins
values(gprobes)$random.loci=!is.na(chris.anno$Random_Loci)
values(gprobes)$M27k.array=!is.na(chris.anno$Methyl27_Loci)
values(gprobes)$fantom=!(chris.anno$Phantom)
values(gprobes)$fantom.nfo=chris.anno$Phantom
values(gprobes)$enhancer=!is.na(chris.anno$Enhancer)
values(gprobes)$dnase.hypersensitivity=!is.na(chris.anno$DHS)

##Type of probe info
values(gprobes)$probe.type=chris.anno$Infinium_Design_Type
values(gprobes)$probe.nxt.base=chris.anno$Next_Base
values(gprobes)$probe.color.chan=chris.anno$Color_Channel



  

but altered by me:
##Identify Sex chromosome probes
x <- IlluminaHumanMethylation450kanno$Chr_37=="X" #11232 probes
y <- IlluminaHumanMethylation450kanno$Chr_37=="Y" #416 probes
autosomal <- (!x&!y)
sex.probes = IlluminaHumanMethylation450kanno$Name[(x|y)] #just the probe names


##From Chris - array annotation
#---Load the array annotation files--#
int <- read.delim("~/Data/Infinium/100411_analysis/Probe_intersect_SNP_20111006",
                  header=FALSE, stringsAsFactors=FALSE) #probe coordinates intersecting All132SNP table
anno <-  read.table("~/Data/Infinium/100411_analysis/450kanno.txt",header=TRUE, stringsAsFactors=FALSE) #450k annotation

rownames(anno)=anno$Name  #-add annotation info
colnames(int)=c("chr", "start", "end", "Name") #header to intersectio file

450kprobes=data.frame(name=IlluminaHumanMethylation450kanno$Name, build=IlluminaHumanMethylation450kanno$Genome_Build,
  chr=IlluminaHumanMethylation450kanno$Chr_37, cpgpos=anno$pos,
  probe.snp=anno$Probe_SNPs, probe.snp.10=anno$Probe_SNPs_10
  probe.type=anno$Infinium_Design_Type, probe.nextbase=anno$Next_Base

