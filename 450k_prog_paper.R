codedir="~/Code/timp_illumina"
plotdir="~/Dropbox/Data/Genetics/Infinium/081313_analysis"
filedir="/mithril/homes/timp/LData/Genetics/Infinium/071313_analysis"
expdatapath="/mithril/Data/Infinium"

source(file.path(codedir,"450k_general_init.R"))

require(doMC)
registerDoMC()
cores=4
options(cores=cores)

if (file.exists(file.path(filedir, "cancer.rda"))) {
  load(file.path(filedir, "cancer.rda"))
  
} else {
  
  ##Gets only plate files and loads in basic stuff
  source(file.path(codedir,"450k_cancer_loadin.R"))
  plates=plates[!(plates$Tissue %in% c("cell.line", "mix", "urine", "saliva", "blood", "placenta", "na")),]
  ##Filter data we don't have
  plates=plates[plates$Basename!="character(0)",]
  plates=plates[plates$Phenotype!="bad",]
  
  dat <- read.450k.exp(targets=plates, verbose=TRUE) 
  dat=preprocessQuantile(dat)

  ##Add CpG Island anno info
  dat=timp.probeanno(dat)
  
  save(list="dat", file=file.path(filedir, "cancer.rda"), compress="gzip")
}

pd=colData(dat)

load("~/Dropbox/Data/Genetics/Infinium/091013_writing/rafa_blocks.rda")
load("~/Dropbox/Data/Genetics/MethSeq/072111_blocks/gene_island.rda")

##First, let's get hugo symbols for all refseq


##Ok - find all genes in blocks/genes out of rafablocks
breast.block=blocks[[1]]$table
breast.block.gr=GRanges(seqnames=breast.block$chr, ranges=IRanges(start=breast.block$start, end=breast.block$end),
    delta=breast.block$value)

colon.block=blocks[[5]]$table
colon.block.gr=GRanges(seqnames=colon.block$chr, ranges=IRanges(start=colon.block$start, end=colon.block$end),
    delta=colon.block$value)

lung.block=blocks[[14]]$table
lung.block.gr=GRanges(seqnames=lung.block$chr, ranges=IRanges(start=lung.block$start, end=lung.block$end),
    delta=lung.block$value)

pancreas.block=blocks[[15]]$table
pancreas.block.gr=GRanges(seqnames=pancreas.block$chr, ranges=IRanges(start=pancreas.block$start, end=pancreas.block$end),
    delta=pancreas.block$value)

thyroid.block=blocks[[18]]$table
thyroid.block.gr=GRanges(seqnames=thyroid.block$chr, ranges=IRanges(start=thyroid.block$start, end=thyroid.block$end),
    delta=thyroid.block$value)


##Trying synapse
library(synapseClient)
library(GenomicRanges)
synapseLogin('wtimp@jhu.edu', '1!HGUy0v20tn')
all.mut <- syn1710680 <- synGet(id='syn1710680')

##BRCA is breast file 2
##COADREAD is colorectal carcinoma file 4
##LUAD is Lung adenocarcinoma file 11
##PAAd is pancreas adenocarcinoma  file 14
##THCA is thyroid carcinoma file 19

breast.mut <- read.delim(file.path(all.mut$cacheDir, all.mut$files[[2]]), check.names=F)
breast.mut.gr=GRanges(seqnames=paste0('chr', breast.mut$Chromosome),
    ranges=IRanges(start=breast.mut$Start_Position, end=breast.mut$End_Position),
    strand=breast.mut$Strand)
        
colon.mut <- read.delim(file.path(all.mut$cacheDir, all.mut$files[[4]]), check.names=F)
colon.mut.gr=GRanges(seqnames=paste0('chr', colon.mut$Chromosome),
    ranges=IRanges(start=colon.mut$Start_Position, end=colon.mut$End_Position),
    strand=colon.mut$Strand)

lung.mut <- read.delim(file.path(all.mut$cacheDir, all.mut$files[[11]]), check.names=F)
lung.mut.gr=GRanges(seqnames=paste0('chr', lung.mut$Chromosome),
    ranges=IRanges(start=lung.mut$Start_Position, end=lung.mut$End_Position),
    strand=lung.mut$Strand)

pancreas.mut <- read.delim(file.path(all.mut$cacheDir, all.mut$files[[14]]), check.names=F)
pancreas.mut.gr=GRanges(seqnames=paste0('chr', pancreas.mut$Chromosome),
    ranges=IRanges(start=pancreas.mut$Start_Position, end=pancreas.mut$End_Position),
    strand=pancreas.mut$Strand)

thyroid.mut <- read.delim(file.path(all.mut$cacheDir, all.mut$files[[19]]), check.names=F)
thyroid.mut.gr=GRanges(seqnames=paste0('chr', thyroid.mut$Chromosome),
    ranges=IRanges(start=thyroid.mut$Start_Position, end=thyroid.mut$End_Position),
    strand=thyroid.mut$Strand)

##Mut table


mut.tab=data.frame(row.names=c("colon.mut", "breast.mut", "lung.mut", "pancreas.mut", "thyroid.mut"))

mut.tab$colon.block=c(sum(colon.mut.gr %over% colon.block.gr),
    sum(breast.mut.gr %over% colon.block.gr),
    sum(lung.mut.gr %over% colon.block.gr),
    sum(pancreas.mut.gr %over% colon.block.gr),
    sum(thyroid.mut.gr %over% colon.block.gr))

mut.tab$breast.block=c(sum(colon.mut.gr %over% breast.block.gr),
    sum(breast.mut.gr %over% breast.block.gr),
    sum(lung.mut.gr %over% breast.block.gr),
    sum(pancreas.mut.gr %over% breast.block.gr),
    sum(thyroid.mut.gr %over% breast.block.gr))

mut.tab$lung.block=c(sum(colon.mut.gr %over% lung.block.gr),
    sum(breast.mut.gr %over% lung.block.gr),
    sum(lung.mut.gr %over% lung.block.gr),
    sum(pancreas.mut.gr %over% lung.block.gr),
    sum(thyroid.mut.gr %over% lung.block.gr))

mut.tab$pancreas.block=c(sum(colon.mut.gr %over% pancreas.block.gr),
    sum(breast.mut.gr %over% pancreas.block.gr),
    sum(lung.mut.gr %over% pancreas.block.gr),
    sum(pancreas.mut.gr %over% pancreas.block.gr),
    sum(thyroid.mut.gr %over% pancreas.block.gr))

mut.tab$thyroid.block=c(sum(colon.mut.gr %over% thyroid.block.gr),
    sum(breast.mut.gr %over% thyroid.block.gr),
    sum(lung.mut.gr %over% thyroid.block.gr),
    sum(pancreas.mut.gr %over% thyroid.block.gr),
    sum(thyroid.mut.gr %over% thyroid.block.gr))

mut.tab[1,]=mut.tab[1,]/length(colon.mut.gr)
mut.tab[2,]=mut.tab[2,]/length(breast.mut.gr)
mut.tab[3,]=mut.tab[4,]/length(lung.mut.gr)
mut.tab[4,]=mut.tab[5,]/length(pancreas.mut.gr)
mut.tab[5,]=mut.tab[5,]/length(thyroid.mut.gr)

write.csv(mut.tab, file.path('~/Dropbox/Data/Genetics/Infinium/091013_writing', 'tcga_block_mut.csv'))

