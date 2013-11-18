codedir="~/Code/timp_illumina"
plotdir="~/Dropbox/Data/Genetics/Infinium/081313_analysis"
filedir="/mithril/homes/timp/LData/Genetics/Infinium/071313_analysis"
expdatapath="/mithril/Data/Infinium"

source(file.path(codedir,"450k_general_init.R"))

library(qvalue)
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

collapse.cluster <- cpgCollapse(dat)
probey=rowData(dat)


load("~/Dropbox/Data/Genetics/Infinium/091013_writing/rafa_blocks.rda")
load("~/Dropbox/Data/Genetics/MethSeq/072111_blocks/gene_island.rda")

##Ok - find all genes in blocks/genes out of rafablocks

##Rafa blocks are sometimes of length 0, I assume because he has set the end to same as start?  Let's see if using end+1 helps
breast.block=blocks[[1]]$table
breast.block$q.value=qvalue(breast.block$p.value)$qvalues
breast.block.gr=foreach(i=1:dim(breast.block)[1], .combine='c') %dopar% {
    GRanges(seqnames=breast.block$chr[i], ranges=IRanges(
                                              start=min(start(probey[collapse.cluster$blockInfo$indexes[[breast.block$indexStart[i]]]])),
                                              end=max(end(probey[collapse.cluster$blockInfo$indexes[[breast.block$indexEnd[i]]]]))),
            delta=breast.block$value[i], q.value=breast.block$q.value[i])
}
sum(values(breast.block.gr)$q.value<.05)
length(breast.block.gr)
breast.block.gr=breast.block.gr[values(breast.block.gr)$q.value<.05]



#Rafa method
#ind <- colon.block$indexStart
#st=granges(collapse.cluster$object)[ind,] ##this is the left side of the block

colon.block=blocks[[5]]$table
colon.block$q.value=qvalue(colon.block$p.value)$qvalues
colon.block.gr=foreach(i=1:dim(colon.block)[1], .combine='c') %dopar% {
    GRanges(seqnames=colon.block$chr[i], ranges=IRanges(
                                             start=min(start(probey[collapse.cluster$blockInfo$indexes[[colon.block$indexStart[i]]]])),
                                             end=max(end(probey[collapse.cluster$blockInfo$indexes[[colon.block$indexEnd[i]]]]))),
                        delta=colon.block$value[i], q.value=colon.block$q.value[i])
}
sum(values(colon.block.gr)$q.value<.05)
length(colon.block.gr)
colon.block.gr=colon.block.gr[values(colon.block.gr)$q.value<.05]


pdf(file.path(plotdir, "pvalhist.pdf"))
print(ggplot(colon.block, aes(x=p.value))+geom_histogram())
dev.off()


lung.block=blocks[[14]]$table
lung.block$q.value=qvalue(lung.block$p.value)$qvalues
lung.block.gr=foreach(i=1:dim(lung.block)[1], .combine='c') %dopar% {
    GRanges(seqnames=lung.block$chr[i], ranges=IRanges(
                                           start=min(start(probey[collapse.cluster$blockInfo$indexes[[lung.block$indexStart[i]]]])),
                                           end=max(end(probey[collapse.cluster$blockInfo$indexes[[lung.block$indexEnd[i]]]]))),
                        delta=lung.block$value[i], q.value=lung.block$q.value[i])
}
sum(values(lung.block.gr)$q.value<.05)
length(lung.block.gr)
lung.block.gr=lung.block.gr[values(lung.block.gr)$q.value<.05]


pancreas.block=blocks[[15]]$table
pancreas.block$q.value=qvalue(pancreas.block$p.value)$qvalues
pancreas.block.gr=foreach(i=1:dim(pancreas.block)[1], .combine='c') %dopar% {
    GRanges(seqnames=pancreas.block$chr[i], ranges=IRanges(
                                           start=min(start(probey[collapse.cluster$blockInfo$indexes[[pancreas.block$indexStart[i]]]])),
                                           end=max(end(probey[collapse.cluster$blockInfo$indexes[[pancreas.block$indexEnd[i]]]]))),
            delta=pancreas.block$value[i], q.value=pancreas.block$q.value[i])
}
sum(values(pancreas.block.gr)$q.value<.05)
length(pancreas.block.gr)
pancreas.block.gr=pancreas.block.gr[values(pancreas.block.gr)$q.value<.05]


thyroid.block=blocks[[18]]$table
thyroid.block$q.value=qvalue(thyroid.block$p.value)$qvalues
thyroid.block.gr=foreach(i=1:dim(thyroid.block)[1], .combine='c') %dopar% {
    GRanges(seqnames=thyroid.block$chr[i], ranges=IRanges(
                                           start=min(start(probey[collapse.cluster$blockInfo$indexes[[thyroid.block$indexStart[i]]]])),
                                           end=max(end(probey[collapse.cluster$blockInfo$indexes[[thyroid.block$indexEnd[i]]]]))),
            delta=thyroid.block$value[i], q.value=thyroid.block$q.value[i])
}
sum(values(thyroid.block.gr)$q.value<.05)
length(thyroid.block.gr)
thyroid.block.gr=thyroid.block.gr[values(thyroid.block.gr)$q.value<.05]


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

library(BSgenome.Hsapiens.UCSC.hg19)
chrnames=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
      "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
      "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
hg19.len=foreach(i=1:24, .combine="sum") %dopar% {
    len=nchar(Hsapiens[[chrnames[i]]])
    return(as.numeric(len))
}

mut.tab$tot.muts=c(length(colon.mut.gr), length(breast.mut.gr), length(lung.mut.gr), length(pancreas.mut.gr), length(thyroid.mut.gr))

mut.tab[6,]=c(sum(width(colon.block.gr))/hg19.len, sum(width(breast.block.gr))/hg19.len, sum(width(lung.block.gr))/hg19.len,
           sum(width(pancreas.block.gr))/hg19.len, sum(width(thyroid.block.gr))/hg19.len, 0)


write.csv(mut.tab, file.path('~/Dropbox/Data/Genetics/Infinium/091013_writing', 'tcga_block_mut.csv'))



## Block variation

##Simple Colon
sub=dat[,(pd$Tissue=="colon")]
subpd=colData(sub)

subpd$anno=NA
subpd$anno[subpd$Tissue=="colon" & subpd$Phenotype=="normal"]="colon.normal"
subpd$anno[subpd$Tissue=="colon" & subpd$Phenotype=="adenoma"]="colon.adenoma"
subpd$anno[subpd$Tissue=="colon" & subpd$Phenotype=="cancer"]="colon.cancer"
subpd$anno[grepl("metastasis", subpd$Notes)]="colon.metastasis"

sub=sub[,!is.na(subpd$anno)]
colData(sub)=subpd[!is.na(subpd$anno),]

colon.stats=range.stats(sub, colon.block.gr, logit=F)
write.csv(colon.stats, file.path('~/Dropbox/Data/Genetics/Infinium/091013_writing', 'colon_block_stats.csv'))

compname="pancreassimple"
sub=dat[,which(pd$Tissue=="pancreas"|grepl("pancreas", pd$Notes))]

subpd=colData(sub)

subpd$anno=NA
subpd$anno[subpd$Notes=="adenocarcinoma"]="pancreas.cancer"
subpd$anno[subpd$Phenotype=="metastasis"]="metastasis"
subpd$anno[subpd$Phenotype=="normal"&subpd$Tissue=="pancreas"]="pancreas.normal"
subpd$anno[subpd$Phenotype=="adenoma"]="IPMN"

sub=sub[,!is.na(subpd$anno)]
colData(sub)=subpd[!is.na(subpd$anno),]

pancreas.stats=range.stats(sub, pancreas.block.gr, logit=F)
write.csv(pancreas.stats, file.path('~/Dropbox/Data/Genetics/Infinium/091013_writing', 'pancreas_block_stats.csv'))

compname="breast"
sub=dat[,pd$Tissue=="breast"]

breast.stats=range.stats(sub, breast.block.gr, logit=F)
write.csv(breast.stats, file.path('~/Dropbox/Data/Genetics/Infinium/091013_writing', 'breast_block_stats.csv'))

compname="lung"
sub=dat[,pd$Tissue=="lung"]

lung.stats=range.stats(sub, lung.block.gr, logit=F)
write.csv(lung.stats, file.path('~/Dropbox/Data/Genetics/Infinium/091013_writing', 'lung_block_stats.csv'))

compname="thyroid"
sub=dat[,pd$Tissue=="thyroid"]

thyroid.stats=range.stats(sub, thyroid.block.gr, logit=F)
write.csv(thyroid.stats, file.path('~/Dropbox/Data/Genetics/Infinium/091013_writing', 'thyroid_block_stats.csv'))
