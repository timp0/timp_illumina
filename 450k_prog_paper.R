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

cancer.census=read.csv("~/Dropbox/Data/Genetics/Infinium/091013_writing/cancer_gene_census1.csv")

load("~/Dropbox/Data/Genetics/MethSeq/072111_blocks/gene_island.rda")

z=match(cancer.census$Symbol, values(refseq.genes)$gene.name)

y=grepl(cancer.census$Symbol[4], hugo$all.symbol)

x=match(cancer.census$GeneID, hugo$Entrez.Gene.ID)

aaa=match(cancer.census$Symbol[is.na(x)], hugo$Approved.Symbol)




mutdir="~/Dropbox/Data/Genetics/Mutations/092512_dl"
#Mutation tables
load(file.path(mutdir, "muts.rda"))

