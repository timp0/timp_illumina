##Run using R CMD BATCH findseq1.R
##Finding chromosomal locations of sequences




    





##START MAIN

setwd("~/Big_Data/Illumina_Bead/Analysis/")

##Load libraries in
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)

##
source("~/Work/Analysis/Illumina/Illumina_density1.R")
source("~/Work/Analysis/Illumina/genetic_plots1.R")
source("~/Work/Analysis/Illumina/Illumina_controls1.R")
source("~/Work/Analysis/Illumina/Illumina_sample1.R")

##Load CSV file
probes<-read.csv("New_norm_mat_targer_ucscisl_hmmisl.csv",stringsAsFactors=FALSE)
old_probes<-read.csv("Old_Illumina_pri_mat_targer_ucscisl_hmmisl.csv",stringsAsFactors=FALSE)
##Load UCSC Isl file
ucsc_isl<-read.delim("ucsc_cpgisl.txt",stringsAsFactors=FALSE)
##Load HMM Isl file
hmm_isl<-read.csv("hmm_isl1.txt",stringsAsFactors=FALSE)  
##Load RefSeq Gene file
genes<-read.delim("ref_genes.txt",stringsAsFactors=FALSE)
##Load Data file
raw<-read.csv("New_norm_mat_data.csv",stringsAsFactors=FALSE)
raw_old<-read.csv("Old_Illumina_pri_mat_data.csv",stringsAsFactors=FALSE)
##Get Avg Data
avg_data<-raw[,seq(2,ncol(raw)-1,by=5)]
old_avg_data<-raw_old[,seq(2,ncol(raw_old)-1,by=5)]
##Load Sample file
sampy<-read.csv("New_norm_mat_sample.csv",stringsAsFactors=FALSE)
old_sampy<-read.csv("Old_Illumina_pri_mat_sample.csv",stringsAsFactors=FALSE)

##CHARM probe loc load
##load("CHARM_nfo.rda")

##Load CHARM Data
##load("CHARM_data.rda")

##Define Classes/Colors
class_stuff=data.frame(nums=c(1,2,3,4,5,6,7,8,9,10,11),
  nam=c("Controls","Breast Tumor","Breast Normal","Colon Tumor", "Colon Normal","Lung Tumor","Lung Normal","Ovary Tumor","Ovary Normal","Wilms Tumor","Wilms Normal"),
  coloring=c("black","coral1","coral4","cadetblue1","cadetblue4","goldenrod1","goldenrod4","seagreen1","seagreen4","plum1","plum4"),stringsAsFactors=FALSE)


##Plot new data all sets

##plotsubdTN(probes,genes,avg_data,sampy,"Movie/full_big_hist.pdf")


##Plot Hector Regions
topprobes=read.csv("From_Hector_20methylationRegions.csv",stringsAsFactors=F)
goody=is.element(probes$Illumina_ID,topprobes$TargetID)
##plotsubdTN(probes,genes,avg_data,sampy,"Movie/topprobes_hist.pdf",interest=goody)



plotClass(probes,genes,avg_data,sampy,c(6:7),"Movie/lung_R.pdf",interest=goody)
plotClass(probes,genes,avg_data,sampy,c(2:11),"Movie/apple1.pdf",interest=goody)
