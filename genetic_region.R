##Run using R CMD BATCH findseq1.R
##Finding chromosomal locations of sequences




    





##START MAIN

setwd("~/Big_Data/Illumina_Bead/Analysis/")

##Load libraries in
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)

##
source("~/Work/Analysis/Illumina/Illumina_density1b.R")
source("~/Work/Analysis/Illumina/genetic_plots1.R")
source("~/Work/Analysis/Illumina/Illumina_controls1.R")
##source("~/Work/Analysis/Illumina/Illumina_sample1b.R")

##Load CSV file
##Custom set of probes and old set
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
raw2<-read.csv("Second_norm_data_mat_data.csv",stringsAsFactors=FALSE)



##Get Avg Data
avg_data<-raw[,seq(2,ncol(raw)-1,by=5)]
avg_data2<-raw2[,seq(2,ncol(raw2)-1,by=5)]
old_avg_data<-raw_old[,seq(2,ncol(raw_old)-1,by=5)]
##Combine runs 1 and 2 of the Custom Array
fuseavg<-cbind(avg_data,avg_data2)
##Load Sample file
sampy2<-read.csv("Second_norm_data_mat_sample.csv",stringsAsFactors=FALSE)
sampy<-read.csv("New_norm_mat_sample.csv",stringsAsFactors=FALSE)
##Again - combine runs 1 and 2
fusesampy=rbind(sampy,sampy2)
old_sampy<-read.csv("Old_Illumina_pri_mat_sample.csv",stringsAsFactors=FALSE)

##CHARM probe loc load
##load("CHARM_nfo.rda")

##Load CHARM Data
##load("CHARM_data.rda")

class_stuff=data.frame(tissue=c(0,1,rep(2,4), rep(3,5),rep(4,2), rep(5,2), rep(6,2), rep(7,10),rep(8,5)),
  type=c(0,0,-2, -1, 1, 2, -3,-1,1,2,3,-1,1,-1,1,-1,1,-7,-4, -2, 1, 2, 3, 4, 5, 6, 7,-2,-1, 1, 2, 3),
  nam=c("Controls","Cross_controls","Breast DCIS Normal","Breast Normal","Breast Tumor", "Breast DCIS Tumor",
    "Colon Seq Normal", "Colon Normal", "Colon Tumor", "Tubular Adenoma", "Colon Seq Tumor",
    "Lung Normal", "Lung Tumor", "Ovary Normal", "Ovary Tumor", "Wilms Normal", "Wilms Tumor",
    "Thyroid PTC Normal", "Thyroid FVPTC Normal", "Thyroid FA Normal",
    "Thyroid Adenomatoid", "Thyroid FA Tumor", "Thryoid FC Tumor", "Thyroid FVPTC Tumor", "Thryoid HA Tumor", "Thyroid HC Tumor", "Thyroid PTC Tumor",
    "Pancreas NET Normal", "Pancreas IPMN Normal", "Pancreas IPMN Tumor", "Pancreas NET Tumor", "Pancreas LCC Tumor"),
  coloring=c("black","orange", rep("coral1",2), rep("coral4",2), rep("cadetblue1",2), rep("cadetblue4",3), "goldenrod1","goldenrod4","seagreen1","seagreen4","plum1","plum4",
    rep("firebrick1", 3), rep("firebrick4", 7), rep("lightsteelblue1", 2), rep("lightsteelblue4", 3)),stringsAsFactors=FALSE)


##Plot new data all sets

##plotsubdTN(probes,genes,avg_data,sampy,"Movie/full_big_hist.pdf")


##Plot Hector Regions
topprobes=read.csv("From_Hector_20methylationRegions.csv",stringsAsFactors=F)
goody=is.element(probes$Illumina_ID,topprobes$TargetID)


plotsubdTN(probes,genes,fuseavg,fusesampy,"Movie/new_topprobes_hist.pdf",interest=goody)



#plotClass(probes,genes,avg_data,sampy,c(6:7),"Movie/lung_R.pdf",interest=goody)
#plotClass(probes,genes,avg_data,sampy,c(2:11),"Movie/apple1.pdf",interest=goody)
