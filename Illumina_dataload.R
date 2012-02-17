##To make the CSV file
##perl ~/Work/Analysis/Illumina/Illumina_MATLAB_CSV_Converter5.pl New_unnorm.csv 
##perl ~/Work/Analysis/Illumina/Illumina_MATLAB_CSV_Converter5.pl Second_unnorm_data.csv 
##perl ~/Work/Analysis/Illumina/ucsc_island_dist1.pl New_unnorm_mat_targer.csv
##perl ~/Work/Analysis/Illumina/rafa_island_dist1.pl New_unnorm_mat_targer_ucscisl.csv


##Load CSV file
##Custom set of probes and old set
probes<-read.csv("New_norm_mat_targer_ucscisl_hmmisl.csv",stringsAsFactors=FALSE)
##Original GoldenGate Methylation assay
##old_probes<-read.csv("Old_Illumina_pri_mat_targer_ucscisl_hmmisl.csv",stringsAsFactors=FALSE)


##Load Data file
raw<-read.csv("New_norm_mat_data.csv",stringsAsFactors=FALSE)
##raw_old<-read.csv("Old_Illumina_pri_mat_data.csv",stringsAsFactors=FALSE)
raw2<-read.csv("Second_norm_data_mat_data.csv",stringsAsFactors=FALSE)
##Get rid of last anomolous line and first line of second set - to merge smoothly
fuseraw=cbind(raw[-length(raw)],raw2[,-1])


##Load Sample file
##sampy2<-read.csv("Second_norm_data_mat_sample_anno.csv",stringsAsFactors=FALSE)
##sampy<-read.csv("New_norm_mat_sample_anno.csv",stringsAsFactors=FALSE)
##Again - combine runs 1 and 2
##fusesampy=rbind(sampy,sampy2)
##New manual sample file
fusesampy=read.csv("all_sample_anno1.csv", stringsAsFactors=FALSE)

##Original Goldengate
##old_sampy<-read.csv("Old_Illumina_pri_mat_sample.csv",stringsAsFactors=FALSE)

batch=read.csv("sample_plate_map1.csv", stringsAsFactors=FALSE)

fusesampy$Plate=batch$Plate[match(fusesampy$Sample_ID, batch$Sample_ID)]

data_norm=list(avg=as.matrix(fuseraw[,seq(2,ncol(fuseraw)-1,by=5)]), green=as.matrix(fuseraw[,seq(3,ncol(fuseraw)-1, by=5)]), red=as.matrix(fuseraw[,seq(4, ncol(fuseraw)-1, by=5)]),
  beads=as.matrix(fuseraw[,seq(5,ncol(fuseraw)-1, by=5)]), samp=fusesampy, probes=probes)

##CHARM probe loc load
##load("CHARM_nfo.rda")

##Load CHARM Data
##load("CHARM_data.rda")

##Unnormalized
##Load CSV file

##Load Data file
raw<-read.csv("New_unnorm_mat_data.csv",stringsAsFactors=FALSE)
##raw_old<-read.csv("Old_Illumina_pri_mat_data.csv",stringsAsFactors=FALSE)
raw2<-read.csv("Second_unnorm_data_mat_data.csv",stringsAsFactors=FALSE)
##Get rid of last anomolous line and first line of second set - to merge smoothly
fuseraw=cbind(raw[-length(raw)],raw2[,-1])


data=list(avg=as.matrix(fuseraw[,seq(2,ncol(fuseraw)-1,by=5)]), green=as.matrix(fuseraw[,seq(3,ncol(fuseraw)-1, by=5)]), red=as.matrix(fuseraw[,seq(4, ncol(fuseraw)-1, by=5)]),
  beads=as.matrix(fuseraw[,seq(5,ncol(fuseraw)-1, by=5)]), samp=fusesampy, probes=probes)

save(data,data_norm,file="good_data.rda")





