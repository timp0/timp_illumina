##Put us in the right dir for output
setwd("~/Data/Infinium/122511_analysis/")

##Experimental file location
expdatapath="/thumper2/feinbergLab/core/arrays/illumina/"

##Current 450k package name
##Kasper says ignore warnings about "contains no R code"
##This is the local(expanded) version of minfi
library(minfi)
library(minfiLocal)

##This will list commands
##ls("package:minfi")

##Needs this data object apparently?  Let's test!!
##data(IlluminaHumanMethylation450kmanifest)


##Read in plates
plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))

##Read in data
##RGset=read.450k.exp(base=expdatapath, targets=plates)

dna.org=read.csv(file="DNA_Org2.csv", stringsAsFactors=F)

##The met plates have only 24 samples which aren't in the dna.org - samples
##from Yun which are from Marcia - don't care about them - they are terrible
##Samples anyway

z=plates$Sample.ID %in% dna.org$X450k_ID

y=plates$Sample.ID[!z]
##Only samples which don't match are the stupid percentage samples and Yun
##samples

##ok - am going to regenerate a better DNA_Org chart

dna.base=data.frame(sample.id=dna.org$ID., matched.id=dna.org$Paired.With,
  state=dna.org$State, tissue=dna.org$Tissue, type=dna.org$Type,
  tissue.location=dna.org$Freezer.Location, 
  sample.group=dna.org$Collection, sample.source=dna.org$Source,
  sample.notes=dna.org$Notes,
  path.report=dna.org$Path.report, path.qc=dna.org$Path.QC,
  flag.p53=dna.org$Flag_P53,
  rafa.type=dna.org$Rafa_Type, rafa.progression=dna.org$Rafa_Progression,
  patient.stats=dna.org$Patient.Stats, patient.age=dna.org$Age,
  patient.surgery.date=dna.org$Date.of.Surgery,
  sample.received.date=dna.org$Date.Received,
  dna.purification.date=dna.org$DNA.Date, dna.location=dna.org$DNA.Location,
  dna.well=dna.org$DNA.Well, dna.260.280=dna.org$DNA.260.280,
  dna.260.230=dna.org$DNA.260.230,
  dna.concentration=dna.org$Stock.DNA.Concentration,
  dna.left=dna.org$DNA_Left,
  custom.illumina.id=dna.org$Illumina_Data_ID,
  kun.id=dna.org$Kun_Id, kun.random=dna.org$Kun_Random,
  agilent.id=dna.org$Agilent_ID, 
  minfi.id=dna.org$X450k_ID, minfi.run=dna.org$X450k_Run,
  stringsAsFactors=F)

samp.corr=match(plates$Sample.ID, dna.base$minfi.id)

##Remove non matching samples - aka Yun's stuff, and % stuff
good.corr=!is.na(samp.corr)

dna.base$minfi.run[samp.corr[good.corr]]=plates$Plate[good.corr]

##a couple other corrections

dna.base$kun.id[dna.base$sample.id=="517N"]=""
dna.base$kun.id[dna.base$sample.id=="517T"]=""

dna.base$agilent.id[grepl("Control",dna.base$sample.id)]=""

##And I need to add a column for samples which were BS-Sequenced and those
##IDs

dna.base$bsseq.id=character(dim(dna.base)[1])

dna.base$bsseq.id[dna.base$custom.illumina.id=="Adenoma_T17"]="Adn1"
dna.base$bsseq.id[dna.base$custom.illumina.id=="Adenoma_T2"]="Adn2"

dna.base$bsseq.id[dna.base$sample.id=="499N"]="4N"
dna.base$bsseq.id[dna.base$sample.id=="499T"]="4T"
dna.base$bsseq.id[dna.base$sample.id=="536N"]="5N"
dna.base$bsseq.id[dna.base$sample.id=="536T"]="5T"
dna.base$bsseq.id[dna.base$sample.id=="614N"]="6N"
dna.base$bsseq.id[dna.base$sample.id=="614T"]="6T"


##Good - now I have a table
write.table(dna.base, file="dnabase2a.txt", quote=F, row.names=F, sep="\t")
save(file="dnabase2a.rda", list="dna.base")

##

