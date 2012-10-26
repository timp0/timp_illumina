##Initialize 450k Software in general - general case stuff gets loaded, no data
if (!exists("codedir")) {
  codedir=getwd()
}

source(file.path(codedir,"450k_timp_plotting.R"))
source(file.path(codedir,"450k_timp_tests.R"))
##source(file.path(codedir, "Rafa_functions.R"))


##Experimental file location
expdatapath="/thumper2/feinbergLab/core/arrays/illumina"

##Current 450k package name
##Kasper says ignore warnings about "contains no R code"
##This is the local(expanded) version of minfi
library(minfi)
library(minfiLocal)

##This will list commands
##ls("package:minfi")

##Get probe/snp map
load(file.path(codedir, "timp_illumina_data", "probe_obj_final.rda"))

##Probes with no problems SNPs
good.probes=values(gprobes)$minfi.idx[(!values(gprobes)$sbe.snp.boo)&
  (!values(gprobes)$boo.snps)&(values(gprobes)$single.hyb)&
  (!values(gprobes)$g.site.snp.boo)]

