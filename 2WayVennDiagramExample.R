install.packages("VennDiagram",repos="http://cran.us.r-project.org")
library(VennDiagram)

source("/thumper2/feinbergLab/personal/mmulthau/Diabetes/functions/find2GeneListOverlap.R")
source("/home/bst/student/pmurakam/feinberg/CharmFiles/cis-genome.R")

#Get two gene lists
#Here, my two gene lists are the top 200 asthma results at regionfinder cutoff 2 and the top 100 results at regionfinder cutoff 4
load("/thumper2/feinbergLab/personal/mmulthau/Asthma/OldRafaCutoff2/foundRegions.rda")
Cutoff2Regions=matchGenes(tab,build="hg19")
Cutoff2Regions=Cutoff2Regions[1:200,]
load("/thumper2/feinbergLab/personal/mmulthau/Asthma/OldRafaCutoff4/foundRegions.rda")
Cutoff4Regions=matchGenes(tab,build="hg19")
Cutoff4Regions=Cutoff4Regions[1:100,]

#Get counts of overlaps and uniques for the two gene lists
overlaps=find2GeneListOverlap(Cutoff2Regions$name,Cutoff4Regions$name)

#Create VennDiagram
pdf("Cutoffs2and4VennDiagram.pdf")
draw.pairwise.venn(overlaps[1], overlaps[2], overlaps[3], c("Cutoff=2", "Cutoff=4"), fill = c("blue", "red"))
dev.off()