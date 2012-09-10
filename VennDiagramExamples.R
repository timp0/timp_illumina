install.packages("VennDiagram",repos="http://cran.us.r-project.org")
library(VennDiagram)



#Two way venn diagram examples
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




#Three way venn diagram example - I don't think there is a way to get this to scale. Apparently the math is hard.
source("/home/bst/student/pmurakam/feinberg/CharmFiles/cis-genome.R")
source("/thumper2/feinbergLab/personal/mmulthau/Diabetes/functions/findGeneListOverlaps.R")

#Get three gene lists
#Here, my two gene lists are the top 200 asthma results at regionfinder cutoff 2 and the top 100 results at regionfinder cutoff 4
load("/thumper2/feinbergLab/personal/mmulthau/Asthma/OldRafaCutoff2/foundRegions.rda")
Cutoff2Regions=matchGenes(tab,build="hg19")
Cutoff2Regions=Cutoff2Regions[1:200,]
load("/thumper2/feinbergLab/personal/mmulthau/Asthma/OldRafaCutoff4/foundRegions.rda")
Cutoff4Regions=matchGenes(tab,build="hg19")
Cutoff4Regions=Cutoff4Regions[1:100,]
load("/thumper2/feinbergLab/personal/mmulthau/Asthma/foundRegions.rda")
Cutoff5Regions=matchGenes(tab,build="hg19")

#Get counts of overlaps and uniques for the two gene lists
overlaps=find3GeneListOverlap(Cutoff2Regions$name,Cutoff4Regions$name,Cutoff5Regions$name)

#Create VennDiagram - note that the areas supplied aren't just the straight readouts of find3GeneListOverlaps
pdf("Cutoffs245VennDiagram.pdf")
draw.triple.venn(
  area1 = nrow(Cutoff2Regions),
  area2 = nrow(Cutoff4Regions),
  area3 = nrow(Cutoff5Regions),
  n12 = overlaps[4]+overlaps[7],
  n13 = overlaps[5]+overlaps[7],
  n23 = overlaps[6]+overlaps[7],
  n123 = overlaps[7],
  category = c("Cutoff2", "Cutoff4", "Cutoff5"),
  fill = c("blue", "red", "green"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "green")
  )
dev.off()

pdf("Cutoffs245VennDiagram.pdf")
draw.triple.venn(overlaps[1], overlaps[2], overlaps[3], overlaps[4]+overlaps[7], overlaps[6]+overlaps[7], overlaps[5]+overlaps[7], overlaps[7], c("First", "Second", "Third"))
dev.off()
