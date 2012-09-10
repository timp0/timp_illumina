#install.packages("VennDiagram",repos="http://cran.us.r-project.org")
#library(VennDiagram)


#This function takes two lists and returns counts of the unique entries in list1, the unique entries in list2, and the entries present in both lists
#Note: this function removes duplicate gene names from both lists. Often, DMR lists have multiple DMRs corresponding to the same gene
find2GeneListOverlap=function(list1,list2)
{
  list1=unique(list1)
  list2=unique(list2)
  
  list1unique=0
  list2unique=0
  listoverlap=0
  
  for(i in 1:length(list1))
  {
    if (list1[i]%in%list2) 
    {
      listoverlap=listoverlap+1
    } else list1unique=list1unique+1
  }
  
  list2unique=length(list2)-listoverlap
  
  return(c(list1unique,list2unique,listoverlap))
}

#Example
#source("/home/bst/student/pmurakam/feinberg/CharmFiles/cis-genome.R")

#Get two gene lists
#Here, my two gene lists are the top 200 asthma results at regionfinder cutoff 2 and the top 100 results at regionfinder cutoff 4
#load("/thumper2/feinbergLab/personal/mmulthau/Asthma/OldRafaCutoff2/foundRegions.rda")
#Cutoff2Regions=matchGenes(tab,build="hg19")
#Cutoff2Regions=Cutoff2Regions[1:200,]
#load("/thumper2/feinbergLab/personal/mmulthau/Asthma/OldRafaCutoff4/foundRegions.rda")
#Cutoff4Regions=matchGenes(tab,build="hg19")
#Cutoff4Regions=Cutoff4Regions[1:100,]

#Get counts of overlaps and uniques for the two gene lists
#overlaps=find2GeneListOverlap(Cutoff2Regions$name,Cutoff4Regions$name)

#Create VennDiagram
#pdf("Cutoffs2and4VennDiagram.pdf")
#draw.pairwise.venn(overlaps[1], overlaps[2], overlaps[3], c("Cutoff=2", "Cutoff=4"), fill = c("blue", "red"))
#dev.off()

#This function finds the number of unique entries and overlaps for three gene lists
find3GeneListOverlap=function(list1,list2,list3)
{
  list1=unique(list1)
  list2=unique(list2)
  list3=unique(list3)
  
  list1unique=0
  list2unique=0
  list3unique=0
  listoverlap12=0
  listoverlap13=0
  listoverlap23=0
  listoverlap123=0
  
  #Determine all the overlaps with list 1
  for(i in 1:length(list1)) {
    
    if (list1[i]%in%list2&list1[i]%in%list3) 
    {
      listoverlap123=listoverlap123+1
    } else if(list1[i]%in%list2&!list1[i]%in%list3) {
        listoverlap12=listoverlap12+1
    } else if(list1[i]%in%list3&!list1[i]%in%list2) {
        listoverlap13=listoverlap13+1
    } else list1unique=list1unique+1    
  }
  
  #Determine all the overlaps between 2 and 3
  for(i in 1:length(list2)) {
    if (list2[i]%in%list3&!list2[i]%in%list1) {
      listoverlap23=listoverlap23+1
    }
  }
  
  #Determine remaining overlaps
  list2unique=length(list2)-listoverlap12-listoverlap23-listoverlap123
  list3unique=length(list3)-listoverlap13-listoverlap23-listoverlap123

      
  return(c(list1unique,list2unique,list3unique,listoverlap12, listoverlap13, listoverlap23, listoverlap123))
}

#Example
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

#Create VennDiagram
pdf("Cutoffs245VennDiagram.pdf")
draw.pairwise.venn(overlaps[1], overlaps[2], overlaps[3], c("Cutoff=2", "Cutoff=4"), fill = c("blue", "red"))
dev.off()