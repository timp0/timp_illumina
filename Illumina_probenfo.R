genes<-read.delim("ref_genes.txt",stringsAsFactors=FALSE)

within=logical(384)
genename=character(384)
inside_genename=character(384)
inside_far=numeric(384)
far=numeric(384)

for (i in 1:384){

  chromy=paste("chr", data$probes$Chromosome[i], sep="")
  loc=data$probes$Start_loc[i]
  
  ingenes=genes[(genes$chrom==chromy),]
  
  between=(((ingenes$txStart<loc)&(ingenes$txEnd>loc))|((ingenes$txStart>loc)&(ingenes$txEnd<loc)))
  
  if (any(between)) {
    within[i]=T
    inside_genename[i]=ingenes$X.geneName[which(between)[1]]
    inside_far[i]=abs(loc-ingenes$txStart)
  } 
  
  genedist=abs(ingenes$txStart-loc)
  
  closest=ingenes[order(genedist)[1],]
  
  far[i]=genedist[order(genedist)[1]]

  genename[i]=closest$X.geneName
}


pdf("Movie/probe_island1.pdf", width=11,height=8)
island_dist=numeric(5)
island_dist[1]=sum(data$probes$UCSC_Dist_to_Island==0)
island_dist[2]=sum((data$probes$UCSC_Dist_to_Island>0)&(data$probes$UCSC_Dist_to_Island<501) )
island_dist[3]=sum((data$probes$UCSC_Dist_to_Island>499)&(data$probes$UCSC_Dist_to_Island<1001) )
island_dist[4]=sum((data$probes$UCSC_Dist_to_Island>999) & (data$probes$UCSC_Dist_to_Island<2001) )
island_dist[5]=sum((data$probes$UCSC_Dist_to_Island>1999))
barplot(island_dist)
dev.off()


