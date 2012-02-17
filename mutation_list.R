getcoord <- function(UCSC)  {
##Really stupid function to get the coordinates from the format from mskcc http://cbio.mskcc.org/CancerGenes/SelectAndDisplay.action

  #Split string using colon
  a=strsplit(as.character(UCSC), ":")
  #Get out chromosomes
  chrom=unlist(lapply(a, function(x) x[1]))

  #Take coordinates out
  c=unlist(lapply(a, function(x) x[2]))

  #Split them

  d=strsplit(c, " - ")

  #Split to begin and end
  begin=as.numeric(unlist(lapply(d, function(x) x[1])))
  finish=as.numeric(unlist(lapply(d, function(x) x[2])))

  coord=data.frame(chr=chrom, start=begin, end=finish)

  return(coord)
  

}

refcoord <- function(refgene)  {
#Using Refseq list get coordinates of transcripts

 
  

}
