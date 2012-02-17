#Run using R CMD BATCH findseq1.R
#Finding chromosomal locations of sequences


#Load libraries in
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

#Load the Fasta File
ill <- read.csv("timp_ordered1.csv",stringsAsFactors=FALSE)

seqnames <- seqnames(Hsapiens)


#Initialize hit table
hits<-data.frame(row.names=c("seqname", "start", "end", "strand", "patternID"))


#Just the first 24 entries - these are the useful chromosomes
for (i in 1:24) {
  seqname <- seqnames[i]
  chromo <- Hsapiens[[seqname]]

  cat(">>> Finding all hits in chromosome", seqname, "...\n")


  for (j in 1:dim(ill)[1] ) {   
    probeID <- ill[j,11]
    cat(">> Finding probe", probeID, "...\n")
    seqy <- strsplit(ill[j,2], ']')[[1]][2]
    seqy=DNAString(paste(seqy, "CG", sep=""))
    plus_matches <- matchPattern(seqy, chromo)
    #Make data set for forward strand match
    plus_hit <- data.frame(probeID = rep.int(probeID, length(plus_matches)),
                           seqname=rep.int(seqname, length(plus_matches)),
                           start=start(plus_matches),
                           end=end(plus_matches),
                           Cloc=start(plus_matches),
                           strand=rep.int("+", length(plus_matches)))

    rev_seqy <- reverseComplement(seqy)
    minus_matches <- matchPattern(rev_seqy, chromo)
    minus_hit <- data.frame(probeID = rep.int(probeID, length(minus_matches)),
                            seqname=rep.int(seqname, length(minus_matches)),
                            start=start(minus_matches),
                            end=end(minus_matches),
                            Cloc=end(minus_matches),
                            strand=rep.int("-",length(minus_matches)))

    hits <- rbind(hits, plus_hit, minus_hit)
    }
}


write.table(hits, file="try.csv", append=, quote=FALSE, sep=",")






                        
