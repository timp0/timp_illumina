#Run using R CMD BATCH findseq1.R
#Finding chromosomal locations of sequences


#Load libraries in
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)

#Load the Fasta File
Il_probes <- read.DNAStringSet("Old_Illumina_probe1.fa", "fasta")


seqnames <- seqnames(Hsapiens)

#Initialize hit table
hits<-data.frame(row.names=c("seqname", "start", "end", "strand", "patternID"))


#Just the first 24 entries - these are the useful chromosomes
for (i in 1:24) {
  seqname <- seqnames[i]
  chromo <- Hsapiens[[seqname]]

  cat(">>> Finding all hits in chromosome", seqname, "...\n")


  for (j in seq_len(length(Il_probes))) {   
    probeID <- names(Il_probes)[j]
    cat(">> Finding probe", probeID, "...\n")
    seqy <- Il_probes[[j]]
    plus_matches <- matchPattern(seqy, chromo)
    #Make data set for forward strand match
    plus_hit <- data.frame(probeID = rep.int(probeID, length(plus_matches)),
                           seqname=rep.int(seqname, length(plus_matches)),
                           start=start(plus_matches),
                           end=end(plus_matches),
                           strand=rep.int("+", length(plus_matches)))

    rev_seqy <- reverseComplement(seqy)
    minus_matches <- matchPattern(rev_seqy, chromo)
    minus_hit <- data.frame(probeID = rep.int(probeID, length(minus_matches)),
                            seqname=rep.int(seqname, length(minus_matches)),
                            start=start(minus_matches),
                            end=end(minus_matches),
                            strand=rep.int("-",length(minus_matches)))

    hits <- rbind(hits, plus_hit, minus_hit)
    }
}


write.table(hits, file="try.csv", append=, quote=FALSE, sep=",")






                        
