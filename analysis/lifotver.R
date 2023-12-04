if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")

library(rtracklayer)
library(data.table)
weightData <- fread("~/Downloads/updated_bmi_PRS_weights.txt")

coords <- strsplit(weightData$id, split="[:-]")
weightData$chrom <- sapply(coords, function(x) x[1])
weightData$start <- as.numeric(sapply(coords, function(x) x[2]))
weightData$end <- weightData$start + 1 # end position is start + 1 for SNPs

chainFile <- import.chain("~/Downloads/hg19ToHg38.over.chain")
gr <- GRanges(seqnames = weightData$chrom,
              ranges = IRanges(start = weightData$start, end = weightData$end))
lifted <- liftOver(gr, chainFile)


weightDataLifted <- cbind(weightData, as.data.frame(lifted))
