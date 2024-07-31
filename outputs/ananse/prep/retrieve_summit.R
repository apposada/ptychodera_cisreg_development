peaks <- read.delim2("outputs/ananse/prep/peaks_normalized_min75.bed", header = FALSE)
peaks$V3 = peaks$V2+100
peaks$V2 = peaks$V2+100
write.table(peaks,file = "outputs/ananse/prep/peaks_normalized_min75_summit.bed", sep = "\t", quote = F, row.names = F, col.names = F)

