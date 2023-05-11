args<-commandArgs(TRUE)

# intersecrt between CAGE ckusters and ATAC-seq
peaks <- cage.table(args[1], header=FALSE)

colnames(peaks)[1] <- "cage.chr"
colnames(peaks)[2] <- "cage.start"
colnames(peaks)[3] <- "cage.end"
colnames(peaks)[5] <- "cage.count"
colnames(peaks)[6] <- "cage.strand"
colnames(peaks)[7] <- "cluster.chr"
colnames(peaks)[8] <- "cluster.start"
colnames(peaks)[9] <- "cluster.end"

# get maximum CAGE peak per intersect
max.peaks <-aggregate(peaks$cage.count, list(peaks$cluster.chr, peaks$cluster.start, peaks$cluster.end, peaks$cage.strand), FUN=max, na.rm=TRUE)
colnames(max.peaks)[1] <- "cluster.chr"
colnames(max.peaks)[2] <- "cluster.start"
colnames(max.peaks)[3] <- "cluster.end"
colnames(max.peaks)[4] <- "cage.strand"
colnames(max.peaks)[5] <- "cage.count"

# select maximum from the intersect
max.peaks.positions <- merge(max.peaks, peaks, by=c("cluster.chr","cluster.start","cluster.end","cage.strand","cage.count"), all.x=TRUE)
max.peaks.positions <- max.peaks.positions[c("cage.chr", "cage.start","cage.end","cage.count","V4","cage.strand")]

# save maximum CAGE clusters 
write.table(max.peaks.positions, args[2], row.names=FALSE, sep="\t", quote=FALSE)




