library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Drerio.UCSC.danRer10.refGene)
library(org.Dr.eg.db)
txdb <- TxDb.Drerio.UCSC.danRer10.refGene

#find overlaps between h2az and h3.3 peaks
overlaps <- findOverlaps(all_peaks.h2, all_peaks.h3, maxgap = 50)  # Allow slight gaps
common_peaks_custom <- pintersect(all_peaks.h2[queryHits(overlaps)], all_peaks.h3[subjectHits(overlaps)])
common_peaks_merged <- reduce(common_peaks_custom)

peakAnno.h2h3 <- annotatePeak(common_peaks_merged,
                              tssRegion = c(-3000, 3000),
                              TxDb = txdb,
                              annoDb = "org.Dr.eg.db")
plotDistToTSS(peakAnno.h2h3) + ggtitle("Distribution of H2A.Z+H3.3 binding loci realtive to TSS")
plotAnnoPie(peakAnno.h2h3)
plotAnnoPie(peakAnno.h2, col = c("lightblue","palegoldenrod","deepskyblue3", "green4", "indianred2", "darkorange1", "plum3", "thistle", "saddlebrown"))
H2H3_annotated <- as.data.frame(peakAnno.h2h3@anno)
write.csv(H2H3_annotated, file = "Common_peaks_h2azh33.csv")

H2_unique <- setdiff(all_peaks.h2, common_peaks_merged)
H3.unique <- setdiff(all_peaks.h3, common_peaks_merged)

peakAnno.h2_unique <- annotatePeak(H2_unique,
                                   tssRegion = c(-1000, 1000),
                                   TxDb = txdb,
                                   annoDb = "org.Dr.eg.db")
peakAnno.h3_unique <- annotatePeak(H3.unique,
                                   tssRegion = c(-1000, 1000),
                                   TxDb = txdb,
                                   annoDb = "org.Dr.eg.db")
H2_unique <- as.data.frame(peakAnno.h2_unique@anno)
H3.unique <- as.data.frame(peakAnno.h3_unique@anno)
write.csv(H2_unique, file = "Unique Peaks H2.csv")
write.csv(H3.unique, file = "unique Peaks H3.csv")
export(common_peaks_merged, ".../common_peaks_h2h3.bed", format = "bed")
