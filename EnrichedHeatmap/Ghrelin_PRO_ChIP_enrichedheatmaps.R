##### Import libraries
library(EnrichedHeatmap)
library(ComplexHeatmap) 
library(plyranges)
library(circlize)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(GenomicRanges)
library(magick)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

###### plyranes library define genes and TSS
hg38.genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
gene_c = mutate(anchor_center(hg38.genes), width = 30)
tss <- promoters(hg38.genes, upstream = 0, downstream = 1)

##### Set file paths for original files
# motiffile = 
# mufile = 
# proseqfile <- 
chipseqfile10 = "/Users/jessicawestfall/Documents/Dowell_lab/Experiment_data/Exp18_GO_offset/GHR_ChIPSeq_10_TRIM_HISAT2_HG38_pruned_PC.bedGraph"
chipseqfile11 = "/Users/jessicawestfall/Documents/Dowell_lab/Experiment_data/Exp18_GO_offset/GHR_ChIPSeq_11_TRIM_HISAT2_HG38_pruned_PC.bedGraph"
chipseqfile12 = "/Users/jessicawestfall/Documents/Dowell_lab/Experiment_data/Exp18_GO_offset/GHR_ChIPSeq_12_TRIM_HISAT2_HG38_pruned_PC.bedGraph"

##### Read in org files and convert to GRanges 
# read the PRO-seq bedgraph and split to positive and negative strand counts
#probg = read.table(proseqfile)
#gr_probg = GRanges(seqnames = probg[[1]], ranges = IRanges(probg[[2]], probg[[3]]), coverage = probg[[4]])
#pos_gr_probg = gr_probg %>% filter(coverage>0)
#neg_gr_probg = gr_probg %>% filter(coverage<0)

# motifs
#motif = read.table(motiffile)
#motif["loc"] = paste(motif$V1,":",  motif$V2,"-",motif$V3, sep="")
#gr_motif <- GRanges(seqnames = motif[[1]], ranges = IRanges(motif[[2]], motif[[3]]), score = motif[[5]], name= motif[[4]],  strand=motif[[6]], loc=motif[,"loc"])

# bidirectionals
#mu = read.table(mufile)
#gr_mu = GRanges(seqnames = mu[[1]], ranges = IRanges(mu[[2]], mu[[3]]), name= mu[[4]])

# chip-seq
chipDMSO = read.table(chipseqfile12)
gr_chip_DMSO = GRanges(seqnames = chipDMSO[[1]], 
                       ranges = IRanges(chipDMSO[[2]], chipDMSO[[3]]), 
                       coverage = chipDMSO[[4]])

chipBRD10min = read.table(chipseqfile10)
gr_chip_BRD10min = GRanges(seqnames = chipBRD10min[[1]], 
                           ranges = IRanges(chipBRD10min[[2]], chipBRD10min[[3]]), 
                           coverage = chipBRD10min[[4]])

chipBRD30min = read.table(chipseqfile11)
gr_chip_BRD30min = GRanges(seqnames = chipBRD30min[[1]], 
                           ranges = IRanges(chipBRD30min[[2]], chipBRD30min[[3]]), 
                           coverage = chipBRD30min[[4]])

##### Find mu(center) of all the files and expand the regions by 1500bp on both sides
# motifs
#gr_motif_c =  mutate(anchor_center(gr_motif), width = 1)
#gr_motif_c_plus1500 =  mutate(anchor_center(grmotif_c), width = 3000)

# bidir
#gr_mu_c <- mutate(anchor_center(gr_mu), width = 1)
#gr_mu_c_plus1500 <- mutate(anchor_center(gr_mu_c), width = 3000)

# ChIP peak
chipBRD10min_c <- mutate(anchor_center(gr_chip_BRD10min), width = 1)
gr_chipBRD10min_c_plus1500 <- mutate(anchor_center(chipBRD10min_c), width = 3000)

chipBRD30min_c <- mutate(anchor_center(gr_chip_BRD30min), width = 1)
gr_chipBRD30min_c_plus1500 <- mutate(anchor_center(chipBRD30min_c), width = 3000)

##### Find overlaps between motifs, bidir, chip
# bidir near ChIP, calculate distance
#mu_near_ChIP <- join_overlap_inner(gr_mu_c, gr_chip_c_plus1500)
#gr_mu_near_ChIP <- granges(mu_near_ChIP)
#dists_from_mu_near_ChIP = distanceToNearest(gr_mu_near_ChIP, gr_chip_c, select="arbitrary") #arbitrary returns nearest
#hist(mcols(dists_from_mu_near_ChIP)$distance, breaks=1500, main="Distance of closest bidir relative to ChIP peak center")

# ChIP peaks near motifs, return peaks
#chip_near_motifs_1500window <- join_overlap_inner(gr_motif_c, gr_chip_c_plus1500)
#gr_chip_near_motifs <- granges(chip_near_motifs_1500window)
#dists_from_mu_near_ChIPwithmotifs = distanceToNearest(gr_mu_near_ChIP, gr_chip_near_motifs)
#hist(mcols(dists_from_mu_near_ChIPwithmotifs)$distance, main="Distance of closest mu relative to ChIP peaks with motif")

##### Association between ChIP-seq peaks and targets (i.e. TSS) by normalizing into a matrix.  Enriched heatmap
matDMSO = normalizeToMatrix(gr_chip_DMSO, tss, value_column = "coverage",
                           extend = 5000, mean_mode = "w0", w = 50) #mean mode w0 accounts for all bases including 0 coverage
quantile(matDMSO, c(0, 0.25, 0.5, 0.75,  0.9, 0.99, 1))
col1_fun <- colorRamp2(quantile(matDMSO, c(0, 0.99)), c("white", "darkgreen"))
EnrichedHeatmap(matDMSO, col = col1_fun, row_km = 3, 
                name = "BDR4 DMSO", 
                column_title = "BDR4 DMSO coverage around TSS") 

mat_BRD4_10min = normalizeToMatrix(gr_chip_BRD10min, tss, value_column = "coverage",
                                   extend = 5000, mean_mode = "w0", w = 50) #mean mode w0 accounts for all bases including 0 coverage
quantile(mat_BRD4_10min, c(0, 0.25, 0.5, 0.75,  0.9, 0.99, 1))
col1_fun <- colorRamp2(quantile(mat_BRD4_10min, c(0, 0.99)), c("white", "darkgreen"))
EnrichedHeatmap(mat_BRD4_10min, col = col1_fun, row_km = 3, 
                name = "BRD4 10 mins", 
                column_title = "BRD4 10 min Ghrelin treatment coverage around TSS")

mat_BRD4_30min = normalizeToMatrix(gr_chip_BRD30min, tss, value_column = "coverage",
                                   extend = 5000, mean_mode = "w0", w = 50) #mean mode w0 accounts for all bases including 0 coverage
quantile(mat_BRD4_30min, c(0, 0.25, 0.5, 0.75,  0.9, 0.99, 1))
col1_fun <- colorRamp2(quantile(mat_BRD4_30min, c(0, 0.99)), c("white", "darkgreen"))
EnrichedHeatmap(mat_BRD4_30min, col = col1_fun, row_km = 3, 
                name = "BRD4 30 mins", 
                column_title = "BRD4 30 min treatment coverage around TSS")

ht_list = EnrichedHeatmap(matDMSO, col = col1_fun, row_km = 3, name = "BRD4 DMSO") +
  EnrichedHeatmap(mat_BDR4_10min, col = col1_fun, row_km = 3, name = "BRD4 10 mins") +
  EnrichedHeatmap(mat_BDR4_30min, col = col1_fun, row_km = 3, name = "BRD4 30 mins")
draw(ht_list)



