library("dplyr")

t21color <- c("#fec44f", "#d95f0e")
d21color <- c("#ece7f2", "#2b8cbe")
short_treatment <- c("BSA", "IFN")

##### RNA-seq assay
rnabam_dir <- data/
rna_bam <- c("RNA-BSA-Ethan", "RNA-IFN-Ethan", "RNA-BSA-Eric", "RNA-IFN-Eric")
rna_name <- c("Ethan", "Ethan", "Eric", "Eric")
rna_genotype <- c(rep("downsyn", 2), rep("typical", 2))
rna_treatment <- rep(short_treatment, 4)
rna_color <- c(rep(t21color, 1), rep(d21color, 1))
rna_assay <- rep("RNA", 4)
rna_metadata <- data.frame(rna_bam, rna_name, rna_genotype, rna_treatment, rna_color, rna_assay)
rna_metadata <- rbind(rna_metadata, rna_metadata, rna_metadata)
rna_metadata$rna_replicates <- c(rep(1, 8), rep(2, 8), rep(3, 8))
rna_metadata$rna_bamdir <- rnabam_dir$path
rna_metadata$rna_bam_fp <- paste(rna_metadata$rna_bamdir,rna_metadata$rna_bam, '-', rna_metadata$rna_replicates, '.sorted.bam', sep = "")
rna_metadata$rna_group <- paste(rna_metadata$rna_genotype, rna_metadata$rna_treatment, sep='_')
rna_metadata$rna_label <- paste(rna_metadata$rna_name, rna_metadata$rna_treatment, rna_metadata$rna_replicates, rna_metadata$rna_assay, sep='__')
rna_metadata <- rna_metadata %>% 
  dplyr::rename("bam" = "rna_bam", "name" = "rna_name", "genotype" = "rna_genotype",
         "treatment" = "rna_treatment", "color" = "rna_color", "assay" = "rna_assay",
         "replicates" = "rna_replicates", "bamdir" = "rna_bamdir", "bam_fp" = "rna_bam_fp",
         "group" = "rna_group", "label" = "rna_label")

##### Write metadata table
write.csv(metadata, 'metadata.txt')

