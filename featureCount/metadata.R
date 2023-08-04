library("dplyr")

pathdir <- read.delim('/scratch/Users/jewe1055/exp_files/exp99/scripts/paths.txt', sep="\t")

t21color <- c("#fec44f", "#d95f0e")
d21color <- c("#ece7f2", "#2b8cbe")
short_treatment <- c("BSA", "IFN")

##### RNA-seq assay
rnabam_dir <- pathdir[pathdir$filename == "RNA_bam_fp", ]
rna_bam <- c("RNA-BSA-Dave", "RNA-IFN-Dave", "RNA-BSA-Ethan", "RNA-IFN-Ethan", "RNA-BSA-Eric", "RNA-IFN-Eric", "RNA-BSA-ChenChao", "RNA-IFN-ChenChao",
         "RNA-BSA-Khaondo", "RNA-IFN-Khaondo", "RNA-BSA-Niyilolawa", "RNA-IFN-Niyilolawa", "RNA-BSA-Pedro", "RNA-IFN-Pedro",
         "RNA-BSA-Sengbe", "RNA-IFN-Sengbe", "RNA-BSA-Srivathani", "RNA-IFN-Srivathani", "RNA-BSA-Ursula", "RNA-IFN-Ursula")
rna_name <- c("Dave", "Dave", "Ethan", "Ethan", "Eric", "Eric", "ChenChao", "ChenChao", "Khaondo", "Khaondo",
              "Niyilolawa", "Niyilolawa", "Pedro", "Pedro", "Sengbe", "Sengbe", "Srivathani","Srivathani", 
              "Ursula",  "Ursula")
rna_genotype <- c(rep("downsyn", 4), rep("typical", 16))
rna_treatment <- rep(short_treatment, 10)
rna_color <- c(rep(t21color, 2), rep(d21color, 8))
rna_assay <- rep("RNA", 20)
rna_metadata <- data.frame(rna_bam, rna_name, rna_genotype, rna_treatment, rna_color, rna_assay)
rna_metadata <- rbind(rna_metadata, rna_metadata, rna_metadata)
rna_metadata$rna_replicates <- c(rep(1, 20), rep(2, 20), rep(3, 20))
rna_metadata$rna_bamdir <- rnabam_dir$path
rna_metadata$rna_bam_fp <- paste(rna_metadata$rna_bamdir,rna_metadata$rna_bam, '-', rna_metadata$rna_replicates, '.sorted.bam', sep = "")
rna_metadata$rna_group <- paste(rna_metadata$rna_genotype, rna_metadata$rna_treatment, sep='_')
rna_metadata$rna_label <- paste(rna_metadata$rna_name, rna_metadata$rna_treatment, rna_metadata$rna_replicates, rna_metadata$rna_assay, sep='__')
rna_metadata <- rna_metadata %>% 
  dplyr::rename("bam" = "rna_bam", "name" = "rna_name", "genotype" = "rna_genotype",
         "treatment" = "rna_treatment", "color" = "rna_color", "assay" = "rna_assay",
         "replicates" = "rna_replicates", "bamdir" = "rna_bamdir", "bam_fp" = "rna_bam_fp",
         "group" = "rna_group", "label" = "rna_label")

##### PRO-seq assay
probam_dir <- pathdir[pathdir$filename == "PRO_bam_fp", ]
pro_bam <- c("PRO-BSA-Dave", "PRO-IFNB-Dave", "PRO-BSA-Ethan", "PRO-IFNB-Ethan", "PRO-BSA-Eric", "PRO-IFNB-Eric", 
             "PRO-BSA-ChenChao", "PRO-IFNB-ChenChao", "PRO-BSA-Khaondo", "PRO-IFNB-Khaondo", "PRO-BSA-Niyilolawa", "PRO-IFNB-Niyilolawa",
             "PRO-BSA-Pedro", "PRO-IFNB-Pedro", "PRO-BSA-Srivathani", "PRO-IFNB-Srivathani")
pro_name <- c("Dave", "Dave", "Ethan", "Ethan", "Eric", "Eric", "ChenChao", "ChenChao", "Khaondo", "Khaondo",
              "Niyilolawa", "Niyilolawa", "Pedro", "Pedro", "Srivathani","Srivathani")
pro_genotype <- c(rep("downsyn", 4), rep("typical", 12))
pro_treatment <- rep(short_treatment, 8)
pro_color <- c(rep(t21color, 2), rep(d21color, 6))
pro_assay <- rep("PRO", 16)
pro_metadata <- data.frame(pro_bam, pro_name, pro_genotype, pro_treatment, pro_color, pro_assay)
pro_metadata <- rbind(pro_metadata, pro_metadata)
pro_metadata$pro_replicates <- c(rep(1, 16), rep(3, 6), rep(2, 10))
pro_metadata$pro_bamdir <- probam_dir$path
pro_metadata$pro_bam_fp <- paste(pro_metadata$pro_bamdir, pro_metadata$pro_bam, '-', pro_metadata$pro_replicates, '.sorted.bam', sep = "")
pro_metadata$pro_group <- paste(pro_metadata$pro_genotype, pro_metadata$pro_treatment, sep='_')
pro_metadata$pro_label <- paste(pro_metadata$pro_name, pro_metadata$pro_treatment, pro_metadata$pro_replicates, pro_metadata$pro_assay, sep='__')
pro_metadata <- pro_metadata %>% 
  dplyr::rename("bam" = "pro_bam", "name" = "pro_name", "genotype" = "pro_genotype",
                "treatment" = "pro_treatment", "color" = "pro_color", "assay" = "pro_assay",
                "replicates" = "pro_replicates", "bamdir" = "pro_bamdir", "bam_fp" = "pro_bam_fp",
                "group" = "pro_group","label" = "pro_label")

##### Master metadata table
metadata <- rbind(rna_metadata, pro_metadata)
write.csv(metadata, '/scratch/Users/jewe1055/exp_files/exp99/output/metadata.txt')

