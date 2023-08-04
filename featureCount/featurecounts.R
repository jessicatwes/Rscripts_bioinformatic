##### R script for featureCounts count table ##### 
##### Set up the environment #####
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("Rsubread")
library(Rsubread)

##### Set up working directory #####
workdir <- '/scratch/Users/jewe1055/exp_files/exp99/'
setwd(workdir)
getwd()
outdir <- paste(workdir, 'output/featureCounts', '/', sep='') ##naming our outdir
dir.create(outdir, showWarnings = FALSE) ###creating the directory

pathdir <- read.delim('scripts/paths.txt', sep="\t")

### Pull in annotation file
hg38_refseq_gtf <- pathdir[pathdir$filename=="hg38_refseq_gtf", ]
hg38_refseq_5trunc_gtf <- pathdir[pathdir$filename == "hg38_refseq_5trunc_gtf", ]

### List of bam files. 
metadata <- read.csv('output/metadata.txt', header=TRUE, sep=",")

############################ SET UP PARAMETERS FOR ANNOTATION AND ASSAYS
### RNA-seq data for featureCounts
whichassay = 'RNA'
filetable <- subset(metadata, assay == whichassay )
filelist <-as.vector(filetable$bam_fp)
if(whichassay == 'RNA'){
  annot_file <- as.character(hg38_refseq_gtf$path)
  coverage <- featureCounts(files=filelist,
                            # annotation
                            annot.ext=annot_file,
                            isGTFAnnotationFile=TRUE,
                            GTF.featureType="exon",
                            GTF.attrType="gene_id",
                            
                            # level of summarization
                            useMetaFeatures=TRUE,
                            
                            # overlaps between reads and features
                            allowMultiOverlap=FALSE,
                            largestOverlap=TRUE,
                            
                            # multi-mapping reads
                            countMultiMappingReads=TRUE,
                            
                            #paired ends parameter
                            isPairedEnd=TRUE,
                            
                            # miscellaneous
                            nthreads=32)

  colnames(coverage$counts) <- filetable$label

  ### Write results
  fileroot<-paste0(outdir, whichassay, "_featureCounts")
  save.image(paste0(fileroot, ".RData"))

  write.csv(x=data.frame(coverage$annotation[,c("GeneID","Length")],
                           coverage$counts,stringsAsFactors=FALSE),
              paste0(fileroot, ".coverage.csv"),
              quote=FALSE, sep="\t",
              row.names=FALSE)
  write.csv(coverage$stat, paste(fileroot,".stat.csv", sep=""))
  write.csv(coverage$annotation, paste(fileroot,".annotation.csv", sep=""))
  write.csv(coverage$targets, paste(fileroot,".targets.csv", sep=""))
} 

### PRO-seq data for featureCounts
whichassay = 'PRO'
filetable <- subset(metadata, assay == whichassay )
filelist <-as.vector(filetable$bam_fp)
if(whichassay == 'PRO'){
  annot_file <- as.character(hg38_refseq_5trunc_gtf$path)
  coverage <- featureCounts(files=filelist,
                            # annotation
                            annot.ext=annot_file,
                            isGTFAnnotationFile=TRUE,
                            GTF.featureType="gene_length",
                            GTF.attrType="gene_id",
                            
                            # level of summarization
                            useMetaFeatures=TRUE,
                            
                            # overlaps between reads and features
                            allowMultiOverlap=FALSE,
                            largestOverlap=TRUE,
                            
                            # multi-mapping reads
                            countMultiMappingReads=TRUE,
                            
                            # paired ends parameter
                            isPairedEnd=FALSE,
                            
                            # strand-specific
                            strandSpecific=2,
                            
                            # miscellaneous
                            nthreads=32)

  colnames(coverage$counts) <- filetable$label

  ### Write results
  fileroot<-paste0(outdir, whichassay, "_featureCounts")
  save.image(paste0(fileroot, ".RData"))

  write.csv(x=data.frame(coverage$annotation[,c("GeneID","Length")],
                           coverage$counts,stringsAsFactors=FALSE),
              paste0(fileroot, ".coverage.csv"),
              quote=FALSE, sep="\t",
              row.names=FALSE)
  write.csv(coverage$stat, paste(fileroot,".stat.csv", sep=""))
  write.csv(coverage$annotation, paste(fileroot,".annotation.csv", sep=""))
  write.csv(coverage$targets, paste(fileroot,".targets.csv", sep=""))
}