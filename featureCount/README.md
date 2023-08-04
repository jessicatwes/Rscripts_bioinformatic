# FeatureCounts 
To perform downstream analysis such as differential gene expression (DEG) of your sequencing data, you need to make a count matrix that will consist of your read counts over your annotation file. The scripts in this directory include making your metadata table and piping that metadata table into your featureCounts script for it to form the matrix. 

## Scripts in this directory
1. metadata.R: This is the script that will create your metadata table. The columns are dependent on your data and what information would be valuable for you to distingish your samples.
2. featureCounts.R: This script uses the featureCount module from Rsubreads. You will need to provide it an annotation file to count over.
3. featureCount.sbatch: This is the bash script to run  metadata.R and featureCounts.R
