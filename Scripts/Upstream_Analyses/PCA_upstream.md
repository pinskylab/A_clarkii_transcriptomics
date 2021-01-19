PCA\_upstream
================

Code for generating PCAs (eigenvalues & eigenvectors) using plink. Run on an HPC cluster (Amarel at Rutgers) before transferring output files to local computer for downstream analyses and visualization in R.

Downloaded plink v.1.9 to Amarel and installed in `~/Programs`.

*Note: Process below is for creating eigenvalues/vectors with the full mac2 dataset. For the other datasets, inputs and file names were changed as needed.*

Created `PCA.sbatch` file to run plink for PCAs:

``` bash

#!/bin/bash

#SBATCH --partition=p_mlp195 #p_mlp195 for pinsky lab partition, otherwise use main (or EOAS or E&E)
#SBATCH --job-name=PCA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=2:00:00

#SBATCH --requeue
#SBATCH --export=ALL
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=rene.clark@rutgers.edu

#script for running Plink1.9 for PCA on Amarel
#NOTE: probably don't need to sbatch, runs very quickly (either on log-in node or interactive node)

cd ~/Aclarkii_Transcriptome/PCA

plink --pca --allow-extra-chr --vcf output.hicov2.snps.only.mac2.vcf --out output.hicov2.snps.only.mac2
```

Plink creates 4 files:

1.  `*.eigenvec`: Contains eigenvectors
2.  `*.eigenval`: Contains eigenvalues
3.  `*.log`: Contains information on \# variants and individuals used in analysis
4.  `*.nosex`: Not relevant

Repeated process for full mac2 dataset, only SNPs in HWE, only non-outlier SNPS, and only outlier SNPs.

Copied `*.eigenvec` & `*.eigenval` files to local computer and opened in Excel to create .csv files that were then read into R for downstream analyses (`Scripts/PCAs.R`).
