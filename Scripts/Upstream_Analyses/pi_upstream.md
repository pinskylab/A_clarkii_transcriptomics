pi\_upstream
================

Code for calculating per-site pi using VCFtools. Run on an HPC cluster (Amarel at Rutgers) before transferring output files to local computer for downstream analyses and visualization in R.

VCFtools v.1.16 was already installed as a module on Amarel. To load in working environment, use following code:

``` bash
module load VCFtools
```

First, created sub-setted VCF files for each sampling site (code for Japan below):

``` bash
#created J_individuals.txt (column of individuals in Japan)

#used VCFtools to filter output.hicov2.snps.only.mac2.vcf to only include individuals from Japan
vcftools --vcf output.hicov2.snps.only.mac2.vcf --keep J_individuals.txt --recode --recode-INFO-all --out Jmac2

#repeated for Philippines & Indonesia
```

Calculated site pi, both within each population and with the full dataset (all individuals included in VCF):

``` bash
#NOTE: don't grab interactive node for VCFtools jobs bc run in seconds, BUT if is more data-intensive, should get off log-in node to do this.

vcftools --vcf Jmac2.vcf --site-pi --out Jmac2

#repeated for Philippines, Indonesia, and the original vcf file w/all individuals
```

Copied `*sites.pi` files to local computer and opened in Excel to create .csv files that were then read into R for downstream analyses (`Scripts/pi.R`).
