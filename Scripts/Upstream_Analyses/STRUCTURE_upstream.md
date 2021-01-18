STRUCTURE\_upstream
================

Code for STRUCTURE analyses run on an HPC cluster (Amarel at Rutgers) before transferring output files to local computer for downstream analyses and visualization in R.

Downloaded STRUCTURE v.2.3.4 to Amarel and unzipped/installed in `~/Programs`.

*Note: Process below is for running STRUCTURE with the full mac2 dataset. For the other datasets, inputs and file names were changed as needed.*

Created STRUCTURE input file with PGDSpider v.2.1.1.5.

-   Took `output.hicov2.snps.only.mac2.vcf` file and converted to `STRUCTURE_mac2.txt`.

Used following code to run PGDSpider v.2.1.1.5 on Amarel:

``` bash
#grabbed an interactive node with following line of code
srun --partition=main --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=2000 --time=00:30:00 --export=ALL --pty bash -i

#called #PGDSpider in ~/Programs/PGDSpider_2.1.1.5
java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar -inputfile output.hicov2.snps.only.mac2.vsc -inputformat VCF -outputfile STRUCTURE_mac2.txt -outputformat STRUCTURE -spid VCF_STRUCTURE.spid

#NOTE: if don't have spid file made, or know proper format, running the previous line of code without the -spid argument will generate a template spid file to modify as needed
```

Spid format for PGDSpider:

``` bash
# VCF Parser questions
PARSER_FORMAT=VCF

# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION= Pop_assignments.txt #txt file with indv names in one column and pop assignment (Pop_1, etc.) in next (no headers, tab-delimited)
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=  DIPLOID
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION= TRUE
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION= TRUE
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION= FALSE
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION= FALSE

# STRUCTURE Writer questions
WRITER_FORMAT=STRUCTURE

# Specify the locus/locus combination you want to write to the STRUCTURE file:
STRUCTURE_WRITER_LOCUS_COMBINATION_QUESTION=
# Do you want to include inter-marker distances?
STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=FALSE
# Specify which data type should be included in the STRUCTURE file  (STRUCTURE can only analyze one data type per file):
STRUCTURE_WRITER_DATA_TYPE_QUESTION=SNP
# Save more specific fastSTRUCTURE format?
STRUCTURE_WRITER_FAST_FORMAT_QUESTION=FALSE
```

Moved `STRUCTURE_mac2.txt` to `~/Programs/console`.

Created `STRUCTURE.sbatch` file to run STRUCTURE:

``` bash
#!/bin/bash

#SBATCH --partition=p_mlp195 #p_mlp195 for pinsky lab partition, otherwise use main (or EOAS or E&E)
#SBATCH --job-name=STRUCURE
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

#script for running STRUCTURE on Amarel

structure -m mainparams_K1 -o mac2_K1_run1
```

`mainparams_K1` file structure below. Ran with 100,000 burn-in iterations & 10,000 MCMC reps after burn-in:

``` bash
KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
FILE extraparams.


"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean
        (ie insert 1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes!)


Basic Program Parameters

#define MAXPOPS    1      // (int) number of populations assumed
#define BURNIN    100000   // (int) length of burnin period
#define NUMREPS   10000   // (int) number of MCMC reps after burnin

Input/Output files

#define INFILE   STRUCTURE_mac2.txt   // (str) name of input data file
#define OUTFILE  mac2_K1_run1.out  //(str) name of output data file

Data file format

#define NUMINDS    25    // (int) number of diploid individuals in data file
#define NUMLOCI    4212    // (int) number of loci in data file
#define PLOIDY       2    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line


#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   1     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says
                              whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier

#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data
                             before the genotype data start.

#define MARKERNAMES      1  // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                            // and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances
                            // between loci


Advanced data file options

#define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line
                                  indicating phase (linkage model)
#define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -999 // (int) for use in some analyses of polyploid data
```

`extraparams` file structure (top lines). Ran with the admixture model, correlated allele frequencies and incorporated prior location information:

``` bash
EXTRA PARAMS FOR THE PROGRAM structure.  THESE PARAMETERS CONTROL HOW THE
PROGRAM RUNS.  ATTRIBUTES OF THE DATAFILE AS WELL AS K AND RUNLENGTH ARE
SPECIFIED IN mainparams.

"(int)" means that this takes an integer value.
"(d)"   means that this is a double (ie, a Real number such as 3.14).
"(B)"   means that this variable is Boolean
        (ie insert 1 for True, and 0 for False).

PROGRAM OPTIONS

#define NOADMIX     0 // (B) Use no admixture model (0=admixture model, 1=no-admix)
#define LINKAGE     0 // (B) Use the linkage model model
#define USEPOPINFO  0 // (B) Use prior population information to pre-assign individuals to clusters

#define LOCPRIOR    1 //(B)  Use location information to improve weak data (LOCISPOP = 1)

#define FREQSCORR   1 // (B) allele frequencies are correlated among pops

#everything else left unchanged
```

Submitted job:

``` bash
sbatch STRUCTURE.sbatch
```

Ran 5 replicates of each value of K (1-5). Done for full mac2 dataset, only SNPs in HWE, only non-outlier SNPS, and only outlier SNPs.
