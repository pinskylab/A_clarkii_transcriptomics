Fastsimcoal2\_upstream
================

Code for *fastsimcoal2* analyses run on an HPC cluster (Amarel at Rutgers) before transferring output files to local computer for downstream analyses.

***Note:** A lot of this code was first written by Jennifer Hoey, and modified for this project. For more information, check out her repo: <https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/>. *

For *fastsimcoal2*, had to first generate the observed multidimensional site frequency spectrum (MSFS) with `easySFS.py` from (<https://github.com/isaacovercast/easySFS>). To main neutrality and retain rare variants, built the observed MSFS with a dataset unfiltered for minor allele counts but without outlier SNPs (5,662 SNPs). All individuals were included.

**To generate MSFS**

Opened up conda on Amarel and created environment for easySFS:

``` bash
module use /projects/community/modulefiles #to access community modules not in core set
module load py/data-science/stack/5.1.0-kp807 #module with anaconda

conda create -n easySFS #create easySFS environment
conda config --add channels bioconda #channel with bioinformatics packages
conda install -c bioconda dadi pandas #add dadi & pandas packages

#downloaded easySFS repo fron GitHub and unzipped in ~/Programs directory
unzip easySFS-master.zip
mv easySFS-master easySFS
cd easySFS
chmod 777 easySFS.py #to make executable
```

Copied appropriate VCF file (`output.hicov2.snps.only.mac1.nooutliersfinal.reordered.vcf`) to `~/Programs/easySFS` directory. This VCF was made by filtering original vcf `output.hicov2.snps.only.vcf` to keep only SNPs that were polymorphic (had at least 1 individual with at least 1 minor allele) and were not outlier SNPs. The individuals were also re-ordered to match the geographic order (Japan, Philippines, Indonesia). Filtering was done with VCFtools.

Previewed SFS with following code:

``` bash
conda activate easySFS

./easySFS.py -i output.hicov2.snps.only.mac1.nooutliersfinal.reordered.vcf -p popfile.txt --preview -a 

#popfile.txt is a 2 column, tab-separated file with individuals in first column and pop name in second. 
#Sometimes need to recreate this as doesn't always properly read all populations.

#NOTE: -a flag forces easySFS to keep all SNPs within each RAD locus (doesn't randomly sample 1 SNP per locus)
```

Preview shows \# segregating sites depending on \# individuals sampled in each pop.

-   Because, if missing sites, downsampling can increase \# of segregating sites. Since we don't have missing data, used full dataset (all individuals).

Created SFS with following code:

``` bash
./easySFS.py -i output.hicov2.snps.only.mac1.nooutliersfinal.reordered.vcf -p popfile.txt --proj 16,20,14 -a -o output -f
```

Above code creastes all possible SF (each pop, 2D-joint SFS and MSFS) formatted for input into *dadi* and *fastsimcoal2*. Used the MSFS for *fastsimcoal2* (`output_MSFS.obs` renamed to `apcl_MSFS.obs`).

**To run *fastsimcoal2***

First, specified the model with the `*.tpl` and `*.est` files. The `*.tpl` file contains the model structure (sampling scheme information, number/types of demographic events, etc.). The `*.est` file contains the prior distributions for each parameter estimated in the model.

-   **NOTE:** It's important to keep naming consistent! All of these files (along with the `*.obs` file) should have the same prefix, so *fastsimcoal2* knows what to read in.

The `apcl.tpl` file (also included in the `Data` directory):

``` bash
#Any place where there is a word is a parameter that fsc is going to estimate
#Any place where there is a number is a set parameter
#fsc works BACKWARDS in time (coalescent model). So the first demographic event happened most recently in time, then moves backwards from there

//Amphiprion connectivity: Japan, Philippines, Indonesia
3 samples
//Population effective sizes (=2N, number of alleles)
POPONE
POPTWO
POPTHREE
//Sample sizes: # samples (=2N, number of alleles)
16
20
14
//Growth rates
0
0
0
//Number of migration matrices : If 0 : No migration between demes
3
//migration matrix 0: sources in rows, sinks in columns (backwards in time)
0 DISPONE 0
DISPONE 0 DISPTWO
0 DISPTWO 0
//migration matrix 1: At 3600 gen (when J & P merge to P), migration only P-I
0 0 0
0 0 DISPTWO
0 DISPTWO 0
//migration matrix 2: At TDIVTWO (when P & I merge to I), no migration
0 0 0
0 0 0
0 0 0
//Historical event format: time, src, sink, % mig, new Ne scaling factor for sink, new r, new mig matrix:
2 historical events
3600 0 1 1 1 0 1
TDIVTWO 1 2 1 1 0 2
//Number of chromosomes
1 0
//Number of linkage blocks on Chromosome 1
1
//Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters
FREQ 1 0 1E-8 OUTEXP
```

The `apcl.est` file (also included in the `Data` directory):

``` bash
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt?       #name   #dist.  #min    #max
//all N are in number of haploid individuals
1 POPONE logunif 100 100000 output
1 POPTWO logunif 100 100000 output
1 POPTHREE logunif 100 100000 output
0 DISPONE unif 0 0.3 output
0 DISPTWO unif 0 0.3 output
1 TDIVTWO unif 3600 6000 output

[RULES]

[COMPLEX PARAMETERS]
```

Once the `.tpl`, `.est`, and MSFS files were created, downloaded *fastsimcoal26* to Amarel and unzipped/installed in `~/Programs`.

``` bash
wget http://cmpg.unibe.ch/software/fastsimcoal2/downloads/fsc26_linux64.zip
unzip fsc26_linux64.zip
cd fsc26_linux64
chmod 777 fsc26 #make executable

#copied .tpl, .est, and .obs files to this fsc26_linux64 as well
```

Created a directory in `/scratch` to read, write, and execute fsc jobs from (should go faster and don't have to worry about storage space in home directory).

Created `run_fscmodel_singlethread.sbatch` script for *fastsimcoal2* to perform 100,000 coalescent simulations to estimate the SFS, with 40 optimization (ECM) cycles. This script was modified from <https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/run_model6_singlethread.sh>.

-   This script sets up file structure that will be important for gleaning information from the replicate runs and easily grabbing maximum likelihoods and parameter estimates. Briefly, it makes a `slurm-out` directory that contains the output from each replicate run. In this directory, are the same number of directories as replicate runs, each named by the slurm job id the run was submitted under. The script then copies all input files into each `slurm-out/SLURM_JOBID` directory and runs *fastsimcoal2* from there each time.
    -   It also writes several settings files that help with troubleshooting.

``` bash
#!/bin/bash
#SBATCH --job-name=fsc
#SBATCH --partition=p_mlp195,main #p_mlp195 for pinsky lab partition, otherwise use main (or EOAS or E&E)
#SBATCH -N1
#SBATCH -n1                        #28 cores on amarel
#SBATCH --cpus-per-task=1         #28 cores on amarel
############################
#SBATCH --mem=5000               #128GB /node of amarel
#SBATCH --time=10:00:00            #max time is 3 days:  3-00:00:00  or  72:00:00
#SBATCH --export=ALL
#SBATCH --requeue       #optional: this prevent the job from restarting if the node fails or the jobs is preempted.
#SBATCH -o /scratch/rdc129/Aclarkii_Demo/slurm-%j.out   ## the default file name is "slurm-%j.out"
#SBATCH -e /scratch/rdc129/Aclarkii_Demo/slurm-%j.err
#SBATCH --mail-type=END                # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=rene.clark@rutgers.edu # Email to which notifications will be sent

#############################################################
# standard output is saved in a file:  slurm-$SLURM_JOBID.out
#############################################################

mkdir -p /scratch/$USER/Aclarkii_Demo/slurm-out/$SLURM_JOBID
cp /home/$USER/Programs/fsc26_linux64/fsc26 /scratch/$USER/Aclarkii_Demo/apcl.tpl /scratch/$USER/Aclarkii_Demo/apcl.est /scratch/$USER/Aclarkii_Demo/apcl_MSFS.obs /scratch/$USER/Aclarkii_Demo/slurm-out/$SLURM_JOBID
cd /scratch/$USER/Aclarkii_Demo/slurm-out/$SLURM_JOBID

#create a personal dirctory on the node scratch disk if needed.
#mkdir -p /mnt/scratch/$USER/$SLURM_JOBID

#this will tell you when the job started and the host.
date=`date "+%Y.%m.%d-%H.%M.%S"`
hostname=`hostname`

#to print the variable -- echo
echo $date
echo $hostname

#obtain the environment variables
#this is useful for reference and troubleshooting issues.
env >               /scratch/$USER/Aclarkii_Demo/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-all.out
env | grep SLURM >  /scratch/$USER/Aclarkii_Demo/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-SLURM.out

#start time
date

#start the simulation
srun /home/rdc129/Programs/fsc26_linux64/fsc26 -t apcl.tpl -e apcl.est -n 100000 -m -u -M -L 40 -0 --numBatches 1 --cores 1 > /scratch/$USER/Aclarkii_Demo/slurm-out/$SLURM_JOBID/$SLURM_JOBID-fastsimcoal26.$SLURM_NNODES.$SLURM_NODELIST.$date.a.txt

#end time
date
```

Once the model was running well, wrote a shell script (`run_sbatch.sb`) to submit the previous script many times in a job array (need multiple replicate runs to pick best parameter estimates from). This script was modified from <https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/run_sbatch_Ne.sh>.

``` bash
#submits previous script 50x

#!/bin/bash
#SBATCH --partition=p_mlp195,main

for i in {1..50}
do
sbatch run_fscmodel_singlethread.sbatch
sleep 1

done
```

Once all 50 replicate runs finished running, wrote `cat_ML.sbatch` script to concatenate the results (`.bestlhoods` files) into a single file.

-   This script was modified from <https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/cat_bestlhoods>.

``` bash
#!/bin/bash

for file in /scratch/$USER/Aclarkii_Demo/slurm-out/*/apcl/apcl.bestlhoods
do
sed '1d' "$file" > "$file.trimmed"

cat /scratch/$USER/Aclarkii_Demo/slurm-out/*/apcl/apcl.bestlhoods.trimmed > /scratch/$USER/Aclarkii_Demo/apcl_raw.bestlhoods

done
```

Copied `apcl_raw.bestlhoods` file to local computer and read into Excel to sort by maximum likelihood (ML) and find the parameters from the best ML run. These are the best point estimates.

**To create 95% CIs**

To create the 95% CIs, first had to create 100 bootstrapped SFS from the raw SFS to run in *fastsimcoal2*. Used the `Bootstrap_forSFS.R` script to bootstrap `output.hicov2.snps.only.mac1.nooutliersfinal.reordered.vcf` 100X. Then used easySFS to create MSFS from each of the bootstrapped VCFs (with the same code as was used to create the original MSFS).

Moved bootstrapped MSFS to *fastsimcoal2* `/scratch` directory. Each MSFS was placed in its own directory, all of which were placed in a `bootstrapSFS` directory (ex. the first bootstrapped MSFS was placed in `/scratch/rdc129/Aclarkii_Demo/bootstrapSFS/apcl_boot1`).

Wrote `run_fscmodel_CI_singlethread.sbatch` script to run *fastsimcoal2* for each of the bootstrapped directories. This script was modified from <https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/run_fsc_Ne_CI_singlethread.sh>.

-   The structure of this script is very similar to the one for the original MSFS (creates `slurm-out/SLURM_JOBID` directories in each bootstrapped directory). One major difference is this script uses the parameter values from the best-fit model with the observed dataset as the starting values for each replicate run (to speed things up). This is the `apcl.pv` file from the best-fit model, and should be copied to the starting directory (here, `scratch/rdc129/Aclarkii_Demo`) before running the script.
-   To run the script, need `dir_names.txt` file, which is a list of all the bootstrapped SFS directory names (apcl\_boot1, apcl\_boot2, apcl\_boot3, etc.).

``` bash
#!/bin/bash
#SBATCH --job-name=fsc_boot
#SBATCH --partition=p_mlp195,p_eoas_1,p_deenr_1,main #p_mlp195 for pinsky lab partition, p_eoas_1 for EOAS partition, p_deenr_1 for E&E partition
#SBATCH -N1
#SBATCH -n1                        #28 cores on amarel
#SBATCH --cpus-per-task=1         #28 cores on amarel
############################
#SBATCH --mem=5000               #128GB /node of amarel
#SBATCH --time=33:00:00            #max time is 3 days:  3-00:00:00  or  72:00:00
#SBATCH --export=ALL
#SBATCH --array=1-100
#SBATCH --requeue       #optional: this prevent the job from restarting if the node fails or the jobs is preempted.
#SBATCH -o /scratch/rdc129/Aclarkii_Demo/slurm-%A_%a.out   #%A is the jobid and %a is the array index
#SBATCH -e /scratch/rdc129/Aclarkii_Demo/slurm-%A_%a.err   #/dev/null
#SBATCH --mail-type=END                #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=rene.clark@rutgers.edu #Email to which notifications will be sent

#############################################################
# standard output is saved in a file:  slurm-$SLURM_JOBID.out
#############################################################

echo "Starting task $SLURM_ARRAY_TASK_ID"
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" dir_names.txt)
cd $DIR

mkdir -p /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/slurm-out/$SLURM_JOBID
cp /home/$USER/Programs/fsc26_linux64/fsc26 /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/slurm-out/$SLURM_JOBID
cp /scratch/$USER/Aclarkii_Demo/apcl.tpl /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/slurm-out/$SLURM_JOBID/apcl_boot.tpl
cp /scratch/$USER/Aclarkii_Demo/apcl.est /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/slurm-out/$SLURM_JOBID/apcl_boot.est
cp /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/apcl_boot_MSFS.obs /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/slurm-out/$SLURM_JOBID
cp /scratch/$USER/Aclarkii_Demo/apcl.pv /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/slurm-out/$SLURM_JOBID/apcl_boot.pv
cd /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/slurm-out/$SLURM_JOBID

#Create a personal dirctory on the node scratch disk if needed.
#mkdir -p /mnt/scratch/$USER/$SLURM_JOBID


#this will tell you when the job started and the host.
date=`date "+%Y.%m.%d-%H.%M.%S"`
hostname=`hostname`

#to print the variable -- echo
echo $date
echo $hostname

#obtain the environment variables
#this is useful for reference and troubleshooting issues.
env >               /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-all.out
env | grep SLURM >  /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-SLURM.out


#start time
date

#start the simulation
srun ./fsc26 -t apcl_boot.tpl -e apcl_boot.est -n 100000 -m -u -M -L 40 -0 --numBatches 1 --cores 1 --initValues apcl_boot.pv  > /scratch/$USER/Aclarkii_Demo/bootstrapSFS/$DIR/slurm-out/$SLURM_JOBID/$SLURM_JOBID-fastsimcoal26.$SLURM_NNODES.$SLURM_NODELIST.$date.a.txt
```

Once the model was running well, modified `run_sbatch.sb` to submit the previous script 20 times in a job array (need multiple replicate runs to pick best parameter estimates from).

-   This results in 20 replicate runs for each bootstrapped MSFS.

Once all replicate runs for the bootstrapped data finished, wrote `cat_CIs.sbatch` script to concatenate the results (`.bestlhoods` files) from each bootstrapped MSFS into a single file.

-   This script was modified from <https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/cat_cis>.

``` bash
#!/bin/bash

cat /scratch/rdc129/Aclarkii_Demo/bootstrapSFS/apcl_boot1/slurm-out/*/apcl_boot/apcl_boot.bestlhoods.trimmed > /scratch/rdc129/Aclarkii_Demo/bootstrapSFS/apcl_boot1/cis_sfs_summary

#repeat lines for each boot directory (apcl_boot1:apcl_boot100)
```

Then, wrote `cat_max_summary_CIs.sbatch` script to find the parameters from the ML run for each of the 100 bootstrapped MSFS and concatenate into `fsc_maxLhood_CI_summary.txt`.

-   This script was modified from <https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/cat_max_summary>.

``` bash
#!/bin/bash

for file in /scratch/rdc129/Aclarkii_Demo/bootstrapSFS/*/cis_sfs_summary
do

#Sorts the file by column 6 (likelihood) in descending order and then prints the last row (maximum since it's negative)
sort -nk 7 "$file" | tail -n 1 > "$file.max"

#Concatenates the run with the maximum likelihood for each SFS into a file
cat /scratch/rdc129/Aclarkii_Demo/bootstrapSFS/*/cis_sfs_summary.max > /scratch/rdc129/Aclarkii_Demo/bootstrapSFS/fsc_maxLhood_CI_summary.txt

done
```

Copied `fsc_maxLhood_CI_summary.txt` file to local computer, converted to `.csv`, and read into R to calculate 95% CIs with `fsc_CIs.R`.
