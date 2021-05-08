StairwayPlot\_upstream
================

Code for Stairway Plot v.2 analyses run on an HPC cluster (Amarel at Rutgers) before transferring output files to local computer for downstream analyses.

For Stairway Plot, had to first generate the observed folded site frequency spectrum (SFS). To main neutrality and retain rare variants, built the observed SFS with a dataset unfiltered for minor allele counts but without outlier SNPs (5,662 SNPs). All individuals were included.

**To generate SFS**

The folded SFS for each population was calculated by hand based on the allele frequency of the reference allele for each SNP in each sampling site

To calculate the folded allele frequency at each SNP in each sampling site (see `Data/allele_freqs_all.xlsx`):

1.  If the reference allele frequency was &gt; 0.5, it was subtracted from 1.
2.  The folded allele frequency was then multiplied by 2N (N = sample size of each site).
3.  The number of observations of each allele count was then calculated (ex. how many SNPs had folded allele count of 0, 1, 2, 3, etc.). These counts were then used as the SFS in each population.

**To run Stairway Plot**

Downloaded Stairway Plot v.2 from <https://github.com/xiaoming-liu/stairway-plot-v2> to Amarel and unzipped/installed in `~/Programs`.

``` bash
unzip stairway_plot_v2.1.1.zip
```

Created blueprint file (config/input file) to run Stairway Plot for each population. Example for Japan (`Japan_fold.blueprint`) is below:

``` bash
#Japan blueprint file
#input setting
popid: Japan_fold # id of the population (no white space)
nseq: 16 # number of sequences (2N)
L: 1199553 # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: true # whether the SFS is folded (true or false)
SFS: 1434 586 361 289 175 223 208 188 # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space) (should be N/2 values)
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2
#largest_size_of_SFS_bin_used_for_estimation: 15 # default is nseq/2 for folded SFS
pct_training: 0.67 # percentage of sites for training
nrand: 4        7       11      14 # number of random break points for each try (separated by white space)
project_dir: Japan_fold # project directory
stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
#random_seed: 6
#output setting
mu: 1e-8 # assumed mutation rate per site per generation
year_per_generation: 5 # assumed generation time (in years)
#plot setting
plot_title: Japan_fold # title of the plot
xrange: 0.1,100 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size
```

Notes:

-   **L:** length of sequence (total number of obesrved nucleic sites both polymorphic and monomorphic). Here, this was the total \# of bp in the contigs included in the analysis.
-   **nrand:** \# of random break points for each try. Used 4 numbers ((nseq-2)/4, (nseq-2)/2, ((nseq-2)\*3)/4, nseq-2), suggested in manual.
-   **mu:** mutation rate per site per generation, same as used in *fastsimcoal2*.
-   **year\_per\_generation:** assumed generation time (years), same as used in *fastsimcoal2*.

Ran Stairway Plot:

``` bash
#first make batch file
java -cp stairway_plot_es Stairbuilder Japan_fold.blueprint

#run batch file
bash Japan_fold_blueprint.sh
```

Once done, made output folder `Japan_fold` with diff estimations of demographic histories (one for each of the rand values put in the blueprint file). The `Japan_fold_final.summary.png` contains the best estimate.

Repeated for each sampling site. Copied `*_fold_final.summary.png` to local computer for visualization.
