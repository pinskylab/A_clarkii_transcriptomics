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
