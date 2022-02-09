#!/bin/bash -l
#SBATCH -p gorilla
#SBATCH -n 15
#SBATCH -t 20-00:00:0
#SBATCH -J angsd
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err
​
DIRang=/usr/local/sw/angsd/angsd/angsd
DIRout=~/workshop/kinship/ngsRelate
​
​
#Go to the angsd directory
cd $DIRout/angsd
​
#Prepare "sites" that are SNP positions we want to use for genotype likelihood estimation with angsd
#"sites" file contains 4 columns : chr pos major minor
cut -f1,4-6 $DIRout/dataset/yri.simon_Anatolia.bim > $DIRout/angsd/angsd.sites
​
#Indexing the angsd sites file
​
$DIRang sites index $DIRout/angsd/angsd.sites
​
#Output:
#angsd.sites.idx
#angsd.sites.bin
​
#Run angsd to calculate genotype likelihoods from the preferred sites
​
#list of bam files                                            --> -bam bamlist
#sites                                                        --> -sites angsd.sites
#mapping quality filter >30                                   --> -minMapQ 30
#base quality filter >20                                      --> -minQ 20
#15 threads                                                   --> -nThreads 15
#Genotype likelihood estimation using SAMtools method         --> -GL 1
#Output the log genotype likelihoods to beagle binary file    --> -doGlf 3
#Use major and minor from a file (requires -sites file)       --> -doMajorMinor 3
#Use known major and minor alleles to calculate MAF           --> -doMAF 1
​
​
$DIRang -bam $DIRout/angsd/bam.list -GL 1 -doGlf 3 -doMajorMinor 3 -doMAF 1 -sites $DIRout/angsd/angsd.sites -nThreads 15 -minMapQ 30 -minQ 20 -out $DIRout/angsd/angsd.out
​
#Output :
#Positions and major/minor alleles          : angsd.out.glf.pos.gz
#Genotype likelihoods                       : angsd.out.glf.gz
#Calculated MAF                             : angsd.out.mafs.gz
#Arguments and parameters for all analysis  : angsd.out.arg
​
# -> Total number of sites analyzed: 2786164931
# -> Number of sites retained after filtering: 5764784
