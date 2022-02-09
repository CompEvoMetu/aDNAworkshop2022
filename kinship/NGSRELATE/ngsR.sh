#!/bin/bash -l
#SBATCH -p gorilla
#SBATCH -n 1
#SBATCH -t 20-00:00:0
#SBATCH -J ngsR
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err
​
​
DIRngs=/usr/local/sw/NgsRelate/ngsRelate/ngsRelate
DIRout=~/workshop/kinship/ngsRelate
DIRangsd=/mnt/NEOGENE1/toTransfer/workshopfiles/kinship/ngsRelate/angsd
​
##Go to the ngs directory
cd $DIRout/ngs
​
##Prepare frequency file by matching filtered positions from angsd output to the positions from the plink MAF file of our Anatolia dataset
​
$DIRngs extract_freq_bim $DIRangsd/angsd.out.glf.pos.gz $DIRout/dataset/yri.simon_Anatolia.bim $DIRout/dataset/yri.simon_Anatolia_freq.frq > $DIRout/ngs/freq
​
##Run NgsRelate
​
# Genotype likelihoods     --> -g angsd.out.glf.gz
# Number of individuals    --> -n 11
# Frequency file           --> -f freq
# Pairs IDs                --> -z pairs
​
$DIRngs -g $DIRangsd/angsd.out.glf.gz -n 11 -f $DIRout/ngs/freq -z $DIRout/ngs/pairs -O $DIRout/ngs/ngsRelate.res
​
# -> Frequency file: '~/workshop_2022/kinship/ngsRelate/ngs/freq' contain 5764784 number of sites
# -> nind:11 overall_number_of_sites:5764784
​
​
##Look at the results
​
#column 3 : first individual of the pair
#column 4 : second individual of the pair
#column 5 : number of overlap SNP sites for a pair
#column 6 : J9 = k0
#column 7 : J8 = k1
#column 8 : J7 = k2
#column 18 : theta estimate for a pair
​
cut -f3,4,5,6,7,8,18 $DIRout/ngs/ngsRelate.out.newres
​
#Related Pairs :
#Bon004 - Bon005
#Ash128 - Ash133
#Ash131 - Ash136
​
##Visualization using R
Rscript $DIRout/ngs/workshop_2022_ngs.R
