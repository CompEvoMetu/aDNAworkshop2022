#!/bin/bash -l
#SBATCH -p chimp
#SBATCH -n 3
#SBATCH -t 5-00:00:00
#SBATCH -J dnSNP_angsd
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

bamlist=$1
outname=$2

angsd=/usr/local/sw/angsd/angsd
ref=/mnt/NEOGENE3/share/ref/genomes/hsa/hs37d5.fa

${angsd}/angsd -bam ${bamlist} -out ${outname} -doPlink 2  -doGeno -4 -doPost 1 -doMajorMinor 4 -GL 1 -doCounts 1 -doMaf 1 -postCutoff 0.99 -SNP_pval 1e-6 -geno_minDepth 3 -geno_maxDepth 13 -nthreads 3 -ref ${ref} -uniqueOnly -1
