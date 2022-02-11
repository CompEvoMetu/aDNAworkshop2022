#!/bin/bash -l
#SBATCH -p chimp
#SBATCH -n 3
#SBATCH -t 5-00:00:00
#SBATCH -J slim
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

simfile=$1
outfile=$2

#seed=$RANDOM
/usr/local/sw/SLiM-3.7/build/slim -d seed=$RANDOM ${simfile} > ${outfile}
