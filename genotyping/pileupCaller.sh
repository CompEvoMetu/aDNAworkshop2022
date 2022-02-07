#!/bin/bash -l
#SBATCH -p chimp
#SBATCH -t 5-00:00:00
#SBATCH -J snpcall
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err
#SBATCH -n 1

pileupCaller=/usr/local/sw/sequenceTools-1.4.0/pileupCaller
samtools=/usr/local/sw/samtools-1.9/samtools
mergeit=/usr/local/sw/EIG-7.2.1/bin/mergeit

ref=/mnt/NEOGENE3/share/ref/genomes/hsa/hs37d5.fa
pos=/mnt/NEOGENE3/share/dataset/positionW.bed
eigen=/mnt/NEOGENE3/share/dataset/ho.snp
dataset=/mnt/NEOGENE3/share/dataset

bamlist=$1

##snpcalled with pileupcaller
$samtools mpileup -R -B -q 30 -Q 30 -l ${pos} -f ${ref} $(cat ${bamlist} | xargs) > ${bamlist}.pileup.txt

sampleNames=$(grep -o '[^/]*$' ${bamlist} | cut -d "." -f1 | cut -f1,2 -d "_" | tr '\n' ',' | sed 's/,$//')
echo $sampleNames

$pileupCaller --randomHaploid --sampleNames ${sampleNames}  -f ${eigen} -e ${bamlist}.out. < ${bamlist}.pileup.txt

##merged with dataset
cat > parmerge <<EOF
geno1:            ${dataset}/ho.geno
snp1:             ${dataset}/ho.snp
ind1:             ${dataset}/ho.ind
geno2:            ${bamlist}.out.geno.txt
snp2:             ${bamlist}.out.snp.txt
ind2:             ${bamlist}.out.ind.txt
genooutfilename:  merged_dataset.geno
snpoutfilename:   merged_dataset.snp
indoutfilename:   merged_dataset.ind
strandcheck:      NO
EOF

${mergeit} -p parmerge
