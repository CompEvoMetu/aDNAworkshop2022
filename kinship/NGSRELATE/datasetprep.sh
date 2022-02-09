#!/bin/bash -l
#SBATCH -p gorilla
#SBATCH -n 1
#SBATCH -t 20-00:00:0
#SBATCH -J ngskinship
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err
​
#We will use 11 Anatolian Neolithic individuals to estimate kinship coefficient with NgsRelate.
#Bon001,Bon002,Bon004,Bon005,Ash002,Ash033,Ash128,Ash129,Ash131,Ash133,Ash136
​
DIRplink=/usr/local/sw/plink-1.90/plink
​
mkdir -p ~/workshop/kinship/ngsRelate
DIRout=~/workshop/kinship/ngsRelate
​
cp -r /mnt/NEOGENE1/toTransfer/workshopfiles/kinship/ngsRelate/* $DIRout
​
​
##Go to the dataset directory
cd $DIRout/dataset
​
​
##Anatolian Neolithic individuals (n=59) from Simon Yoruba dataset to calculate population allele frequencies
#Asikli,Boncuklu,Çatalhöyük,Tepecik-Çiftlik,Barcın
​
​
##Get the individual IDs from the dataset and keep it in a file
grep -w -f $DIRout/dataset/indv /mnt/NEOGENE3/share/dataset/hsa/yri.simon.2021.comp.fam > $DIRout/dataset/keep.fam
​
##Create Anatolian dataset in plink format from those individuals to calculate minor allele frequencies for NgsRelate
​
#create bed file only from those individuals         --> --keep keep.fam
#remove all variants with missing call rates > 0.95  --> --geno 0.95
#recode SNP alleles from 1,2,3,4 to A,C,G,T          --> --allele-ACGT
#keep only autosomes                                 --> --autosome
#output format is bed                                --> --make-bed
​
$DIRplink --bfile /mnt/NEOGENE3/share/dataset/hsa/yri.simon.2021.comp --keep $DIRout/dataset/keep.fam --geno 0.95 --allele-ACGT --autosome --make-bed --out $DIRout/dataset/yri.simon_Anatolia
​
#Output:
#.bed
#.bim
#.fam
#.log
#.nosex
​
##Generate a file containing major/minor alleles and MAF values from Yoruba Anatolian dataset to use later for NgsRelate
$DIRplink --bfile $DIRout/dataset/yri.simon_Anatolia --freq --out $DIRout/dataset/yri.simon_Anatolia_freq
​
#Output:
# Minor allele frequency                   : .frq  [CHR SNP A1(minor) A2(major) MAF NCHROBS (# of allele observations) ]
# Log file                                 : .log
# List of samples with ambiguous sex codes : .nosex
