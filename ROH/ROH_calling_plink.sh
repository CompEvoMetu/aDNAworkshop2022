#!/bin/bash -l
#SBATCH -p gorilla
#SBATCH -n 3
#SBATCH -t 360:00:00
#SBATCH -J ROH_calling
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

bamfile=$1
bcftools=/usr/local/sw/bcftools-1.9/bcftools
plink=/usr/local/sw/plink-1.90/plink
bedfile=/mnt/NEOGENE1/toTransfer/workshopfiles/roh/1240K_v2.bed
ref=/mnt/NEOGENE3/share/ref/genomes/hsa/hs37d5.fa
filebase=$(basename ${bamfile} .bam)

# Long streches of autozygous (i.e. homozygous due to inheritance from a recent common ancestor) regions on the genome of an individual are called runs of homozygosity (ROH).
# ROH can be detected using various methods and software (Ceballos et al., 2018).
# Using PLINK is one of the best available options to find ROH. It's observational (i.e., not model based), relatively simple and robust.
# Here we describe a pipeline to call ROH, starting from bam files as the first type of input files.

# SNP calling step: bam to bcf
samtools mpileup -B -q30 -Q30 \ # minimum base (Q) and map (q) qualities set
-l  ${bedfile} \ # a bed file with SNP positions is provided, in this example a 1240K.bed file
-uf ${ref} \ # a reference genome file provided
${bamfile} > ${filebase}.bcf # the available SNPs on the reads of the bam files are called and stored on a bcf

${bcftools} call \
-mV indels \ # to skip indel sites
${filebase}.bcf > ${filebase}.vcf # a vcf is generated from a bcf

#vcf to PLINK: PLINK input files are generatted using the vcfs

${plink} \
--vcf ${filebase}.vcf \
--make-bed \
--const-fid Sample \
--out ${filebase}

${plink}  \
--bfile ${filebase} \
--homozyg-window-snp 30 \ # number of SNPs that the sliding window must contain (30 SNPs)
--homozyg-window-het 1 \ # number of heterozygous SNPs allowed in a window (0 here, but ususally 1 is a better option)
--homozyg-snp 30 \ # minimum number of SNPs that a ROH is required to contain (30 SNPs)
--homozyg-kb 500 \ # length in Kb of the sliding window (500 Kb)
--homozyg-density 30 \ # required minimum density to consider a ROH (1 SNP in 30 Kb)
--out ${filebase}

# These are other parameters that can be set. These are already at the default values PLINK uses, so they are not included in the script
# --homozyg-gap 1000. Length in Kb between two SNPs to be considered in two different segments (1 Mb)
# --homozyg-window-missing 5. Number of missing calls allowed in a window (5 calls)
# --homozyg-window-threshold 0.05. Proportion of overlapping windows that must be called homozygous to define a given SNP as
# in a “homozygous” segment (5%)

# for i in $(cat /mnt/NEOGENE1/toTransfer/workshopfiles/roh/ROHbamlist.txt); do sbatch ROH_calling.sh $i; done
