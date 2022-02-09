#!/bin/bash -l
#SBATCH --partition=macaque1                # write partition name(chimp;bonobo;macaque$i)
#SBATCH --job-name=fst                 # write your job name
#SBATCH --ntasks=1                   # number of task
#SBATCH --output=/mnt/NEOGENE1/home/NAME/slurm_output.txt
#SBATCH --error=/mnt/NEOGENE1/home/NAME/slurm_output.err


echo "SLURM_NODELIST $SLURM_NODELIST"
echo "NUMBER OF CORES $SLURM_NTASKS"

for i in africa eas eurasia; do	cut -f 1 {i}.tsv | awk 'FNR>1 {print$0}'>{i}.txt; done

/usr/local/sw/vcftools-0.1.16/bin/vcftools --gzvcf chr2.1000genome.phase3.v5a.ALL.filtered.recode.vcf.gz --weir-fst-pop eas.txt --weir-fst-pop eurasia.txt --fst-window-size 10000 --fst-window-step 5000 --out eas_euas_norep --exclude-bed chr2_sorted_rep.bed

/usr/local/sw/vcftools-0.1.16/bin/vcftools --gzvcf chr2.1000genome.phase3.v5a.ALL.filtered.recode.vcf.gz --weir-fst-pop africa.txt --weir-fst-pop eurasia.txt --fst-window-size 10000 --fst-window-step 5000  --out af_euas_norep --exclude-bed chr2_sorted_rep.bed

/usr/local/sw/vcftools-0.1.16/bin/vcftools --gzvcf chr2.1000genome.phase3.v5a.ALL.filtered.recode.vcf.gz --weir-fst-pop africa.txt --weir-fst-pop eas.txt --fst-window-size 10000 --fst-window-step 5000 --out af_eas_norep --exclude-bed chr2_sorted_rep.bed

R CMD BATCH PBS.R
exit
