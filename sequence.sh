#!/bin/bash -l
#SBATCH -p bonobo
#SBATCH -n 1
#SBATCH -t 5-00:00:00
#SBATCH -J seq
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err



ref=/mnt/NEOGENE3/share/ref/genomes/hsa/hs37d5.fa
alndir=/mnt/NEOGENE3/share/compevo_rawdata/hsa/mappedlibraries

mergedname=$1
cores=$2
out=$3 # directory for seqdat

# prerequisites
# samtools Version: 1.9 (using htslib 1.9)
scriptdir=/mnt/NEOGENE1/script/adna_swPipe
genomeCoverageBed=/usr/local/sw/bedtools2-2.25.0/bin/genomeCoverageBed
pmdtools=/usr/local/sw/PMDtools/pmdtools.0.60.py
MTname=MT
Xchr=X
Ychr=Y

####
filebase=$(basename $mergedname .fastq.gz)
refbase=$(basename $ref .fa)
####

####
# to create temporary folder
sample="$( cut -d '.' -f 1 <<< "${mergedname}" )"; samplename=${sample##*/}; echo "$samplename"

TMPDIR=${out}/tmp_${samplename}
mkdir ${TMPDIR} -p
## seq stats
####### Number of sequences
seq_num=$(zcat ${mergedname} | wc -l)
let seq_num=${seq_num}/4
####### Number of mapped sequences and Human Proportion
mapped_count=$(samtools view -c -q 30 -F 4 ${alndir}/mapped/${filebase}.${refbase}.bam)
hum_prop=$(echo ${mapped_count}/${seq_num} | bc -l)
####### Number of effective mapped sequences and Human Proportion
mapped_count_eff=$(samtools view -c -q 30 -F 4 ${alndir}/${filebase}.${refbase}.cons.90perc.bam)
hum_prop_eff=$(echo ${mapped_count_eff}/${seq_num} | bc -l)
####### Avg Read Length
########## Raw bam - RL
samtools view -q 30 ${alndir}/mapped/${filebase}.${refbase}.bam | awk '{print length($10)}'| sort -T ${TMPDIR} | uniq -c | sed "s/^ \+//g; s/ \+/\t/g" | cut -f 1,2 > ${out}/${filebase}.${refbase}.length
Rscript ${scriptdir}/readlength_plot.R ${out} ${filebase}.${refbase}.length
########## Final bam - RL
samtools view -q 30 ${alndir}/${filebase}.${refbase}.cons.90perc.bam | awk '{print length($10)}'| sort -T ${TMPDIR} | uniq -c | sed "s/^ \+//g; s/ \+/\t/g" | cut -f 1,2 > ${out}/${filebase}.${refbase}.cons.90perc.length
Rscript ${scriptdir}/readlength_plot.R ${out} ${filebase}.${refbase}.cons.90perc.length
mean_RL=$(awk '{SUM+=$1*$2; N+=$1}END{print SUM/N}' ${out}/${filebase}.${refbase}.cons.90perc.length)
echo "RL plot done"
####### PMD
samtools view -q 30 ${alndir}/${filebase}.${refbase}.cons.90perc.bam | python $pmdtools --deamination > ${out}/${filebase}.${refbase}.pmd.txt
### Add PMD plot script here!!!!!!
tDamage=$(sed -n '2p' ${out}/${filebase}.${refbase}.pmd.txt | cut -f 2) || tDamage=NA
aDamage=$(sed -n '2p' ${out}/${filebase}.${refbase}.pmd.txt | cut -f 6) || aDamage=NA
echo "Damage plot done, calc"
####### Perc short & Clonality
short=$(cat ${TMPDIR}/short.txt)
perc_short=$(echo 100*${short}/${mapped_count} | bc -l)
mapped_count_nodup_noshort=$(samtools view -c -q 30 -F 4 ${alndir}/${filebase}.${refbase}.cons.90perc.bam)
old_clonality=$(echo 100-100*${mapped_count_nodup_noshort}/${mapped_count} | bc -l)
act_clonality=$(echo ${old_clonality}-${perc_short} | bc -l)
echo "Clonality calc done"
####### Coverage
samtools view -b -q 30 -F 4 ${alndir}/${filebase}.${refbase}.cons.90perc.bam | $genomeCoverageBed -ibam - -g ${ref} > ${TMPDIR}/${filebase}.cov
coverage=$(grep genome ${TMPDIR}/${filebase}.cov | awk '{NUM+=$2*$3; DEN+=$3} END {print NUM/DEN}')
coverage_mt=$(grep ${MTname} ${TMPDIR}/${filebase}.cov | awk '{NUM+=$2*$3; DEN+=$3} END {print NUM/DEN}')
echo "Coverage done"
####### samtools idxstats
samtools idxstats ${alndir}/${filebase}.${refbase}.cons.90perc.bam > ${TMPDIR}/${filebase}.idxstats
hum_seqs=$(cat ${TMPDIR}/${filebase}.idxstats | awk '{ sum+=$3} END {print sum}')
mt_seqs=$(grep ${MTname} ${TMPDIR}/${filebase}.idxstats | awk '{ sum+=$3} END {print sum}') || mt_seqs=NA
x_seqs=$(grep ${Xchr} ${TMPDIR}/${filebase}.idxstats | awk '{ sum+=$3} END {print sum}') || x_seqs=NA
y_seqs=$(grep ${Ychr} ${TMPDIR}/${filebase}.idxstats | awk '{ sum+=$3} END {print sum}') || y_seqs=NA
echo "num seq done"
####### Sex determination
samtools view -q 30 ${alndir}/${filebase}.${refbase}.cons.90perc.bam | python ${scriptdir}/ry_compute.py > ${TMPDIR}/${filebase}.ry.txt
bsex=$(cut -f 7 ${TMPDIR}/${filebase}.ry.txt | tail -n 1)
echo "Sex determination done"
date=`expr match "${filebase}" '.*\(1[1-9][0-1][0-9][0-3][0-9]\)'` || date="nodate"
echo "Date done"
####### MT haplogroup
MThaplo=NA
####### Y haplogroup
Yhaplo=NA


flock -x ${filebase}.sequence_stats.txt echo -e "${filebase:0:6},${filebase:0:15},${filebase},${date},${seq_num},${hum_seqs},${hum_prop},${hum_prop_eff},${mean_RL},${act_clonality},${perc_short},${coverage},${coverage_mt},${mt_seqs},${x_seqs},${y_seqs},${bsex},${tDamage},${aDamage},${MThaplo},${Yhaplo}" > ${filebase}.sequence_stats.txt

----------
###running 

sbatch sequence.sh /mnt/NEOGENE3/share/compevo_rawdata/hsa/mergedfastqs/zoc007_b1e1l1p1_GATATTG-AACGAAG_L002_ARmerged.211116_A00689_0397_AHM2K5DRXY.all.fastq.gz 1 /mnt/NEOGENE1/home/${username}
