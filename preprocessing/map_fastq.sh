#!/bin/bash -l
#SBATCH -p gorilla
#SBATCH -n 56
#SBATCH -t 5-00:00:00
#SBATCH -J mappping
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

mergedname=$1 #fastq file
ref=$2 #ref genome
alndir=$3 # directory for output
cores=$4 #number of thread

# prerequisites
# samtools Version: 1.9 (using htslib 1.9)
bwa=/usr/local/sw/bwa-0.7.15/bwa
scriptdir=/mnt/NEOGENE1/script/adna_swPipe
MTname=MT
Xchr=X
Ychr=Y

####
filebase=$(basename $mergedname .fastq.gz)
refbase=$(basename $ref .fa)
mkdir -p mapped #to collect non-filtered bam files
####

####
# to create temporary folder
sample="$( cut -d '.' -f 1 <<< "${mergedname}" )"; samplename=${sample##*/}; echo "$samplename"

TMPDIR=tmp_${samplename}
mkdir ${TMPDIR} -p
echo ${TMPDIR}
####

#### align to ref genome
$bwa aln -l 16500 -n 0.01 -o 2 -t ${cores} ${ref} ${mergedname} \
        | $bwa samse ${ref} - ${mergedname} \
        | samtools view -F 4 -h -Su - \
        | samtools sort - -@ ${cores} -o ${alndir}/mapped/${filebase}.${refbase}.bam
echo "bwa align done"
####

#### FilterUniqueSAMCons
samtools view -F 4 -h -@ ${cores} ${alndir}/mapped/${filebase}.${refbase}.bam \
        | python ${scriptdir}/FilterUniqueSAMCons_rand.py | samtools view -h -Su -@ ${cores} - \
        > ${alndir}/mapped/${filebase}.${refbase}.cons.bam  #To remove bias _rand.py
echo "cons done"
sleep 5s
flock -x ${alndir}/mapped/${filebase}.${refbase}.cons.bam samtools index ${alndir}/mapped/${filebase}.${refbase}.cons.bam
echo "cons index done"

#### Percidentity threshold
samtools calmd -@ ${cores} ${alndir}/mapped/${filebase}.${refbase}.cons.bam ${ref} \
        | python ${scriptdir}/percidentity_threshold.py 0.9 35 ${TMPDIR}/short.txt \
        | samtools view -bS - > ${alndir}/${filebase}.${refbase}.cons.90perc.bam
echo "90perc done"
sleep 5s
flock -x ${alndir}/${filebase}.${refbase}.cons.90perc.bam samtools index ${alndir}/${filebase}.${refbase}.cons.90perc.bam
echo "90perc index done"
echo "mapping done"
