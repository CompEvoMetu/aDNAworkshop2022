#!/bin/bash -l
#SBATCH -p bonobo
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J contamination
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

####### prerequisites

ref=/mnt/NEOGENE3/share/ref/genomes/hsa/hs37d5.fa
bwa=/usr/local/sw/bwa-0.7.15/bwa
samtools=/usr/local/sw/samtools-1.9/samtools
angsd=/usr/local/sw/angsd/angsd
muscle=/usr/local/sw/muscle3.8.31_i86linux64
refY=/mnt/NEOGENE3/share/ref/genomes/hsa/hs37d5.ychr.fa
contdeam=/usr/local/sw/schmutzi/src/contDeam.pl
schmutzi=/usr/local/sw/schmutzi/src/schmutzi.pl
refMT=/mnt/NEOGENE3/share/ref/genomes/hsa/rCRS.MT.fasta
freqs=/usr/local/sw/schmutzi/share/schmutzi/alleleFreqMT/eurasian/freqs
estimate=/usr/local/sw/contamMix/exec/estimate.R
green=/mnt/NEOGENE1/script/adna_swPipe/green311.fas
MTname=MT
Xchr=X
Ychr=Y


bam=$1
sample="$( cut -d '.' -f 1 <<< "$bam" )"; samplename=${sample##*/}; echo "$samplename"

mkdir -p $samplename
cd $samplename

####### estimate mt contamination
##create consensus fasta from bamfile for MT (here -doFasta 2 means consensus base)
${angsd}/angsd -i ${bam} -doFasta 2 -doCounts 1 -minQ 30 -minMapQ 30 -setMinDepth 3 -r MT: -out ${samplename}.mt
gzip -d ${samplename}.mt.fa.gz
$bwa index ${samplename}.mt.fa

## mtBAM obtained from bamfile
$samtools view -q 30 -b ${bam} MT > ${samplename}.mt.q30.bam

## mtBAM is alined to consensus fasta
$bwa aln -l 16500 -n 0.01 -o 2 ${samplename}.mt.fa -b ${samplename}.mt.q30.bam | $bwa samse ${samplename}.mt.fa - ${samplename}.mt.q30.bam | $samtools view -F 4 -q 30 -h -Su - | $samtools sort -O bam - > ${samplename}.mt.cons.bam

#contaMix

cat ${samplename}.mt.fa ${green} > green311.${samplename}.fas
##comparing the mapping affinity of each ancient mtDNA read to the reconstructed consensus sequence and to the 311 possible contaminants
$muscle -diags -in green311.${samplename}.fas -out green311.${samplename}.muscle.fas

##the total fraction of non-endogenous reads as the contamination rate with a Bayesian estimate of the posterior probability of the contamination proportion
Rscript $estimate --samFn ${samplename}.mt.cons.bam --alnFn green311.${samplename}.muscle.fas --figure ${samplename}.cont --nIter 100000  --transverOnly | tee ${samplename}.contamMixout.txt

mtMAPauthentic=`grep "MAP authentic" ${samplename}.contamMixout.txt | cut -d":" -f2`
mtErrorRate=`grep "error rate" ${samplename}.contamMixout.txt | awk '{print $9}' | sed 's/).//g'`

### schmutzi
## hasara dayalı kontaminasyon tahmini
${contdeam} --library double --out ${samplename}.conta.bam ${refMT}  ${bam}
## mitokondriyal alel frekanslarının bir veritabanını kullanma
${schmutzi} --ref ${refMT} --out ${samplename}.mtcont.bam ${freqs} ${samplename}.conta.bam

conta=$(cut -f 2 ${samplename}.conta.bam.cont.est)
mtconta=$(cut -f 2 ${samplename}.mtcont.bam_1_cont.est)

####### estimate X contamination
##tek bir BAM dosyası için chrX dosyası oluşturuyoruz
${angsd}/angsd -i ${bam} -r X:5000000-154900000 -doCounts 1  -iCounts 1 -minMapQ 30 -minQ 30 -out ${samplename}.xcont
${angsd}/misc/contamination -a ${samplename}.xcont.icnts.gz -h ${angsd}/RES/HapMapChrX.gz 2> ${samplename}.xcont001.out

##Fisher's exact test ile p değerini bulup, jackknife ile kontaminasyon tahmini
Rscript ${angsd}/R/contamination.R mapFile=${angsd}/RES/chrX.unique.gz hapFile=${angsd}/RES/HapMapChrX.gz countFile=${samplename}.xcont.icnts.gz minDepth=1 mc.cores=72 maxDepth=100 > ${samplename}.xcont002.out

xcont001=`grep new_llh ${samplename}.xcont001.out | head -n 1 | cut -d ":" -f 4-`
xcont002=`grep "Contamination" ${samplename}.xcont002.out | cut -d " " -f 2-3`
