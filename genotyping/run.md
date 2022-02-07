# Genotyping

## De novo SNP Call:

```bash
ls /mnt/NEOGENE1/toTransfer/workshopfiles/genotyping/files/Ash128.trimBAM.chr7.bam > Ash128.bamlist.txt

sbatch dnSNPCall.sh Ash128.bamlist.txt Ash128.dnSNPCall
```

## PileupCaller

```bash
ls /mnt/NEOGENE3/share/dna/hsa/trimmedbams/Ash136_all.merged.hs37d5.fa.cons.90perc.trimBAM.bam > Ash136.bamlist.txt

sbatch pileupCaller.sh Ash136.bamlist.txt 
```
