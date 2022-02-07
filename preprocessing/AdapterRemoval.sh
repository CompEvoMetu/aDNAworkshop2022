#!/bin/bash -l
#SBATCH -p gorilla
#SBATCH -n 4
#SBATCH -t 10-00:00:00
#SBATCH -J adapter
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

Forward=$1
Reverse=$2
Outdir=$3
Library_out=$4
Cores=$5

AdapterRemoval=/usr/local/sw/adapterremoval-2.3.1/build/AdapterRemoval

echo "Merge start"
$AdapterRemoval --file1 $1 --file2 $2 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT --qualitybase 33 --gzip --qualitymax 60 --trimns --collapse --minalignmentlength 11 --threads ${Cores} --basename $3/$4 --settings $4.settings

echo "Merge end"

cat $3/$4.collapsed.gz $3/$4.collapsed.truncated.gz $3/$4.pair1.truncated.gz $3/$4.pair2.truncated.gz > $3/$(basename $4 .fastq.gz).all.fastq.gz
echo "Both files done"
rm $3/$4.collapsed.gz $3/$4.collapsed.truncated.gz $3/$4.pair1.truncated.gz $3/$4.pair2.truncated.gz $3/$(basename $4 .fastq.gz).discarded.gz $3/$(basename $4 .fastq.gz).singleton.truncated.gz
echo "Remove files done"
