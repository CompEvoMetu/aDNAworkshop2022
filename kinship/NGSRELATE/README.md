# NgsRelate
​
NgsRelate is a program used for relatedness inference, inbreeding coefficients and many other summary statistics for pairs of individuals from low coverage Next Generation Sequencing (NGS) data by using population allele frequencies and genotype likelihoods instead of called genotypes.
​
**Links** : 
​
[NgsRelate: a software tool for estimating pairwise relatedness from next-generation sequencing data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4673978/)
​
[Fast and accurate relatedness estimation from high-throughput sequencing data in the presence of inbreeding](https://academic.oup.com/gigascience/article/8/5/giz034/5481763?login=false)
​
[ANGSD/NgsRelate](https://github.com/ANGSD/NgsRelate)
​
[NgsRelate old version](http://www.popgen.dk/software/index.php?title=NgsRelate&oldid=694)
[ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD#Overview)
​
### **Dataset Preparation**
​
---
​
We will use 11 Anatolian Neolithic individuals to estimate kinship coefficient with NgsRelate.
Bon001,Bon002,Bon004,Bon005,Ash002,Ash033,Ash128,Ash129,Ash131,Ash133,Ash136
​
```bash
#Make a directory in your home folder 
mkdir -p ~/workshop/kinship/ngsRelate
​
#Specify which directory you will be working 
DIRout=~/workshop/kinship/ngsRelate
```
​
```bash
#Copy the files and scripts to your home folder 
cp -r /mnt/NEOGENE1/toTransfer/workshopfiles/kinship/ngsRelate/* $DIRout
​
#Go to your working directory  
cd $DIRout
```
​
NgsRelate requires population allele frequencies for kinship estimation. We will use Anatolian Neolithic individuals (n=59) from Simon Yoruba dataset to calculate population allele frequencies. These individuals are from Aşıklı, Boncuklu, Çatalhöyük, Tepecik-Çiftlik, and Barcın.
​
```bash
DIRplink=/usr/local/sw/plink-1.90/plink
​
#Create a file containing IDs of the Anatolian Neolithic individuals (n=59) from yri.simon.2021.comp
head $DIRout/dataset/indv
Ash002
Ash033
Ash040
Ash128
Ash129
Ash131
Ash133
Ash136
cth006
cth728 
​
#Find the individual IDs from the Yoruba Simon dataset and keep them in a fam file
grep -w -f $DIRout/dataset/indv /mnt/NEOGENE3/share/dataset/hsa/yri.simon.2021.comp.fam > $DIRout/dataset/keep.fam
​
head $DIRout/dataset/keep.fam
420   Bar31 0 0 0 1
421   Bar8 0 0 0 1
425   Bon001 0 0 0 1
426   Bon002 0 0 0 1
427   Bon004 0 0 0 1
428   Bon005 0 0 0 1
740   I0707 0 0 0 1
741   I0708 0 0 0 1
742   I0709 0 0 0 1
748   I0736 0 0 0 1
​
#Create Anatolian dataset in plink format from those individuals to calculate minor allele frequencies
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
```
​
### Population Allele Frequency Calculation
​
```bash
#Generate a file containing major/minor alleles and MAF values from Yoruba Anatolian dataset 
​
#minor allele frequencies  --> --freq
 
$DIRplink --bfile $DIRout/dataset/yri.simon_Anatolia --freq --out $DIRout/dataset/yri.simon_Anatolia_freq
​
#Output:
#Minor allele frequency :                   .frq [CHR SNP A1(minor) A2(major) MAF NCHROBS (# of allele observations)]
#Log file :                                 .log
#List of samples with ambiguous sex codes : .nosex
​
head $DIRout/dataset/yri.simon_Anatolia_freq.frq
CHR          SNP   A1   A2          MAF  NCHROBS
   1      1_30923    G    T            0       16
   1      1_52238    T    G            0       12
   1      1_54716    T    C       0.2857       14
   1      1_55545    T    C            0       22
   1      1_63268    C    T            0       10
   1      1_67181    G    A            0        6
   1      1_79137    T    A            0        8
   1      1_83084    T    A            0       10
   1      1_86331    G    A            0       14
```
​
### Genotype Likelihood Estimation
​
NgsRelate also requires genotype likelihoods, so the next step is the estimation of genotype likelihoods by angsd. We can select biallelic SNP positions from our new Yoruba Anatolian dataset to use for the estimation. 
​
```bash
#Prepare "sites" that are SNP positions we want to use for genotype likelihood estimation with angsd
cut -f1,4-6 $DIRout/dataset/yri.simon_Anatolia.bim > $DIRout/angsd/angsd.sites
​
#"sites" file contains 4 columns : chr pos major minor
head $DIRout/angsd/angsd.sites
1       30923   G       T
1       52238   T       G
1       54716   T       C
1       55545   T       C
1       63268   C       T
1       67181   G       A
1       79137   T       A
1       83084   T       A
1       86331   G       A
1       88169   T       C
​
```
​
```bash
DIRang=/usr/local/sw/angsd/angsd/angsd
​
#Create a list of bam files including your preferred bam files with their path for the kinship estimation
#head bam.list
​
#Run angsd to calculate genotype likelihoods from the preferred sites
​
#list of bam files                                            --> -bam bam.list
#sites                                                        --> -sites angsd.sites
#mapping quality filter >30                                   --> -minMapQ 30
#base quality filter >20                                      --> -minQ 20
#15 threads                                                   --> -nThreads 15
#Genotype likelihood estimation using SAMtools method         --> -GL 1
#Output the log genotype likelihoods to beagle binary file    --> -doGlf 3
#Use major and minor from a file (requires -sites file)       --> -doMajorMinor 3
#Use known major and minor alleles to calculate MAF           --> -doMAF 1
​
$DIRang -bam $DIRout/angsd/bam.list -GL 1 -doGlf 3 -doMajorMinor 3 -doMAF 1 -sites $DIRout/angsd/angsd.sites -nThreads 15 -minMapQ 30 -minQ 20 -out $DIRout/angsd/angsd.out
​
#Output :
#Positions and major/minor alleles          : angsd.out.glf.pos.gz
#Genotype likelihoods                       : angsd.out.glf.gz
#Calculated MAF                             : angsd.out.mafs.gz
#Arguments and parameters for all analysis  : angsd.out.arg
​
# -> Total number of sites analyzed: 2786164931
# -> Number of sites retained after filtering: 5764784
```
​
Then, we will run NgsRelate. Before starting the analysis, we will extract the overlapping positions between filtered sites with genotype likelihoods (angsd output) and Yoruba Anatolian sites with minor allele frequencies (plink output) in order to make a frequency file.  
​
```bash
​
#Prepare frequency file by matching filtered positions from angsd output to the positions from the plink MAF file of our Yoruba Anatolian dataset
$DIRngs extract_freq_bim $DIRout/angsd/angsd.out.glf.pos.gz $DIRout/dataset/yri.simon_Anatolia.bim $DIRout/dataset/yri.simon_Anatolia_freq.frq > $DIRout/ngs/freq
​
head $DIRout/ngs/freq
1.000000
1.000000
0.714300
1.000000
1.000000
1.000000
1.000000
1.000000
1.000000
0.600000
```
​
### Running NgsRelate
​
Now, we are ready to run NgsRelate using genotype likelihoods of our Anatolian individuals (n=11) and the population allele frequencies of Anatolian individuals (n=59) from Simon Yoruba dataset.
​
```bash
DIRngs=/usr/local/sw/NgsRelate/ngsRelate/ngsRelate
​
#Create a file with individual IDs 
head pairs 
Bon001
Bon002
Bon004
Bon005
Ash002
Ash033
Ash128
Ash129
Ash131
Ash133
​
#Run NgsRelate
​
# Genotype likelihoods     --> -g angsd.out.glf.gz
# Number of individuals    --> -n 11
# Frequency file           --> -f freq
# Pairs IDs                --> -z pairs
​
$DIRngs -g $DIRout/angsd/angsd.out.glf.gz -n 11 -f $DIRout/ngs/freq -z $DIRout/ngs/pairs -O $DIRout/ngs/ngsRelate.res
​
# -> Frequency file: '~/workshop_2022/kinship/ngsRelate/ngs/freq' contain 5764784 number of sites
# -> nind:11 overall_number_of_sites:5764784
```
​
Let’s look at the results.
​
```bash
#column 3 : first individual of the pair
#column 4 : second individual of the pair
#column 5 : number of overlap SNP sites for a pair
#column 6 : J9 = k0
#column 7 : J8 = k1
#column 8 : J7 = k2
#column 18 : theta estimate for a pair
​
cut -f3,4,5,6,7,8,18 $DIRout/ngs/ngsRelate.out.newres
​
#Related Pairs :
#Bon004 - Bon005
#Ash128 - Ash133
#Ash131 - Ash136
```
​
### Visualization using R
​
We can draw a plot using R to interpret our result. 
​
```bash
#Run R script 
Rscript $DIRout/ngs/workshop_2022_ngs.R
```
