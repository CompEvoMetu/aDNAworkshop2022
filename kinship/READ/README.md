# aDNA Workshop READ Pipeline

📎 **Usefull links —>** [tguenther / read — Bitbucket](https://bitbucket.org/tguenther/read/src/master/) | [PLINK 1.9 (cog-genomics.org)](https://www.cog-genomics.org/plink/1.9/) | [Estimating genetic kin relationships in prehistoric populations (plos.org)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0195491)


**READ: Relationship Estimation from Ancient DNA**

READ is a method to infer the degree of relationship (up to second degree, i.e. nephew/niece-uncle/aunt, grandparent-grandchild or half-siblings) for a pair of low-coverage individuals.

### **Before starting kinship analysis, create empty folders in your home directory.**

```bash
mkdir ~/workshop
mkdir ~/workshop/kinship/
mkdir ~/workshop/kinship/READ
mkdir ~/workshop/kinship/READ/dataset
mkdir ~/workshop/kinship/READ/read_analysis
mkdir ~/workshop/kinship/READ/scripts
```

### Retrieve required files for READ analysis from designated folders.

```bash
cp /mnt/NEOGENE1/toTransfer/workshopfiles/kinship/READ/scripts/* ~/workshop/kinship/READ/scripts
cp ~/workshop/kinship/READ/scripts/* ~/workshop/kinship/READ/read_analysis/
# You can also download READ scripts via typing wget https://bitbucket.org/tguenther/read/get/f541d553247a.zip
cp /mnt/NEOGENE1/toTransfer/workshopfiles/kinship/READ/dataset/keep ~/workshop/kinship/READ/dataset
```

# Preparing dataset for READ analysis

---

Normally, we must create a new dataset by using individuals' bam files but for this analysis, we will create a subset from a dataset called yri.simon.2021.comp (approximately has 6 million SNPs) which includes our individuals.

```bash
#Enter dataset folder
cd ~/workshop/kinship/READ/dataset

#Create a file that includes individual ids you wish to estimate kinship. An example file has been created previously. You can simply use the file called "keep"
less keep

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
Ash136
```

We will create a subset from *yri.simon.2021.comp* dataset. This process needs index numbers of individual ids'

```bash
#Extract index numbers of individuals from the dataset.
grep -w -f keep /mnt/NEOGENE3/share/dataset/hsa/yri.simon.2021.comp.fam | awk {'print $1, $2'} > keepfile

less keepfile

425 Bon001
426 Bon002
427 Bon004
428 Bon005
2220 Ash002
2221 Ash033
2223 Ash128
2224 Ash129
2225 Ash131
2226 Ash133
2227 Ash136
```

Create a new dataset with selected individuals that includes only autosome chromosomes

```bash
/usr/local/sw/plink-1.90/plink  --bfile /mnt/NEOGENE3/share/dataset/hsa/yri.simon.2021.comp --keep keepfile --make-bed -autosome -out workshop.autosome

#MINOR NOTE: Add "|" character to end of sample ids in the fam file. That's because READ creates pairs with no separator. If you could add a separator it could be helpful for further analysis.
awk 'BEGIN{OFS="\t"}$2=$2"|"' workshop.autosome.fam > modified.fam
mv modified.fam workshop.autosome.fam
```

READ requires Plinks’ transposed TPED/TFAM format. So you must transpose the datasets created in the previous step

```bash
/usr/local/sw/plink-1.90/plink --bfile workshop.autosome -recode transpose --out workshop_READ
```

# Running READ

---

In this section, we will run READ and then examine output files.

```bash
#Enter into read_analysis directory
cd ~/workshop/kinship/READ/read_analysis/

#Start READ analysis
python2 READ.py  ~/workshop/kinship/READ/dataset/workshop_READ > read_run.log

#Calculate a total number of SNPs between pair of two indivuduals
python3 countSNPs.py -f READ_output_ordered > SNPcount.log
mv ~/SNPCount.READ_output_ordered ~/workshop/kinship/READ/read_analysis
```

# Post Analysis

---

If READ is executed successfully, files called *“READ_output_ordered*”, *“READ_results”*, *“Read_intermediate_output”*, *“READ_results_plot.pdf”* and *“meansP0_AncientDNA_normalized”* must be located in  the “read_analysis” folder. 

Estimated relationships in *Yaka et. al., 2021* paper:
ash128-ash133 (Sisters)
ash131-ash136 (Sisters)
Bon004-Bon005 (Brother-Sister)

```bash
#Let's examine every output file that READ created
less read_run.log
less READ_output_ordered
less READ_results
less Read_intermediate_output
less READ_results_plot
less meansP0_AncientDNA_normalized

#Let's look total number of SNPs between investigated pairs. You will realise some pairs have lower than 5000 SNPs.
less SNPCount.READ_output_ordered
```

We must filter pairs whose have lower than 5000 SNPs. Below 5000 SNPs, false positive and false negative relationship degree classifications are relatively higher.

```bash
#Execute post_script.R in your current directory
Rscript post_script.R
```

### Organize and examine result files

---

```bash
#Create folder for raw READ results.
mkdir ~/workshop/kinship/READ/read_analysis/raw_results
mv READ_output_ordered READ_results SNPCount.READ_output_ordered Read_intermediate_output READ_results_plot.pdf meansP0_AncientDNA_normalized ~/workshop/kinship/READ/read_analysis/raw_results

#Create folder for log files.
mkdir ~/workshop/kinship/READ/read_analysis/log 
mv read_run.log SNPcount.log ~/workshop/kinship/READ/read_analysis/log

#Create folder for filtered results.
mkdir ~/workshop/kinship/READ/read_analysis/filtered_results
less READ_results_snpfilter_median_brief
less READ_results_snpfilter_mean_brief
mv READ_results_snpfilter_median READ_results_snpfilter_median_brief READ_plot_snpfilter_median.pdf ~/workshop/kinship/READ/read_analysis/filtered_results
mv READ_results_snpfilter_mean READ_results_snpfilter_mean_brief READ_plot_snpfilter_mean.pdf ~/workshop/kinship/READ/read_analysis/filtered_results
rm Rplots.pdf
```
