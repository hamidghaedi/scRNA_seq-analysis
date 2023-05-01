# scRNA sequencing analysis

This repo presents steps needed to make sense of single cell RNA sequencing (scRNA) data. I used a scRNA dataset coming from Zhaohui Chen *et al.* [paper](https://www.nature.com/articles/s41467-020-18916-5#Sec12) published in Nature Communications 11, Article number: 5077 (2020). The cohort consisted of eight primary bladder tumor tissues (2 low-grade bladder urothelial tumors, six high-grade bladder urothelial tumors) along with 3 adjacent normal mucosae. In SRA datasets are under BioProject [PRJNA662018](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA662018) and SRA-explorer can be used to download the data. For practical scRNA-seq analysis I followed this elegant [tutorial](https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html) from Harvard Chan Bioinformatics Core. 


## 1) Data download on ComputeCanada

``` shell

#!/usr/bin/env bash

#!/bin/bash
#SBATCH --account=#
#SBATCH --job-name=scRNA__DL
#SBATCH --qos=privileged
#SBATCH --nodes=1                # number of Nodes
#SBATCH --tasks-per-node=4        # number of MPI processes per node
#SBATCH --mem 8g
#SBATCH --time 12:00:00
#SBATCH --output=scRNA_fastq_DL.%J.out
#SBATCH --error=scRNA_fastq.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#


curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/089/SRR12603789/SRR12603789_1.fastq.gz -o SRR12603789_urinary_bladder_cancer_fresh_sample_of_patient_number_02_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/089/SRR12603789/SRR12603789_2.fastq.gz -o SRR12603789_urinary_bladder_cancer_fresh_sample_of_patient_number_02_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/087/SRR12603787/SRR12603787_1.fastq.gz -o SRR12603787_urinary_bladder_cancer_fresh_sample_of_patient_number_03_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/087/SRR12603787/SRR12603787_2.fastq.gz -o SRR12603787_urinary_bladder_cancer_fresh_sample_of_patient_number_03_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/086/SRR12603786/SRR12603786_1.fastq.gz -o SRR12603786_urinary_bladder_cancer_fresh_sample_of_patient_number_04_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/086/SRR12603786/SRR12603786_2.fastq.gz -o SRR12603786_urinary_bladder_cancer_fresh_sample_of_patient_number_04_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/090/SRR12603790/SRR12603790_1.fastq.gz -o SRR12603790_urinary_bladder_cancer_fresh_sample_of_patient_number_01_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/090/SRR12603790/SRR12603790_2.fastq.gz -o SRR12603790_urinary_bladder_cancer_fresh_sample_of_patient_number_01_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/088/SRR12603788/SRR12603788_1.fastq.gz -o SRR12603788_normal_bladder_mucosa_3_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/088/SRR12603788/SRR12603788_2.fastq.gz -o SRR12603788_normal_bladder_mucosa_3_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/084/SRR12603784/SRR12603784_1.fastq.gz -o SRR12603784_urinary_bladder_cancer_fresh_sample_of_patient_number_06_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/084/SRR12603784/SRR12603784_2.fastq.gz -o SRR12603784_urinary_bladder_cancer_fresh_sample_of_patient_number_06_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/085/SRR12603785/SRR12603785_1.fastq.gz -o SRR12603785_urinary_bladder_cancer_fresh_sample_of_patient_number_05_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/085/SRR12603785/SRR12603785_2.fastq.gz -o SRR12603785_urinary_bladder_cancer_fresh_sample_of_patient_number_05_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/083/SRR12603783/SRR12603783_1.fastq.gz -o SRR12603783_urinary_bladder_cancer_fresh_sample_of_patient_number_07_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/083/SRR12603783/SRR12603783_2.fastq.gz -o SRR12603783_urinary_bladder_cancer_fresh_sample_of_patient_number_07_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/082/SRR12603782/SRR12603782_1.fastq.gz -o SRR12603782_urinary_bladder_cancer_fresh_sample_of_patient_number_08_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/082/SRR12603782/SRR12603782_2.fastq.gz -o SRR12603782_urinary_bladder_cancer_fresh_sample_of_patient_number_08_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/081/SRR12603781/SRR12603781_1.fastq.gz -o SRR12603781_normal_bladder_mucosa_1_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/081/SRR12603781/SRR12603781_2.fastq.gz -o SRR12603781_normal_bladder_mucosa_1_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/SRR12603780/SRR12603780_1.fastq.gz -o SRR12603780_normal_bladder_mucosa_2_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/SRR12603780/SRR12603780_2.fastq.gz -o SRR12603780_normal_bladder_mucosa_2_2.fastq.gz
```

It is always recomended to run QC on samples; so uisng the following code we can run `fastqc` on all samples

``` shell
#!/bin/bash
#SBATCH --account= # add your account name 
#SBATCH --job-name=fastqc
#SBATCH --qos=privileged
#SBATCH --nodes=12                # number of Nodes
#SBATCH --tasks-per-node=5        # number of MPI processes per node
#SBATCH --mem 16g
#SBATCH --time 24:00:00
#SBATCH --output=fastqc.%J.out
#SBATCH --error=fastqc.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=# add your email address 


module load fastqc

for f in ./fastq/*.fastq.gz
do
 echo $f
 fastqc $f --outdir ./qc_raw  --thread 12 --nogroup
done
```

There are diffrences between sequence length of read 1 and read 2 for each samples; read 1 provides data on Cell barcode & UMI and read 2 has the insert sequence. Files coming from SRA usually come with names like something_1\_fastq.gz and something_2\_fastq.gz . This [blog](https://kb.10xgenomics.com/hc/en-us/articles/115003802691) is helpful to see what are the naming requirements of fastq files for cellranger tool. Briefly;

-   incompatible file name: SRR9291388_1.fastq.gz

-   compatible file name: SRR9291388_S1_L001_R1_001.fastq.g

So we need to change file names and also set up the directories for `cellranger count`;

``` shell
for file in *.fastq.gz
do
  dir=${file%_*}
  echo $dir
  mkdir $dir
  mv $file ./$dir
done
```

## 2) Generating feature-sample expression matrix

The next step is to run `cellranger count` on each samples. There are different types of scRNA libs with different fastq files. This [blog](https://www.10xgenomics.com/support/single-cell-gene-expression/documentation/steps/sequencing/sequencing-requirements-for-single-cell-3) provide details on type of libs from 10X genomics.

Running `cellranger` on a cluster with `SLURM` as job schaduler is not an easy task. The following is what I came up with working best for me working on ComputeCanada cluster:

``` shell
#!/bin/bash
#SBATCH --account=#
#SBATCH -J test_cellranger
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --signal=2
#SBATCH --no-requeue
### Alternatively: --ntasks=1 --cpus-per-task={NUM_THREADS}
###   Consult with your cluster administrators to find the combination that
###   works best for single-node, multi-threaded applications on your system.
#SBATCH --mem=250G
#SBATCH --time 48:00:00
#SBATCH --output=test_cellranger.%J.out
#SBATCH --error=test_cellranger.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#


module load cellranger

cd /home/ghaedi/projects/def-gooding-ab/ghaedi/sc/raw
for d in *
do
echo $d
ID=${d:0:11}
cellranger count --id=$ID \
                 --transcriptome=/home/ghaedi/projects/def-gooding-ab/ghaedi/sc/refdata-gex-GRCh38-2020-A \
                 --fastqs=$d
                 --chemistry=SC3Pv2
done
```

The next steps are mainly based on Harvard Chan Bioinformatics Core materials on [scRNA-seq analysis](https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html) training.

## 3) Loading single-cell RNA-seq count data

``` r
# Setup the Seurat Object

library(tidyverse)
library(Seurat)
library(patchwork)
library(RCurl)
library(cowplot)

# create list of samples
samples <- list.files("~/scRNA/filtered/")
#samples <- samples[grepl('^filtered',samples,perl=T)]

# read files inot Seurat objects
for (file in samples){
  print(paste0(file))
  seurat_data <- Read10X(data.dir = paste0("~/scRNA/filtered/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

# now merging all objects inot one Seurat obj

merged_seurat <- merge(x = SRR12603780, 
                       y = c(SRR12603781,
                             SRR12603782,
                             SRR12603783,
                             SRR12603784,
                             SRR12603785,
                             SRR12603786,
                             SRR12603787,
                             SRR12603788,
                             SRR12603789,
                             SRR12603790),
                       #Because the same cell IDs can be used for different samples, we add a sample-specific prefix 
                       # to each of our cell IDs using the add.cell.id argument.
                       add.cell.id = samples)
```

## 4) Quality control

There are columns in the metadata:

-   orig.ident: this column will contain the sample identity if known. It will default to the value we provided for the project argument when loading in the data

-   nCount_RNA: represents the number of UMIs per cell. UMI (unique molecular identifiers) is used to determine whether a read is a biological or technical duplicate (PCR duplicate). Reads with different UMIs mapping to the same transcript were derived from different molecules and are biological duplicates - each read should be counted. Reads with the same UMI originated from the same molecule and are technical duplicates - the UMIs should be collapsed to be counted as a single read.

-   nFeature_RNA: represents the number of genes detected per cell

Recommended features to add to metadata:

-   number of genes detected per UMI (or novelty score): more genes detected per UMI, more complex our data

-   mitochondrial ratio: this metric will give us a percentage of cell reads originating from the mitochondrial genes (coming from dying cells)

``` r
# Explore merged metadata
View(merged_seurat@meta.data)


#Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# adding sample type to metadata. The orginal file could be download from SRA explorer

SampleType <- c("BLCA", "BLCA", "Normal", "BLCA", "BLCA", "BLCA", "BLCA", "BLCA", "BLCA", "Normal", "Normal")
# sample type with grade (Not tested)
#SampleType <- c("BLCA_LG", "BLCA_LG", "Normal", "BLCA_HG", "BLCA_HG", "BLCA_HG", "BLCA_HG", "BLCA_HG", "BLCA_HG", "Normal", "Normal")

names(SampleType) <- c("SRR12603789", "SRR12603790", "SRR12603788", "SRR12603787", "SRR12603786", "SRR12603785", "SRR12603784", "SRR12603783", "SRR12603782", "SRR12603781", "SRR12603780")

metadata$sampleType <- stringr::str_replace_all(metadata$orig.ident, SampleType)


# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA,
                sample = sampleType)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="~/scRNA/merged_filtered_seurat.RData")

#
```

#### Cell counts per sample

``` r
# Visualize the number of cell counts per sample
bqcc <- metadata %>% 
  ggplot(aes(x=seq_folder, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells before QC")
```
![cell_count_plot](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/cell_count_per_sample_before_QC.png)

#### UMI per cell
Typically, we expect the UMI counts per cell to be higher than 500, which is the lower limit of the expected range. If the UMI counts range between 500-1000, the data is still usable, but deeper sequencing may have been beneficial for these cells.
``` r
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  facet_wrap(~seq_folder) +
  geom_vline(xintercept = 1000) +
  labs(fill = "Sample")
  
```
![cell_count_plot](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/nUMI.png)
So the cells have way more than 1K UMI! 

#### Genes detected per cell
In scRNA-seq, the number of genes detected per cell is a crucial quality metric that we expect to be similar to the UMI detection, albeit slightly lower. For high-quality data, the proportional histogram of genes detected per cell should show a single peak that represents encapsulated cells. However, if there is a small shoulder or a bimodal distribution to the left of the main peak, this could indicate a few things. It could be due to some failed cells or biologically different cell types, such as quiescent cell populations or less complex cells of interest. For instance, larger cells or different cell types may have higher gene counts.


``` r
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(x=nGene, fill= sample)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  facet_wrap(~seq_folder) +
  geom_vline(xintercept = 500) +
  labs(fill = "Sample")
  
```
![gene_count_plot](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/nGene_before_QC.png)

#### Novelty score

The novelty score, computed as the ratio of nGenes over nUMI, measures the complexity of RNA species in each cell. A low number of genes detected in a cell with many captured transcripts (high nUMI) indicates low complexity or novelty. This could be due to an artifact, contamination, or represent a specific cell type (e.g. red blood cells). A good quality cell typically has a novelty score above 0.80.

``` r
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, fill=sample)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  facet_wrap(~seq_folder) +
  xlab("Novelty Score") +
  geom_vline(xintercept = 0.8)
```
![gene_count_plot](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/novelt_score_before_QC.png)


#### Mitochondrial gene expression detected per cell

High level of expression from mitochondria indicate dying or dead cells. Basically poor quality samples are those that surpass 0.2 mitochondria ratio mark. 

``` r
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>%
  ggplot(aes(x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  scale_x_continuous(labels = function(x) sprintf("%.1f", x)) + 
  theme_classic() +
  facet_wrap(~seq_folder) +
  geom_vline(xintercept = 0.2)

```
![mitoRatio_plot](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/mitoRatio_before_qc.png)


#### Joint filtering: nUMI, nGene and mitoRatio

```r
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs

metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
  	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 1000) +
  	geom_hline(yintercept = 500) +
  	facet_wrap(~seq_folder)
```
![mitoRatio_plot](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/nUMI_nGene_mitoRatio.png)

There are samples that shows high-quality cells ; high nUMI, high nGene, low number of cells with high mitoRatio and also there are some samples that would clearely benfit from filtering, as they have low quality cells. We expect to see that dying cells to show high level of mitoRatio and low nUMI and nGene . 

Basically, it is not uncommon to observe cells with high numbers of UMIs and nGene with, but also high mitoRatio. These cells may be stressed or damaged, but they could also represent a heterogeneous population of cells with distinct metabolic states.

To investigate the potential cause of high mitochondrial expression ratios, it is important to examine the expression of specific mitochondrial genes and compare them to other genes in the cell. If the expression of mitochondrial genes is elevated relative to other genes, this could suggest mitochondrial dysfunction. Additionally, examining the expression of other stress or damage markers, such as heat shock proteins or cell cycle genes, can also provide insight into the health and state of the cell.


### Filtering

#### Cell-level filtering

-nUMI > 1000

-nGene > 500

-log10GenesPerUMI > 0.8

-mitoRatio < 0.2

``` r
# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(merged_seurat, 
                          subset= nUMI >= 1000 & 
                          nGene >= 500 & 
                          log10GenesPerUMI > 0.80 & 
                          mitoRatio < 0.20)
```

#### Gene-level filtering

Keep only genes which are expressed in 100 or more cells (usually thi is 10)

``` r
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 100 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 100

# Only keeping those genes expressed in more than 100 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Create .RData object to load at any time
save(filtered_seurat, file="seurat_filtered.RData")
```

#### Re-assess QC metrics

``` r

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

# to see drop in filtering cells:

met_before <- data.frame(unclass(table(metadata$seq_folder)))
met_before$QCgroup <- "before"
met_before$cell<- rownames(met_before)
names(met_before)[1] <- 'count'

met_after <- data.frame(unclass(table(metadata_clean$seq_folder)))
met_after$QCgroup <- "after"
met_after$cell<- rownames(met_after)
names(met_after)[1] <- 'count'
# count
cell_count <- data.frame(rbind(met_before, met_after))

                                
# visualization :
cell_count %>% ggplot(aes(x=cell, y=count, fill=QCgroup)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_fill_manual(values = c("#808080", "#FFBF00")) +
  xlab("samples") +
  ggtitle("nCells count before and after QC")
```
![cell_countbefore_after_plot](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/before_after_cell_count.png)

```r
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs

metadata_clean %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
  	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 1000) +
  	geom_hline(yintercept = 500) +
  	facet_wrap(~seq_folder)
  
```
![nUMI_nGene_mitoRatio_plot](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/nUMI_nGene_mitoRatio_after.png)

## 5) Clustering

The ultimate goal is to define clusters of cells and identify cell types in the samples. To achieve this, there are several steps:

1-Identify unwanted variability by exploring data and covariates such as cell cycle and mitochondrial gene expression.
Both biological source of variation (e.g. effect of cell cycle on transcriptome) and technical source should be explored and account for. In early version of Seurat one needs to normalize data, find variable features and then scale data while setting a variable like mitochondrial contamination or cell cycle stage to be regressed out. So here is code for doing these steps to mitigate cell cycle stage effects on the dataset, however a newer function in Seurat has automated all of these steps (`SCTtansform()`).

2-Normalize and remove unwanted variability using Seurat's `sctransform` function. The normalization step is necessary to make expression counts comparable across genes and/or samples. The counts of mapped reads for each gene is proportional to the expression of RNA (“interesting”) in addition to many other factors (“uninteresting” such as sequencing depth and gene length). Normalization is the process of adjusting raw count values to account for the “uninteresting” factors.
For simplicity , normalization is assumed as two step process: scaling and transforming.
In scaling the goal is to multiply each UMI count by a cell specific factor to get all cells to have the same UMI counts.For transformation simple approaches like log-transformation showed to be not that useful, especially in the case of genes with high expression but showing decent performance for low/intreemediate expressed genes. So we cannot treat all genes the same.
The proposed solution for data transformation is Pearson residuals (inmplemented in Seurat's `SCTransform` function), which applies a gene-specific weight to each measurement based on the evidence of non-uniform expression across cells. This weight is higher for genes expressed in a smaller fraction of cells, making it useful for detecting rare cell populations. The weight takes into account not just the expression level but also the distribution of expression.


3- Integrate data using Seurat's method to compare celltype expression between groups.

4-Cluster cells based on similarity of gene expression profiles using Seurat's PCA scores.

5-Evaluate cluster quality by checking for sources of uninteresting variation, principal component influence, and exploring cell type identities using known markers.


### Exploring sources of unwanted variation

#### Evaluating effects of cell cycle and mitochondrial expression

We will score the cells for cell cycle genes, and then determine whether cell cycle is a major source of variation in our dataset using PCA.

``` r
# Normalize the counts
# This normalization method is solely for the purpose of exploring the sources of variation in our data.
seurat_phase <- NormalizeData(filtered_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Load cell cycle markers
load("C:/Users/qaedi/OneDrive - Queen's University/Documents/scRNA/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
#View(seurat_phase@meta.data) 
table(seurat_phase$Phase)
```
Cells in different cell cycle stages:

   G1   G2M     S 
54106 12680 25441 

So most of the cells are in G1 and then S, which make sense.

```r
# Identify the most variable genes and scaling them
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = TRUE)
                     
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_phase), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_phase)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```
![variable_feature](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/variable_features.png)


```r
# Check quartile values for mitoRatio, we will use this variable later to mitigate unwanted source of variation in dataset
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                   breaks=c(-Inf, 0.015, 0.025, 0.045, Inf), 
                   labels=c("Low","Medium","Medium high", "High"))



# Scale the counts
# This step is essential for PCA , clustering and heatmap generation
seurat_phase <- ScaleData(seurat_phase)
#saveRDS(seurat_phase, "seurat_phase.rds")

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)
```
PC_ 1 
Positive:  IGFBP7, A2M, CALD1, SPARC, MGP, SPARCL1, COL4A2, VIM, RBMS3, NNMT 
	   TCF4, COL4A1, COL6A2, LHFPL6, IGFBP4, GNG11, FSTL1, SERPING1, GSN, PALM2-AKAP2 
	   PTPRG, FBXL7, COL6A1, COL1A2, CAV1, ADAMTS9, PMP22, DLC1, LDB2, HTRA1 
Negative:  KRT7, PSCA, DHRS2, C10orf99, CSTB, KRT17, RPS19, S100A14, AQP3, LYPD3 
	   GDF15, LIPH, AKR1C2, FABP5, KRT13, TSPAN1, AC019117.2, UPK1B, SPINK1, SCHLAP1 
	   CRABP2, CXCL17, MACC1, RAB3IP, HPGD, DMKN, UCA1, TMEM97, ADGRF1, FA2H 
PC_ 2 
Positive:  PCAT19, PECAM1, PLVAP, CALCRL, ADGRL4, AQP1, EMCN, ECSCR, RAMP2, MCTP1 
	   VWF, RAMP3, DIPK2B, FLT1, CLDN5, ZNF385D, CD93, PCDH17, SLCO2A1, ERG 
	   RAPGEF4, EGFL7, CYYR1, ACKR1, FAM167B, RHOJ, LDB2, CLEC14A, DOCK4, RNASE1 
Negative:  COL1A2, COL1A1, DCN, COL3A1, COL6A3, RARRES2, COL6A1, PCOLCE, FBLN1, LUM 
	   SERPINF1, C1S, TPM2, COL6A2, C1R, GPC6, SOD3, COL14A1, MFAP4, MEG3 
	   EMILIN1, MXRA8, PRKG1, CTSK, COL5A2, AEBP1, CRISPLD2, FGF7, ANTXR1, LSAMP 
PC_ 3 
Positive:  CD52, TYROBP, CORO1A, ITGB2, FCER1G, RGS1, AIF1, ALOX5AP, HLA-DPB1, IFI30 
	   CRIP1, GPR183, RGS2, LST1, CCL5, MS4A6A, TRAC, CD69, HLA-DPA1, LYZ 
	   CCL4, HLA-DQA1, C1QA, IGKC, SPI1, C1QB, HLA-DRB1, S100A4, C1QC, FCGR2A 
Negative:  KRT7, PSCA, CSTB, C10orf99, DHRS2, AKR1C3, S100A14, TSPAN1, GDF15, SEMA3C 
	   S100A13, AQP3, AKR1C2, KRT17, FABP5, AKR1C1, SLC14A1, AC019117.2, SCHLAP1, CRABP2 
	   CD24, SPTSSB, RAI14, TMEM141, CCND1, HMGCS2, CTSE, LIPH, SPINK1, DMKN 
PC_ 4 
Positive:  RGS5, FRZB, ACTA2, PPP1R14A, MYH11, COX4I2, NDUFA4L2, MYL9, GJA4, TAGLN 
	   HEYL, PTP4A3, MUSTN1, CDH6, MCAM, TPPP3, MYLK, PLN, SEPTIN4, WFDC1 
	   PDGFRB, PRKG1, ADRA2A, MAP1B, MFGE8, CRIP1, HIGD1B, EDNRA, TPM1, CACNA1C 
Negative:  CFD, IGF1, PTGDS, PID1, SFRP2, TMEM176B, AIF1, CD68, C3, SFRP1 
	   SOD2, PDPN, VCAN, MS4A6A, PTGS2, LYZ, LUM, BASP1, HLA-DRA, NEGR1 
	   BNC2, C1QC, IFI30, CD14, LSAMP, C1QA, MMP2, TYROBP, TNFAIP6, TMEM176A 
PC_ 5 
Positive:  LCN15, CTSE, CRTAC1, LINC01088, PLA2G2A, CA9, SLC7A11-AS1, SLC8A1-AS1, CST1, TOX3 
	   SPINK5, FABP4, REG4, ROBO2, GSTM3, SULT1C2, PTPRR, CXCL10, LINC02672, IVL 
	   PNCK, AC096577.1, APOC1, AC078923.1, TCN1, PLAC8, PPP1R3C, UCA1, SPINK1, ANXA10 
Negative:  KRT13, OLFM4, MUC4, SGMS2, CPA6, TRIM31, LYPD3, AGR3, SGPP2, TFPI2 
	   DSG3, SLPI, CXCL17, LIPH, ADGRF1, MACC1, GALNT5, PCDH7, AC022874.1, HBEGF 
	   ITGA2, LMO7, AREG, PLAT, FGFBP1, ALDH1A3, TM4SF1, AL589693.1, CLDN3, NEBL 

```
# Plot the PCA colored by cell cycle phase
no_split <- DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")
        
with_split <- DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by= "Phase")
        
no_split + with_split

```
![cell_cycle_effect](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/cell_cycle_PCA.png)


For mitochondrial expression:
```r
# Plot the PCA colored by mitochondrial expression
no_split <- DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr")
        
with_split <- DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by= "mitoFr")
        
no_split + with_split
```
![mito_expression_effect](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/images/mitoExpression_effect.png)

Based on the above plots, we can see that cells are scattered regardless of their cell cycle phase and mitochondrial genes expression level. So there is no need to regress out the effect of cell cycle and mitochondrial expression in this dataset.


### SCTransform

This function is useful for automatic normalization and regressing out sources of unwanted variation. This method is more accurate method of normalizing, estimating the variance of the raw filtered data, and identifying the most variable genes. In practice `SCTransform` single command replaces `NormalizeData()`, `ScaleData()`, and `FindVariableFeatures()`. Since we have two group of sample we will run SCTransform on each groups after doing "integration".

### Integration

Integrate or align samples across groups using shared highly variable genes If cells cluster by sample, condition, batch, dataset, or even modality (scRNA, scATAC-seq), this integration step can significantly improve the clustering and the downstream analyses. So we want to integrate normal samples together and BLCA sample together , so downstream analysis would make more sense to do. For integration, we have to keep samples as separate objects and transform them as that is what is required for integration.

``` r
# adjust the limit for allowable object sizes within R
options(future.globals.maxSize = 4000 * 1024^2)

# Split seurat object by group
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

# then normalize by SCTansform
for (i in 1:length(split_seurat)) {
    split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio", "S.Score", "G2M.Score"))
    }
    
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 


# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
                                   
# Find best buddies (using canonical correlation analysis or CCA) - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
                                   
# Check assays in the object:
split_seurat$Normal@assays

                                   
```

### Integration check

After normalization and integration, we can proceed to PCA and UMAP/t-SNE to see effect of integration.

``` r
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated, verbose = TRUE)
# PC_ 1
# Positive:  IGFBP7, MGP, SPARC, VIM, SPARCL1, CALD1, IFITM3, A2M, COL4A1, COL1A2
#            COL4A2, TAGLN, COL6A2, BGN, TCF4, COL3A1, CCL2, MT2A, NNMT, COL1A1
#            MYL9, LGALS1, DCN, SELENOM, CCN1, ADAMTS9, CCN2, TIMP1, IGFBP4, LUM
# Negative:  KRT19, KRT13, CD52, S100P, RPS19, CLDN4, FXYD3, PSCA, TACSTD2, CSTB
#            PTPRC, KRT7, RPS18, KRT17, LYPD3, RPL41, CD3D, CCL5, RPS27, SFN
#            FABP5, SPINK1, GDF15, AQP3, HCST, TRAC, SNCG, ADIRF, ELF3, TRBC2
# PC_ 2
# Positive:  PLVAP, CALCRL, PCAT19, MCTP1, AQP1, VWF, CD74, PECAM1, RAMP2, LDB2
#            RAMP3, FLT1, ZNF385D, TCF4, CLDN5, ACKR1, HSPG2, SPARCL1, EMCN, ADGRL4
#            SLCO2A1, HLA-DRA, SELE, DOCK4, ECSCR, CCL14, HLA-DRB1, ERG, PCDH17, GNG11
# Negative:  COL1A2, COL3A1, COL1A1, TAGLN, DCN, LUM, BGN, COL6A2, C1R, SOD3
#            MFAP4, RARRES2, C1S, MYL9, TPM2, PRKG1, CRYAB, ACTA2, COL6A3, LGALS1
#            CALD1, COL6A1, TIMP1, SERPINF1, PCOLCE, AEBP1, C11orf96, MEG3, FBLN1, GPC6
# PC_ 3
# Positive:  CCL5, CD52, B2M, PTPRC, IL32, SRGN, CD3D, HSPA1A, NKG7, GZMA
#            TRAC, ARHGAP15, HCST, RGS1, CXCR4, SAMSN1, CCL4, CD7, RGS2, CORO1A
#            CD2, CRIP1, CST7, CD69, STAT4, FYN, TMSB4X, DNAJB1, PTPN22, CD3E
# Negative:  ADIRF, SPINK1, IFI27, CSTB, S100P, KRT19, FXYD3, CLDN4, CCT2, CCND1
#            PSCA, S100A6, SNCG, UCA1, KRT17, YEATS4, KRT7, KRT13, TACSTD2, GDF15
#            RAB3IP, KRT18, HES1, GAPDH, S100A2, PLVAP, RAMP2, S100A14, FABP5, TM4SF1
# PC_ 4
# Positive:  IL32, CD3D, CCL5, CRIP1, TRAC, FYN, COL4A1, COL4A2, CD2, GZMA
#            CD7, CALD1, IGFBP7, TRBC2, MCAM, NDUFA4L2, CD3E, NKG7, COL18A1, RGS5
#            SKAP1, MYL9, ITGA1, CYTOR, CD247, SPARC, TRBC1, MAP1B, ACTA2, CACNA1C
# Negative:  HLA-DRA, CD74, TYROBP, HLA-DRB1, HLA-DPB1, FCER1G, HLA-DPA1, AIF1, HLA-DQA1, FTL
#            LYZ, IFI30, C1QA, HLA-DQB1, C1QB, C1QC, MS4A6A, LST1, APOE, CD14
#            TMEM176B, S100A9, HLA-DMA, FCGR2A, CXCL8, SPI1, FTH1, CD68, PSAP, CST3
# PC_ 5
# Positive:  LUM, MMP2, PTGDS, EMP1, DCN, RARRES2, KRT13, SERPINF1, CXCL1, CXCL8
#            LSAMP, COL8A1, TM4SF1, AREG, C1S, PDPN, CFD, APOD, SOD2, CTSK
#            CLMP, KRT17, MFAP4, LYPD3, VCAN, PLAUR, TSHZ2, PLAT, PDGFRA, C1R
# Negative:  RGS5, NDUFA4L2, ACTA2, PPP1R14A, MYL9, TAGLN, FRZB, CRIP1, CALD1, GJA4
#            MCAM, COL18A1, MYH11, TYROBP, TPPP3, COX4I2, PRKG1, IGFBP7, COL4A1, COL4A2
#            CDH6, HIGD1B, AIF1, TPM2, PTP4A3, FCER1G, WFDC1, HEYL, MYLK, HLA-DRA



# Plot PCA
png(filename = "PCA_integrated.png", width = 16, height = 8.135, units = "in", res = 300)
PCAPlot(seurat_integrated,
        split.by = "sample")
dev.off()

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca",
                             verbose = TRUE)

# Plot UMAP 
png(filename = "UMAP_integrated.png", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(seurat_integrated, split.by = "sample")
dev.off()
```

### Clustering cells based on top PCs (metagenes)

#### Identify significant PCs

For new method like SCTransform it is not needed to calculate the number of PCs for clustering. However older methods could not efficiently removed technical biases , so using them it was necessary to have some idea about the number of PCs that can capture most of information in the dataset.

``` r
# Explore heatmap of PCs
png(filename = "PCA_integrated_2.png", width = 16, height = 8.135, units = "in", res = 300)
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
dev.off()

# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
      
# To determine how many Pcs should be considered for clustering:
# Plot the elbow plot
png(filename = "elbow.png", width = 16, height = 8.135, units = "in", res = 300)
ElbowPlot(object = seurat_integrated, 
          ndims = 40)
dev.off()

# to make it more quantitative :
# Determine percent of variation associated with each PC
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation is the optimal number of PC to pick.
pcs <- min(co1, co2)

pcs
```

#### Cluster the cells

``` r
# to check what is active assay
DefaultAssay(object = seurat_integrated)

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                dims = 1:18)
                                
Find clusters
# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
                               

# Explore resolutions
head(seurat_integrated@meta.data)

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
png(filename = "umap_cluster_with_label.png", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()
```

#### Clustering quality control

After clustering, we need to make sure that the assigned clusters are true representative of biological clusters (cell clusters) not due to technical or unwanted source of variation (like cell cycle stages). Also , in this step we need to identify cell type for each cluster based on the known cell type markers.

-   Segregation of clusters by sample

``` r
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample

library(dplyr)
library(tidyr)

# n_cells <- FetchData(seurat_integrated, 
#                      vars = c("ident", "orig.ident")) %>%
#         dplyr::count(ident, orig.ident) %>%
#         tidyr::spread(ident, n)

n_cells <- FetchData(seurat_integrated, 
                      vars = c("ident", "orig.ident"))
n_cells <- dplyr::count(n_cells, ident, orig.ident)
n_cells <- tidyr::spread(n_cells, ident, n)



#Ading sample data from paper; we expect to see samples from same group have more or less similar number of cells in each cluster. 
#So normal samples should show similar patters: SRR12603780, SRR12603781, and SRR12603788.



sampleData<- data.frame(tibble::tribble(
     ~sample_id, ~gender, ~age, ~Grade, ~Invasiveness, ~Surgery_Type, ~Tumor_size_cm,
  "SRR12603790",     "M",  67L,  "low", "Noninvasive",       "TURBT",          "1.9",
  "SRR12603789",     "M",  70L,  "low", "Noninvasive",       "TURBT",          "2.5",
  "SRR12603787",     "M",  63L, "high", "Noninvasive",  "Cystectomy",          "3.5",
  "SRR12603786",     "F",  59L, "high", "Noninvasive",  "Cystectomy",          "4.7",
  "SRR12603785",     "M",  57L, "high",    "Invasive",  "Cystectomy",          "5.1",
  "SRR12603784",     "M",  75L, "high",    "Invasive",  "Cystectomy",          "4.3",
  "SRR12603783",     "M",  77L, "high",    "Invasive",  "Cystectomy",          "4.5",
  "SRR12603782",     "F",  72L, "high",    "Invasive",  "Cystectomy",          "4.1",
  "SRR12603781",     "M",  67L, "normal",   "normal",       "TURBT",            "-",
  "SRR12603788",     "M",  75L, "normal",   "normal",  "Cystectomy",            "-",
  "SRR12603780",     "M",  63L, "normal",   "normal",  "Cystectomy",            "-"
  ))


# View table
head(n_cells)
# saving objects (to mark where and when we stored the file)
#seurat_cluster <- seurat_integrated
#saveRDS(seurat_cluster, "seurat_cluster.RDS")


# UMAP of cells in each cluster by sample
# This would allow us to see condition specefic clusters
png(filename = "umap_cluster_sample.png", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
dev.off()





# Segregation of clusters by cell cycle phase (unwanted source of variation) 
# Explore whether clusters segregate by cell cycle phase

png(filename = "umap_cluster_cell_cucle.png", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
dev.off()
```

-   Segregation of clusters by various sources of uninteresting variation

We expect to see a uniform coluring for all variables in all clusters. Sometimes this is not the case. Like here `nUMI` and `nGene` showing higher value is some clusters. We have to watch these cluster and inspect them in terms of type of cell therein. So that may explain some of the variation that we are seeing.

``` r
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
png(filename = "umap_unwanted_source_clustering.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
dev.off()
```

-   Exploration of the PCs driving the different clusters

We hope that the defined PCs could separate clusters well.We can see how the clusters are represented by the different PCs.Then we could look back at our genes driving this PC to get an idea of what the cell types might be in each cluster.

``` r
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:18),
            "ident",
            "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
                     
                     
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))



# Plotting a UMAP plot for each of the PCs
library(cowplot)
library(tidyverse)
library(HGNChelper)

png(filename = "umap_on_pcs.png", width = 16, height = 8.135, units = "in", res = 300)
map(paste0("PC_", 1:18), function(pc){
        ggplot(pc_data, 
               aes(UMAP_1, UMAP_2)) +
                geom_point(aes_string(color=pc), 
                           alpha = 0.7) +
                scale_color_gradient(guide = FALSE, 
                                     low = "grey90", 
                                     high = "blue")  +
                geom_text(data=umap_label, 
                          aes(label=ident, x, y)) +
                ggtitle(pc)
}) %>% 
        plot_grid(plotlist = .)
dev.off()


# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)
```

``` r
# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

png(filename = "umap_fibroblast.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("IGFBP7", "MGP"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()

png(filename = "umap_endothelial.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("PLVAP", "CALCRL"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()


png(filename = "umap_t_cells.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CCL5", "CD52", "IL32"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
```

### Marker identification

This can be the last step in our pipeline which aims to determine the gene markers for each of the clusters and identify cell types of each cluster using markers. Also this step helps to determine whether there's a need to re-cluster based on cell type markers, or maybe clusters need to be merged or split.

-   Identification of all markers for each cluster

``` r
# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = seurat_integrated, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)     
DefaultAssay(seurat_integrated) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                              ident.1 = 0,
                              grouping.var = "sample",
                                only.pos = TRUE,
                             logfc.threshold = 0.60)
```

To add more annotation to the results

``` r
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
        mcols() %>%
        rownames() %>%
        tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
        dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
```

To find conserved markers for all clusters:

``` r
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
               by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
  }
  
# this function can be an argument for 'map_dfr' function :
# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0:21), get_conserved)

# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (Normal_avg_log2FC + BLCA_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)
        
head(top10)

data.table::fwrite(top10, "top10_conserved_markers.csv")
```

There are a number of tools that one may use to assign cell type to a cluster. However almost non of them at the time of writing can help with bladder tissue. So I have to looked up for genes in the [PanglaoDB](https://panglaodb.se/index.html) database manually and assign cell type to each cluster. So to do this job in a more efficient way, lets first identify which markers are associated to more clusters, then assign cell type to those clusters.

``` r
top10_mod <- data.frame(unclass(table(top10$gene, top10$cluster_id)))


data.table::fwrite(top10_mod, "top10_mod_conserved_markers.csv")

# markers with highest frequency
M <- c("KRT7", "KRT19", "FCER1G", "AIF1", "AQP3", "CCL5", "CD24", "CD3D", "CD52", "CLDN4", "COL1A1", "COL1A2", "CRTAC1", "CXCL8", "DCN", "FABP4", "FABP5", "FXYD3", "GZMA", "HLA-DRA", "IGHA1", "IGHG1", "IGHG3", "IGHG4", "IGKC", "IGLC1", "IGLC2", "IGLC3", "JCHAIN")
```
So cell type for top markers:

| cell type                | genes                                                 |
|--------------------------|-------------------------------------------------------|
| Basal cells              | KRT7, KRT19, AQP3, CD24,CXCL8,FXYD3                   |
| Dendritic cells          | FCER1G,AIF1,FABP4                                     |
| Gamma delta T cells      | CCL5,GZMA                                             |
| NK cell                  | CCL5,GZMA                                             |
| T cells                  | CD3D, CD52                                            |
| Macrophages              | CD52                                                  |
| Luminal epithelial cells | CLDN4                                                 |
| Fibroblasts              | COLA1, COLA2,CXCL8, DCN                               |
| Epithelial cells         | CRTAC1,FXYD3                                          |
| Endothelial cells        | FABP4, FABP5                                          |
| Plasma cell              | IGHA1,IGHG1,IGHG3,IGHG4,IGKC,IGLC1,IGLC3,IGLC3,JCHAIN |

Lets visualize some of the genes and see in which cluster they show expression.

```R

# Plot interesting marker gene expression 
png(filename = "umap_high_freq_basal_cells.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
                        features = c("KRT7", "KRT19", "AQP3", "CD24", "FXYD3", "CXCL8"),
                         sort.cell = TRUE,
                         min.cutoff = 'q10', 
                         label = TRUE,
                         repel = TRUE)
dev.off()


# Vln plot - cluster 0
png(filename = "violin_high_freq_basal_cells.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = c("KRT7", "KRT19", "AQP3", "CD24", "FXYD3", "CXCL8"))
dev.off() 
```
So according to "basal cells" visualization , following clusters may show clusters of basal cells:
0, 2,3,9,10,14,17,18,19(?),20 and 21.


``` r
# Plot interesting marker gene expression 
png(filename = "umap_high_freq_Plasma_cells.png", width = 26, height = 15.135, units = "in", res = 600)
FeaturePlot(object = seurat_integrated, 
                        features = c("IGHA1","IGHG1","IGHG3","IGHG4","IGKC","IGLC1","IGLC3","IGLC3","JCHAIN"),
                         sort.cell = TRUE,
                         min.cutoff = 'q10', 
                         label = TRUE,
                         repel = TRUE)
dev.off()


# Vln plot - cluster 0
png(filename = "violin_high_freq_basal_cells.png", width = 26, height = 10.135, units = "in", res = 600)
VlnPlot(object = seurat_integrated, 
        features = c("IGHA1","IGHG1","IGHG3","IGHG4","IGKC","IGLC1","IGLC3","IGLC3","JCHAIN"))
dev.off() 
```
Cluster 19 looks to be hard to assign it a cell type. Will keep an eye on it.

```R

# Plot interesting marker gene expression 
png(filename = "umap_high_freq_dc.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
                        features = c("FCER1G","AIF1","FABP4"),
                         sort.cell = TRUE,
                         min.cutoff = 'q10', 
                         label = TRUE,
                         repel = TRUE)
dev.off()


# Vln plot - cluster 0
png(filename = "violin_high_freq_dc.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = c("FCER1G","AIF1","FABP4"))
dev.off() 
```

For cluster1 and 4

```R
d <- data.table::fread("top10_conserved_markers.csv")

png(filename = "umap_cluster4_markers.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
                        features = d$gene[d$cluster_id == "4"],
                         sort.cell = TRUE,
                         min.cutoff = 'q10', 
                         label = TRUE,
                         repel = TRUE)
dev.off()


# Vln plot - cluster 0
png(filename = "violin_cluster4_markers.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = d$gene[d$cluster_id == "4"])
dev.off() 

png(filename = "umap_cluster12_markers.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
                        features = d$gene[d$cluster_id == "12"],
                         sort.cell = TRUE,
                         min.cutoff = 'q10', 
                         label = TRUE,
                         repel = TRUE)
dev.off()

```
So it seems clusters 1,4,8,12,13 and 16 are the same cell type with minor subtypes.However cluster 5 and 7 showing some level of expression for markers in cluster 1 + 4, but they need separate inspection. 
Conserved markers in these two clusters are as follow

```
[1] "LDB2"    "SPARCL1" "GNG11"   "FLT1"    "IFI27"   "RAMP2"   "PECAM1"
[8] "TCF4"    "PLVAP"   "PCAT19"  "ACKR1"   "SELE"    "CALCRL"  "ZNF385D"
[15] "ADAMTS9" "MCTP1"   "AQP1"    "VWF"     "CCL14"   "PCAT19"
```

| cell type                | genes                                                 |
|--------------------------|-------------------------------------------------------|
| Endothelial cells        | LDB2, SPARCL1, FLT1,IFI27,RAMP2,PECAM1,TCF4,PLVAP,PCAT19,ACKR1, SELE,CALCRL,ZNF385D,ADAMTS9,MCTP1,AQP1,VWF,                     |
| Fibroblast              | SPARCL1,TCF4,PLVAP,ADAMTS9                                |



```R
png(filename = "umap_cluster5_7_markers.png", width = 26, height = 18.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
                        features = d$gene[d$cluster_id %in% c("5", "7")],
                         sort.cell = TRUE,
                         min.cutoff = 'q10', 
                         label = TRUE,
                         repel = TRUE)
dev.off()


# Vln plot - cluster 0
png(filename = "violin_cluster5_7_markers.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = d$gene[d$cluster_id == "4"])
dev.off()
```
So as expected cluster 7 and 5 representing same cell types. 
In almost all cases cluster 17 showed some level of expression for EC markers. 

Lets now have a look at markers for clusters 6,15 and 11:

```R
png(filename = "umap_cluster11_15_6_markers.png", width = 26, height = 18.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
                        features = unique(d$gene[d$cluster_id %in% c("6", "11","15")]),
                         sort.cell = TRUE,
                         min.cutoff = 'q10', 
                         label = TRUE,
                         repel = TRUE)
dev.off()
```
According to the markers such ass different collagen, clusters 11 and 15 seem to be fibroblasts while cluster 6 showing a mix of fibroblast and non-fibroblast markers.

Now we can generate UMAP with cell type as labels;

```R
# Rename all identities
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                               "0" = "Basal cells",
                               "1" = "Gamma delta T cell", # NKG7 and GZMK
                               "2" = "Basal cells",
                               "3" = "Basal cells",
                               "4" = "T cells",
                               "5" = "Endothelial cells",
                               "6" = "Fibroblast",
                               "7" = "Endothelial cells",
                               "8" = "DC",
                               "9" = "Basal cells",
                               "10" = "Basal cells",
                               "11" = "Fibroblast",
                               "12" = "T cells",
                               "13" = "T cells",
                               "14" = "Basal cells",
                               "15" = "Fibroblast",
                               "16" = "T cells",
                               "17" = "?", 
                               "18" = "Basal cells", 
                               "19" = "Plasma cells",
                               "20" = "Basal cells",
                               "21" = "Basal cells")

# Plot the UMAP withy new labells
png(filename = "umap_with_label.png", width = 16, height = 8.135, units = "in", res = 600)
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
dev.off()
```
OK then it's time to write the Seurat object for later analysis , if any:

```R
# Save final R object
write_rds(seurat_integrated,
          path = "seurat_labelled.rds")      
```
As final note I tried to assign cell types to clusters using scType tool, however there was no bladder tissue profile in their database. So the generated lables using this apporach should be only considered as a guide.

```R
# Cell maker assignemnt using scType

# load libraries
library(HGNChelper)


# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
#tissue = c("Immune system") # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
tissue = c("Immune system", "Pancreas", "Liver","Kidney","Intestine","Placenta","Spleen",
           "Stomach") 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = seurat_integrated[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call(
  "rbind", lapply(unique(seurat_integrated@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_integrated@meta.data[seurat_integrated@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_integrated@meta.data$seurat_clusters==cl)), 10)
}))


sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])



# overlay the identified cell types on UMAP plot:

seurat_integrated@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_integrated@meta.data$customclassif[seurat_integrated@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

png(filename = "umap_with_label_scType.png", width = 16, height = 8.135, units = "in", res = 600)
DimPlot(seurat_integrated, 
        reduction = "umap", 
        label = TRUE, 
        repel = TRUE, 
        group.by = 'customclassif')        

dev.off()

```
