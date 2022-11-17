# scRNA sequencing analysis

Steps toward making sense of single cell RNA sequencing data.

#### Update;
The repo has been update by the new data set coming from Zhaohui Chen *et al.* [paper](https://www.nature.com/articles/s41467-020-18916-5#Sec12) published in 11, Article number: 5077 (2020). Sampels ;

" Eight primary bladder tumor tissues (2 low-grade bladder urothelial tumors, 6 high-grade bladder urothelial tumors) along with 3 adjacent normal mucosae, were involved in this cohort." 


## 1) Data download on ComputeCanada 

```shell

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

```shell
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
There are diffrences between sequence length of read 1 and read 2 for each samples; read 1 provides data on Cell barcode & UMI and read 2 has the insert sequence. Files coming from SRA usually come with names like something_1_fastq.gz and something_2_fastq.gz . This [blog](https://kb.10xgenomics.com/hc/en-us/articles/115003802691) is  helpful to see what are the naming requirements of fastq files for cellranger tool. Briefly; 

- incompatible file name: SRR9291388_1.fastq.gz

- compatible file name: SRR9291388_S1_L001_R1_001.fastq.g

So we need to change file names and also set up the directories for `cellranger count`;

```shell
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

```shell
#!/bin/bash
#SBATCH --account=def-gooding-ab
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
#SBATCH --mail-user=qaedi.65@gmail.com


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

The next steps are mainly based on  Harvard Chan Bioinformatics Core materials on [scRNA-seq analysis](https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html) training. 

## 3) Loading single-cell RNA-seq count data

```R
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
                       y = c(SRR12603782,
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

- orig.ident: this column will contain the sample identity if known. It will default to the value we provided for the project argument when loading in the data

- nCount_RNA: represents the number of UMIs per cell. UMI (unique molecular identifiers) is used to determine whether a read is a biological or technical duplicate (PCR duplicate). Reads with different UMIs mapping to the same transcript were derived from different molecules and are biological duplicates - each read should be counted. Reads with the same UMI originated from the same molecule and are technical duplicates - the UMIs should be collapsed to be counted as a single read.


- nFeature_RNA: represents the number of genes detected per cell

Recommended features to add to metadata:

- number of genes detected per UMI (or novelty score): more genes detected per UMI, more complex our data

- mitochondrial ratio: this metric will give us a percentage of cell reads originating from the mitochondrial genes (coming from dying cells)

```R
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

```R
# Visualize the number of cell counts per sample
bqcc <- metadata %>% 
  ggplot(aes(x=seq_folder, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells before QC")

```
#### UMI per cell

```R
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=seq_folder, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
  
```
#### Genes detected per cell

```R
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=seq_folder, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
```
#### Novelty score

```R
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
```

#### Mitochondrial gene expression detected per cell

```R
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=seq_folder, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
```
#### Checking for cells with low numbers of genes/UMIs

```R
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
```
### Filtering

#### Cell-level filtering

-nUMI > 500

-nGene > 250

-log10GenesPerUMI > 0.8

-mitoRatio < 0.2

```R
# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(merged_seurat, 
                          subset= nUMI >= 500 & 
                          nGene >= 250 & 
                          log10GenesPerUMI > 0.80 & 
                          mitoRatio < 0.20)

```

#### Gene-level filtering

Keep only genes which are expressed in 10 or more cells

```R
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Create .RData object to load at any time
save(filtered_seurat, file="seurat_filtered.RData")
```

#### Re-assess QC metrics

```R

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

# Visualize the number of cell counts per sample
aqcc <- metadata_clean %>% 
  ggplot(aes(x=seq_folder, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells after QC")

# cell count before and after QC
bqcc + aqcc


# Visualize the number UMIs/transcripts per cell
metadata_clean %>% 
  ggplot(aes(color=seq_folder, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata_clean %>% 
  ggplot(aes(color=seq_folder, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata_clean %>% 
  ggplot(aes(color=seq_folder, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
  
```
## 5) Normalization and regressing out unwanted variation

### Explore sources of unwanted variation

Both biological source of variation (e.g. effect of cell cycle on transcriptome) and technical source should be explored and account for.

```R
# Normalize the counts
# This normalization method is solely for the purpose of exploring the sources of variation in our data.
seurat_phase <- NormalizeData(filtered_seurat)
```
### Evaluating effects of cell cycle

```R
# Load cell cycle markers
load("data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)    
```


