# scRNA sequencing analysis

Steps toward making sense of single cell RNA sequencing data.

#### Update;
The repo has been update by the new data set coming from Zhaohui Chen *et al.* [paper](https://www.nature.com/articles/s41467-020-18916-5#Sec12) published in 11, Article number: 5077 (2020). Sampels ;

" Eight primary bladder tumor tissues (2 low-grade bladder urothelial tumors, 6 high-grade bladder urothelial tumors) along with 3 adjacent normal mucosae, were involved in this cohort." 


### 1) Data download on ComputeCanada 

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
The next step is to run `cellranger count` on each samples. Running `cellranger` on a cluster with `SLURM` as job schaduler is not an easy task. The following is what I came up with working best for me working on ComputeCanada cluster:

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
