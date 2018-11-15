#!/bin/bash

### Job name
#PBS -N HGTID

### Join queuing system output and error files into a single output file
#PBS -j oe

### Send email to user when job ends or aborts
#PBS -m ae

### email address for user
#PBS -M m.gauthier@garvan.org.au

### Queue name that job is submitted to
#PBS -q normal

### Provide project to be billed
#PBS -P jp48

### Make the data readable to everyone in the group
#PBS -W umask=027

### Set out (Use job ID so it's unique)
#PBS -e /g/data1a/jp48/scripts/transposon/ocscc/MakeVCF_LINE1.${ERRFILE}_${TASK_ID}

### Request nodes, memory, walltime. NB THESE ARE REQUIRED.
##### run, change these entries
#PBS -l ncpus=4
#PBS -l mem=12gb
#PBS -l walltime=24:00:00

# Load required modules
module load java
#module load bowtie2
#module load samtools
#module load java

echo Running on host `hostname`
echo Time is `date`

# This jobs working directory

#output folder should not already exist
echo Working directory is ${PBS_O_WORKDIR}
cd ${PBS_O_WORKDIR}
#sample=183410_T
#sample=193958_T


#perl /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/src/hgt_step2.pl -b /g/data1a/jp48/trafic/${sample}_recal_reads_sub.bam -c /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/config.txt -o /g/data1a/jp48/scripts/hgtid/${sample}

#perl /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/src/hgt_step3_onwards.pl -b /g/data1a/jp48/trafic/${sample}_recal_reads_sub.bam -c /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/config.txt -o /g/data1a/jp48/scripts/hgtid/${sample}

##bash /g/data1a/jp48/scripts/hgtid/193958_T/.scripts/extract.ed.h.sh

#/g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/bin/samtools/samtools index /g/data1a/jp48/gdc-client.1.3/bams/https\:/api.gdc.cancer.gov/data/53f0ce53-95fe-4d6c-930d-f8427c0c46ec/TCGA-CV-5442-01A-01D-2266-10_Illumina.bam

perl /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/src/hgt_adapted2.pl -b /g/data1a/jp48/gdc-client.1.3/bams/https:/api.gdc.cancer.gov/data/53f0ce53-95fe-4d6c-930d-f8427c0c46ec/TCGA-CV-5442-01A-01D-2266-10_Illumina.bam -s TCGA-CV-5442 -c  /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/config.txt -o /g/data1a/jp48/scripts/hgtid/TCGA-CV-5442
#perl /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/src/hgt_adapted.pl -b /g/data1a/jp48/young_cohort/TCGA-BA-5557/bb531fb4-064f-47dd-8231-bd6131866968/TCGA-BA-5557-01A-01D-1509_120420_SN1120_0135_BD0T5FACXX_s_4_rg.sorted.bam -s TCGA-BA-5557 -c /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/config.txt -o /g/data1a/jp48/scripts/hgtid/TCGA-BA-5557
#'@RG\tID:foo\tSM:bar'
#SAMPLE=TCGA-BA-5557
#/g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0_original/bin/bwa/bwa mem -t 2 -M -R '@RG\tID:TCGA-BA-5557\tSM:TCGA-BA-5557' /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0_original/resources/virus.fa /g/data1a/jp48/scripts/hgtid/TCGA-BA-5557/.human_again/read1.fq /g/data1a/jp48/scripts/hgtid/TCGA-BA-5557/.human_again/read2.fq | /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0_original/bin/samtools/samtools view -u -bS - >  /g/data1a/jp48/scripts/hgtid/TCGA-BA-5557/.virus/virus.ed.bam

#/g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/bin/samtools/samtools view -b -f2 -F 256 -F 1024 -L /g/data1a/jp48/scripts/hgtid/TCGA-CV-5442/regions.2keep.bed /g/data1a/jp48/scripts/hgtid/TCGA-CV-5442/.human_again/2_human.bam > /g/data1a/jp48/scripts/hgtid/TCGA-CV-5442/human.again.bam
#/g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/bin/samtools/samtools merge -f /g/data1a/jp48/scripts/hgtid/TCGA-CV-5442/merged.bam /g/data1a/jp48/scripts/hgtid/TCGA-CV-5442/VIRUS.HUMAN.flt.sort.bam /g/data1a/jp48/scripts/hgtid/TCGA-CV-5442/human.org.bam /g/data1a/jp48/scripts/hgtid/TCGA-CV-5442/human.again.bam /g/data1a/jp48/scripts/hgtid/TCGA-CV-5442/virus.org.bam
