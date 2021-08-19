#!/bin/bash

### run using the following command:
# bash STAR_pw.sh <SRA_accession_list> <data_directory> <alignment_out_directory> <threads> <outFilterMatchNmin_val> <outFilterScoreMinOverLread_val> <outFilterMatchNminOverLread_val>

### outline of parsed args
#1=<SRA_accession_list>
#2=<DATADIR> / data_directory
#3=<OUTDIR> / alignment_out_directory
#4=<THREADS>
#5=<outFilterMatchNmin_val>
#6=<outFilterScoreMinOverLread_val>
#7=<outFilterMatchNminOverLread_val>

# required parameters:
SRR_FILE=${1}
DATADIR${2}
OUTDIR=${3}

# optional parameters:
THREADS=${4:-16}
outFilterMatchNmin_val=${5:-37} # min read length after soft-clipping
outFilterScoreMinOverLread_val=${6:-0}
outFilterMatchNminOverLread_val=${7:-0}

# these presets are from the preprocessing of the scRNA-seq data from:
# A single-cell RNA-seq survey of the developmental landscape of the human prefrontal cortex
# Zhong et al., Nature volume 555, pages524–528(2018)
# FASTQ files downloaded from GEO accession: GSE104276

SRR_LIST=$(cat ${SRR_FILE})

for SRR in ${SRR_LIST}
do
  echo ${SRR}
  OUT=${OUTDIR}/${SRR}
  mkdir ${OUT}

  READ1=${DATADIR}/${SRR}/${SRR}_1.fastq
  READ2=${DATADIR}/${SRR}/${SRR}_2.fastq

  STAR --runThreadN ${THREADS} \
  --genomeDir /home/mvinyard/ref/ \
  --readFilesIn ${READ1} ${READ2} \
  --outFileNamePrefix ${OUT}/aligned_${SRR}_ \
  --outFilterScoreMinOverLread ${outFilterScoreMinOverLread_val} \
  --outFilterMatchNminOverLread ${outFilterMatchNminOverLread_val} \
  --outFilterMatchNmin ${outFilterMatchNmin_val} \
  --twopassMode Basic
done

