#!/bin/bash

### run using the following command:
# bash star_gg.sh <genomeDir> <genomeFA> <GTF> <sjdbOverhang> <THREADS>

### outline of parsed args
#1=<genomeDir> 
#1=<genomeFastaFiles> / genomeFA
#2=<sjdbGTFfile> / GTF
#3=<sjdbOverhang> / splice junction overhang (N-1 of the read length)
#4=<THREADS>

# required parameters:
genomeDir=${1}
genomeFA=${2}
GTF=${3}

overhang=${4}
# splice junction overhang

# optional parameters:
THREADS=${5:-16}

STAR --runThreadN ${THREADS} \
--runMode genomeGenerate \
--genomeDir ${genomeDir} \
--genomeFastaFiles ${genomeFA} \
--sjdbGTFfile ${GTF} \
--sjdbOverhang ${overhang} #149 for a 150 bp PE reads (N-1)
