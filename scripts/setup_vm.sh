#!/bin/bash 

# sequence read analysis tools 
conda install -c bioconda samtools -y
conda install -c bioconda bedtools -y
conda install -c bioconda star -y

# python essentials
conda install -c anaconda numpy -y
conda install -c anaconda pandas -y
conda install -c anaconda matplotlib -y

# packages needed for jupyter notebooks
conda install -c conda-forge port-for -y
conda install -c conda-forge nb_conda_kernels -y

# cleanup
rm Anaconda3-2020.11-Linux-x86_64.sh 
