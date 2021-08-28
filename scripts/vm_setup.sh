#!/bin/bash

########################################
# vm_setup.sh
#
# simple bash script to deploy a
# reproducable project structure
########################################

# move to home directory and create directories
cd ~
mkdir data github notebooks ref results scripts software .ARCHIVE
mkdir data/raw data/preprocessed
mkdir notebooks/exploration notebooks/publication
mkdir ref/hg38

# this part is optional but useful for me
cd ~/software
git clone https://github.com/mvinyard/vintools.git
git clone https://github.com/mvinyard/vintools.wiki.git

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
rm Anaconda3-2021.05-Linux-x86_64.sh # updated June 29, 2021
mv vm_setup.sh scripts


# -----------------------------
# Swap Space
# -----------------------------
sudo fallocate -l4G /swapfile
sudo chmod 700 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
echo "/swapfile none swap sw 0 0" | sudo tee -a /etc/fstab
