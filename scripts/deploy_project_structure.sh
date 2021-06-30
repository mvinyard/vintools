#!/bin/bash

########################################
# deploy_project_structure.sh 
#
# simple bash script to deploy a
# reproducable project structure
########################################

cd ~

mkdir data github notebooks ref results scripts software .ARCHIVE

mkdir data/raw data/preprocessed

mkdir notebooks/exploration notebooks/publication

mkdir ref/hg38

# this part is optional

cd ~/software
git clone https://github.com/mvinyard/vintools.git
git clone https://github.com/mvinyard/vintools.wiki.git

# clean up
mv deploy_project_structure.sh scripts
