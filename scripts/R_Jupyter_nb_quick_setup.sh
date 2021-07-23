#!/bin/bash
# R_Jupyter_nb_quick_setup.sh

env_name=$1

echo "Creating environment: $env_name"
conda create -n ${env_name} -y

source /home/mvinyard/.anaconda3/etc/profile.d/conda.sh

conda activate ${env_name}

echo "Now installing notebook..."
conda install -c conda-forge notebook --quiet -y

echo "Now installing r-base..."
conda install -c conda-forge r-base=4.1.0 --quiet -y

echo "Now installing python package port_for using pip..."
pip install port_for

sudo apt-get update
sudo apt-get install libzmq3-dev libcurl4-openssl-dev libssl-dev jupyter-core jupyter-client

Rscript='Renv_setup_Rcommands.R'

R_mirror="options('repos' = c(CRAN = 'http://cran.ma.imperial.ac.uk/'))"
R_commands1="install.packages(c('repr', 'IRdisplay', 'IRkernel'), type = 'source')"
R_commands2="IRkernel::installspec()"
echo $R_mirror >> $Rscript
echo $R_commands1 >> $Rscript
echo $R_commands2 >> $Rscript

chmod a+x $Rscript
Rscript $Rscript

echo "Now installing r-irkernel..."
conda install -n ${env_name} r-irkernel --quiet -y
echo "Now installing nb_conda_kernels..."
conda install nb_conda_kernels --quiet -y

echo "Exporting YAML file to ~/.anaconda3/${env_name}.yml for future use..."
conda env export --name ${env_name} > ~/.anaconda3/${env_name}.yml

# clean up
rm Rscript

echo "Done!"
