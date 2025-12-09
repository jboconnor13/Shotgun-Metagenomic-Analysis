#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --job-name=env_creation
#SBATCH --nodes=1 # use 1 node 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=16
#SBATCH --time=00:10:00 # Time limit days-hrs:min:sec
#SBATCH --qos=normal
#SBATCH --mem=10gb # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=john.2.oconnor@cuanschutz.edu

module purge 
module load miniforge
module load python/3.10.2

#The environment is created
conda create -n shotgun_analysis

#The environment is activated
conda activate shotgun_analysis

#Snakemake is installed
conda install -c conda-forge snakemake

#The environment proiorities are set to strict for snakemake
conda config --set channel_priority strict