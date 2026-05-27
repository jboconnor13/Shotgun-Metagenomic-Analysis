#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --output=/scratch/alpine/joconnor@xsede.org/Shotgun-Metagenomic-Analysis/slurm_outputs/slurm-%j.out
#SBATCH --job-name=metaphlan_db_download
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=05:00:00
#SBATCH --qos=normal
#SBATCH --mem=150gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=john.2.oconnor@cuanschutz.edu

module purge 
module load miniforge
module load python/3.10.2

#A MetaPhlan Conda Environment is Created
conda create --name MetaPhlAn

#The environment is activated
conda activate MetaPhlAn

#Metaplan is installed into the environment
conda install metaphlan

#$1 for the index and $2 for the directory
metaphlan --install --index "$1" --db_dir $2/metaphlan
