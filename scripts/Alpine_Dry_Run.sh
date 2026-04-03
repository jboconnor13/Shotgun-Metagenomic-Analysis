#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --output=/scratch/alpine/joconnor@xsede.org/Shotgun-Metagenomic-Analysis/slurm_outputs/slurm-%j.out
#SBATCH --job-name=SG
#SBATCH --nodes=1 # use 1 node 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=16
#SBATCH --time=00:15:00 # Time limit days-hrs:min:sec
#SBATCH --qos=normal
#SBATCH --mem=150G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=john.2.oconnor@cuanschutz.edu

module purge 
module load miniforge
module load python/3.10.2

conda activate shotgun_analysis

#First we go to the main directory to run the snakefile
cd ..

#This script runs the entire shotgun pipeline

snakemake -s snakefile \
    --latency-wait 3000 \
    --configfile config/config.yaml \
    --cores 8 \
    --use-conda \
    --conda-prefix "/projects/.xsede.org/joconnor/.snakemake/conda/" \
    --dry-run
