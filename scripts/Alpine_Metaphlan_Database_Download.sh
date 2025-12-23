#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --output=/scratch/alpine/joconnor@xsede.org/Shotgun-Metagenomic-Analysis/slurm_outputs/slurm-%j.out
#SBATCH --job-name=metaphlan_db_download
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --qos=normal
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=john.2.oconnor@cuanschutz.edu

module purge 
module load miniforge
module load python/3.10.2

#Load the conda humann environment
conda activate /projects/.xsede.org/joconnor/.snakemake/conda/881d1b15caa13c225b3068b60e42350f_

metaphlan --install --bowtie2db ../ref_databases/metaphlan
