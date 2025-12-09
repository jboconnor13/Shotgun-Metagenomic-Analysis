#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --output=/scratch/alpine/joconnor@xsede.org/Shotgun-Metagenomic-Analysis/slurm_outputs/slurm-%j.out
#SBATCH --job-name=humann_db_download
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

#We first install humann
pip install humann


#Now we create a directory for humann database
mkdir -p ../ref_databases/humann

#We create a subdirectories for Uniref and Chocophlan
mkdir -p ../ref_databases/humann/chocophlan
mkdir -p ../ref_databases/humann/uniref

#We download Chocophlan
humann_databases --download chocophlan full ../ref_databases/humann/chocophlan

#We download Uniref
humann_databases --download uniref uniref90_diamond ../ref_databases/humann/uniref

#We also create a subdirectory for the utility mapping
mkdir -p ../ref_databases/humann/utility_mapping

#We download the utlity mapping
humann_databases --download utility_mapping full ../ref_databases/humann/utility_mapping