#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --output=/scratch/alpine/joconnor@xsede.org/Shotgun-Metagenomic-Analysis/slurm_outputs/slurm-%j.out
#SBATCH --job-name=kraken_db_download
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

mkdir -p ../ref_databases/kraken

cd ../ref_databases/kraken

#The standard database is downloaded from the kraken website
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20251015.tar.gz

#The downloaded file is unpacked
tar -xzvf k2_standard_20251015.tar.gz

#The old .tar.gz is removed for cleanup
rm k2_standard_20251015.tar.gz
