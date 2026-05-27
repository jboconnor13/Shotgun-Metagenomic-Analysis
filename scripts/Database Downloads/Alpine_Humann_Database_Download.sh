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
#pip install humann


#Now we create a directory for humann database
mkdir -p $1/humann

#We download Chocophlan
humann_databases --download chocophlan full $1/humann

#We download Uniref
humann_databases --download uniref uniref90_diamond $1/humann

#We download the utlity mapping
humann_databases --download utility_mapping full $1/humann

cd $1/humann

for f in *.bt2l; do 
    ln -s "$f" "${f%.bt2l}.bt2"
done


# Add this before the humann call in your script
mkdir -p /scratch/alpine/joconnor@xsede.org/metaphlan_clean

# Link only the bowtie indices and the pkl file (ignore the .tar, .md5, and .bz2)
ln -sf /pl/active/LozuponeLab/JOC_ref_databases/humann/metaphlan/*.bt2* /scratch/alpine/joconnor@xsede.org/metaphlan_clean/
ln -sf /pl/active/LozuponeLab/JOC_ref_databases/humann/metaphlan/*.pkl /scratch/alpine/joconnor@xsede.org/metaphlan_clean/