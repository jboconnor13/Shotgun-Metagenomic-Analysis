#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --job-name=env_creation
#SBATCH --nodes=1 # use 1 node 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00 # Time limit days-hrs:min:sec
#SBATCH --qos=normal
#SBATCH --mem=10gb # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=john.2.oconnor@cuanschutz.edu

module purge 
module load miniforge
module load python/3.10.2

#The environment is created
conda create -n assembly

#The environment is activated
conda activate assembly

#Install metabat2 (MAG assembly)
conda install bioconda::metabat2

#Install DAS tool (Bin refinement)
conda install bioconda::das_tool

#Install CheckM (Bin quality)  
conda install bioconda::checkm-genome

#Install dREP (MAG dereplication)
conda install bioconda::drep

conda install -c conda-forge -c bioconda -c defaults prokka -y
conda install -c bioconda megahit trim-galore -y
conda install -c bioconda blast bwa diamond -y
conda install -c bioconda hmmer -y
conda install -c bioconda samtools bedtools seqkit -y
conda install -c bioconda kraken2 -y
conda install -c agbiome bbtools
conda install -c bioconda seqtk flye minimap2
conda install -c conda-forge -c bioconda mmseqs2
conda install dbcan -c conda-forge -c bioconda

#Export to yaml
conda env export > ../envs/assembly.yaml
