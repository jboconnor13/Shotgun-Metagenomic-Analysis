# Shotgun-Metagenomics-Pipeline

This repository contains the entire shotgun metagenomic analysis pipeline first created by J. O'Connor in 2025.  

## Description of Directories/Files:
- **config**: contains config specifying directories of inputs and outputs  
- **envs**: contains the conda environment YAMLs used for each rule  
- **ref_databases**: directory with reference databases  
- **Snakefile**: contains the script for all the analysis  

## Conda Environment Setup

First create a new envirnment and activate it 
```bash
conda create -n shotgun_analysis
conda activate shotgun_analysis
conda config --set channel_priority strict
```

## Run the Snakemake pipeline
```bash
cd scripts

chmod +x Run.sh

bash Run.sh
```
