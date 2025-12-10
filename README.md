# Shotgun-Metagenomics-Pipeline

This repository contains the entire shotgun metagenomic analysis pipeline first created by J. O'Connor in 2025.  

## Description of Directories/Files:
- **config**: contains config specifying directories of inputs and outputs  
- **envs**: contains the conda environment YAMLs used for each rule  
- **ref_databases**: directory with reference databases
- **scripts**: Includes scripts for running the pipeline locally or on alpine
- **snakefile**: contains the script for all the analysis  

The steps below indicate how to run the pipeline on ALpine

## Conda Environment Setup
First create a new envirnment and activate it 
```bash
cd scripts

chmod +x Alpine_Environment_Creation.sh

sbatch Alpine_Environment_Creation.sh
```

##Database download
The hostile, kraken, and humann databases need to be downloaded

```bash
chmod +x Alpine_Kraken_Database_Download.sh
sbatch Alpine_Kraken_Database_Download.sh

chmod +x Alpine_Hostile_Database_Download.sh
sbatch Alpine_Hostile_Database_Download.sh

chmod +x Alpine_Humann_Database.Download.sh
sbatch Alpine_Humann_Database.Download.sh
```

## Config, metadata, and adapter file adjustment 

The following need to be adjusted
- **config/config.yaml**: specify inputs and parameters specific to project  
- **example_data/example_metadata**: adjust to include your samples and the locations of forward and reverse reads  
- **scripts**:  for alpine bash files specify your username and email for sbatch submissions  

## Run the Snakemake pipeline
Then run the pipeline

```bash
chmod +x Alpine_Run.sh

sbatch Alpine_Run.sh
```
