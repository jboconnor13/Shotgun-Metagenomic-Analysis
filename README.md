# Shotgun-Metagenomics-Pipeline

This repository contains the entire shotgun metagenomic analysis pipeline first created by J. O'Connor in 2025.  

## Description of Directories/Files:
- **config**: contains config specifying directories of inputs and outputs  
- **envs**: contains the conda environment YAMLs used for each rule  
- **ref_databases**: directory with reference databases  
- **Snakefile**: contains the script for all the analysis  

## Conda Environment Setup

```bash

conda create -n shotgun_analysis
conda activate shotgun_analysis
conda config --set channel_priority strict

```
## Reference Database Downloads

Reference databases need to downloaded for host filtering, taxonomic classification, and functional classification

## Metaphlan

```bash

#Navigate to the directory
cd scripts/Database\ Downloads/

#Make the Metaphlan Database download script executable
chmod +x Alpine_Metaphlan_Database.Download.sh

#Run the metaphlan download script specifying the index/version and directory for the database
sbatch Alpine_Metaphlan_Database_Download.sh \
    mpa_vJan21_CHOCOPhlAnSGB_202103 \
    /pl/active/LozuponeLab/JOC_ref_databases
```

## Hostile

```bash

#Make the Hostile Database download script executable
chmod +x Alpine_Hostile_Database_Download.sh

#Run the Hostile download script specifying the directory for the database
sbatch Alpine_Hostile_Database_Download.sh \
    /pl/active/LozuponeLab/JOC_ref_databases
```

## HuMAnn

```bash

#Make the Hostile Database download script executable
chmod +x Alpine_Humann_Database_Download.sh

#Run the Hostile download script specifying the directory for the database
sbatch Alpine_Hostile_Database_Download.sh \
    /pl/active/LozuponeLab/JOC_ref_databases
```

## Kraken

```bash

#Make the Hostile Database download script executable
chmod +x Alpine_Kraken_Database_Download.sh

#Run the Hostile download script specifying the directory for the database
sbatch Alpine_Kraken_Database_Download.sh \
    /pl/active/LozuponeLab/JOC_ref_databases
``


