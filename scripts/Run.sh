#!/bin/bash

#First we go to the main directory to run the snakefile
cd ..

#This script runs the entire shotgun pipeline

snakemake -s Snakefile \
    --configfile config/config.yaml \
    --cores 8 \
    --use-conda \
    --conda-prefix /Users/johnoconnor/snakemake_envs

#-s specifes the snakefile to use
#--configfile secifies the config file to use
#--use-conda specifies to sue conda environements put in each rule and to create them if they are not already present
#--specify the conda environment file where you want environments saved
