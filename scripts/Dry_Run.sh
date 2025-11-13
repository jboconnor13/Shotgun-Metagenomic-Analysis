#!/bin/bash

#First we go to the main directory to run the snakefile
cd ..

#This script runs the entire shotgun pipeline

snakemake -s Snakefile \
    --configfile config/config.yaml \
    -n \
    --use-conda

#-s specifes the snakefile to use
#--configfile secifies the config file to use
#-n is for the dry-run
#--use-conda specifies to sue conda environements put in each rule and to create them if they are not already present
