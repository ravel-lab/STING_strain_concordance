#!/bin/bash

#use the current working directory and current modules
#$ -cwd -V

#$ -b y -l mem_free=2G -P jravel-lab -N mkdir -j y -o /local/scratch/kaylacarter2/logs/ -e /local/scratch/kaylacarter2/logs/

mkdir /local/scratch/kaylacarter2/STING/reads/
mkdir /local/scratch/kaylacarter2/STING/MAGs/
mkdir /local/scratch/kaylacarter2/STING/01_indexDB/
mkdir /local/scratch/kaylacarter2/STING/01_indexDB/genomeDB/
mkdir /local/scratch/kaylacarter2/STING/01_indexDB/geneFiles/
mkdir /local/scratch/kaylacarter2/STING/01_indexDB/magKeys/
mkdir /local/scratch/kaylacarter2/STING/01_indexDB/indexB2/
mkdir /local/scratch/kaylacarter2/STING/02_mapping/
mkdir /local/scratch/kaylacarter2/STING/02_mapping/bamFiles/
mkdir /local/scratch/kaylacarter2/STING/03_inStrainCompare/
