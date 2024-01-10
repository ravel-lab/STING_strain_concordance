#!/bin/bash

#use the current working directory and current modules
#$ -cwd -V

#$ -b y -l mem_free=256G -P jravel-lab -q threaded.q -pe thread 8 -N inStrain_compare -j y -o /local/scratch/kaylacarter2/logs/ -e /local/scratch/kaylacarter2/logs/

cd /local/scratch/kaylacarter2/STING/03_inStrainCompare

source /home/mfrance/software/inStrain/.venv/bin/activate
source /home/kaylacarter2/.bashrc

ind_dir=/local/scratch/kaylacarter2/STING/01_indexDB

inStrain compare -i MG* -p 8 -s $ind_dir/magKeys/STING_magKey.txt --clusterAlg ward -o STING_inStrain_compare
