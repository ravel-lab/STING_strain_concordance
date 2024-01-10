#!/bin/bash

#use the current working directory and current modules
#$ -cwd -V

#$ -b y -l mem_free=48G -P jravel-lab -q threaded.q -pe thread 4 -N indexDB -j y -o /local/scratch/kaylacarter2/logs/ -e /local/scratch/kaylacarter2/logs/

DB_dir=/local/scratch/kaylacarter2/STING/01_indexDB/genomeDB

cd /local/scratch/kaylacarter2/STING/01_indexDB/genomeDB

source /home/mfrance/software/inStrain/.venv/bin/activate
source /home/kaylacarter2/.bashrc

cat /local/scratch/kaylacarter2/STING/MAGs/*.fasta > $DB_dir/STING_mags.fa

prodigal -i $DB_dir/STING_mags.fa -d ../geneFiles/STING_mags.fna

parse_stb.py --reverse -f /local/scratch/kaylacarter2/STING/MAGs/*.fasta -o ../magKeys/STING_magKey.txt

bowtie2-build -f $DB_dir/STING_mags.fa ../indexB2/STING_mags
