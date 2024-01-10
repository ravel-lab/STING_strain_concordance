#!/bin/bash

#use the current working directory and current modules
#$ -cwd -V

#$ -b y -l mem_free=250G -P jravel-lab -q threaded.q -pe thread 8 -N mapping -j y -o /local/scratch/kaylacarter2/logs/ -e /local/scratch/kaylacarter2/logs/

#setting the number of jobs to be executed
#$ -t 1-145

cd /local/scratch/kaylacarter2/STING/02_mapping/bamFiles

ind_dir=/local/scratch/kaylacarter2/STING/01_indexDB
map_dir=/local/scratch/kaylacarter2/STING/02_mapping
profile_dir=/local/scratch/kaylacarter2/STING/03_inStrainCompare
seq_dir=/local/scratch/kaylacarter2/STING/reads

infile=`sed -n -e "$SGE_TASK_ID p" /local/scratch/kaylacarter2/STING/manifest.lst`

bowtie2 --no-unal -p 4 -N 0 -X 1000 -x $ind_dir/indexB2/STING_mags -1 $seq_dir/${infile}.R1.fq.gz -2 $seq_dir/${infile}.R2.fq.gz -U $seq_dir/${infile}.unpaired.fq.gz -S ${infile}.sam

samtools view -bT $ind_dir/genomeDB/STING_mags.fa ${infile}.sam > ${infile}.bam

rm ${infile}.sam

samtools sort -o ${infile}_s.bam -T ${infile}-TEMP ${infile}.bam

rm ${infile}.bam

mv ${infile}_s.bam ${infile}.bam

source /home/mfrance/software/inStrain/.venv/bin/activate
source /home/kaylacarter2/.bashrc

cd $infile
inStrain profile --pairing_filter non_discordant --use_full_fasta_header -p 4 -g $ind_dir/geneFiles/STING_mags.fna -s $ind_dir/magKeys/STING_magKey.txt -o ${profile_dir}/${infile} $map_dir/bamFiles/${infile}.bam $ind_dir/genomeDB/STING_mags.fa
