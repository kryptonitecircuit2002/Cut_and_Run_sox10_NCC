#!/bin/bash

#SBATCH --job-name=cutandrun_ncc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=cutandrun_ncc_%j.out
#SBATCH --error=cutandrun_ncc_%j.err
#SBATCH --partition=compute

#system modules
module load anaconda3

cd ../Cut_and_Run/

#Initialise conda
eval "$(conda shell.bash hook)"

#activate conda environment
conda activate Nextflow


nextflow run nf-core/cutandrun -profile singularity \
  --input ../Samplesheet_igg1.csv \
  --peakcaller 'seacr,MACS2' \
  --genome .../.gtf \
  --outdir .../outputfiles/
