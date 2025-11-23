#!/bin/bash

#SBATCH --job-name dev_mags
#SBATCH --output devMAG-%j.error
#SBATCH --error devMAG-%j.out
#SBATCH -N 1
#SBATCH --cpus-per-task 2
#SBATCH --partition=cpu
#SBATCH --mem=50G
#SBATCH --time=2:00:00

eval "$(mamba shell hook --shell bash)"
mamba activate mag
MAGBuild run --input testReads/paired --extn fq.gz --sequencing paired --host_seq /users/bnalaga1/scratch/reference_db --profile slurm


