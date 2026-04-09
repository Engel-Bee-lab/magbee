#!/bin/bash

#SBATCH --job-name dev_mags
#SBATCH --output devMAG-%j.error
#SBATCH --error devMAG-%j.out
#SBATCH -N 1
#SBATCH --cpus-per-task 16
#SBATCH --partition=cpu
#SBATCH --mem=150G
#SBATCH --time=8:00:00

#eval "$(mamba shell hook --shell bash)"
#mamba activate mag

magbee run --input //users/bnalaga1/scratch/AF_bee_reads --extn fq.gz --sequencing paired  \
    --host_seq /users/bnalaga1/scratch/reference_db/Apis_mellifera_genome/GCA_003254395.2_Amel_HAv3.1_genomic.fna \
    --conda-frontend mamba --output /users/bnalaga1/scratch/african_bee_project \
    --profile slurm -k host 

#magbee run --input testReads/paired --extn fq.gz --sequencing paired \
#    --host_seq /users/bnalaga1/scratch/reference_db --conda-frontend mamba -k \
#    --output dev_test_run binning

#this is to fix the script from 6_bin_quality onwards
#magbee run --input testReads/paired --extn fq.gz --sequencing paired \
#    --host_seq /users/bnalaga1/scratch/reference_db/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
#    --conda-frontend mamba --output output -k binning

