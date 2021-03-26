#!/bin/bash -l
#SBATCH
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --output=.index_bt.out
#SBATCH --error=.index_bt.err
#SBATCH --time=5:00:00
#SBATCH --ntasks-per-node=4

bowtie-build grch38_1kgmaj.fa grch38_1kgmaj

