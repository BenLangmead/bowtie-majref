#!/bin/bash -l
#SBATCH
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --output=.index_bt2.out
#SBATCH --error=.index_bt2.err
#SBATCH --time=5:00:00
#SBATCH --ntasks-per-node=4

bowtie2-build h38_1kgmaj.fa grch38_1kgmaj

