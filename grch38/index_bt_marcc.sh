#!/bin/bash -l
#SBATCH
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --output=.index_bt.out
#SBATCH --error=.index_bt.er
#SBATCH --time=4:00:00
#SBATCH --ntasks-per-node=4

bowtie-build h38_1kgmaj.fa grch38_1kgmaj

