#!/bin/bash -l
#SBATCH
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --output=.index_bt.out
#SBATCH --error=.index_bt.err
#SBATCH --time=5:00:00
#SBATCH --ntasks-per-node=4

bowtie-build hg19_1kgmaj.fa hg19_1kgmaj

