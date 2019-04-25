#!/bin/bash -l
#SBATCH
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --output=.index_bt.out
#SBATCH --error=.index_bt.er
#SBATCH --time=4:00:00
#SBATCH --ntasks-per-node=4

bowtie-build h37_1kgmaj.fa hg19_1kgmaj

