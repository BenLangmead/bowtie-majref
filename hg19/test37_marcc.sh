#!/bin/bash -l
#SBATCH
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --output=.test37.out
#SBATCH --error=.test37.err
#SBATCH --time=4:00:00
#SBATCH --ntasks-per-node=2

./test37.sh
