#!/bin/bash -l
#SBATCH
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --output=.build38.out
#SBATCH --error=.build38.err
#SBATCH --time=4:00:00
#SBATCH --ntasks-per-node=4

./build38.sh
