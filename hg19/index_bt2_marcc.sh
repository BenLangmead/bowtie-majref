#!/bin/bash -l
#SBATCH
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --output=.index_bt2.out
#SBATCH --error=.index_bt2.err
#SBATCH --time=5:00:00
#SBATCH --ntasks-per-node=4

module load bowtie2
bowtie2-build hg19_1kgmaj_snvs.fa hg19_1kgmaj_snvs
bowtie2-build hg19_1kgmaj_snvindels.fa hg19_1kgmaj_snvindels

