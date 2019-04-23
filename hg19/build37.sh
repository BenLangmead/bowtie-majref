#!/bin/bash

set -ex

# 1. Download GRCh37 genome
mkdir -p hg19
if [[ ! -f hg19.fa ]] ; then
    cd hg19
    wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*.fa.gz
    gunzip *.fa.gz
    cat *.fa >> ../hg19.fa
    cd ..
fi
sed -iE 's/^>chr/>/' hg19.fa

# 2. Download WGS vcf from 1000 Genomes Project
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz

# 3. build major allele ref
sh ../scripts/build_major_allele_ref.sh h37_1kgmaj ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz hg19.fa
