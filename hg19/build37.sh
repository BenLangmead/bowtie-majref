#!/bin/bash

set -ex

# 1. Download GRCh37 genome
mkdir -p hg19
if [[ ! -f hg19.fa ]] ; then
    cd hg19
    wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*.fa.gz
    gunzip -f *.fa.gz
    for f in `ls -1 | sort -k1,1 -V`; do
        cat $f >> ../hg19.fa
    done
    cd ..
fi

# 2. Download WGS vcf from 1000 Genomes Project
if [[ ! -f ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz ]] ; then
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
fi


filter_major_snp() {
    # Get major-alleles (SNVs) from a VCF.
    # $1: the prefix of a VCF.GZ file
    OUT="${1}.major.vcf.gz"
    bcftools view -O z -m 2 -M 2 -v snps -q 0.5 -e 'AF = 0.5' -o ${OUT} ${1}.vcf.gz
    bcftools index ${OUT}
}

filter_major_snp_indel() {
    # Get major-alleles (SNVs and indels) from a VCF.
    # $1: the prefix of a VCF.GZ file
    OUT="${1}.major.vcf.gz"
    bcftools view -O z -m 2 -M 2 -v snps,indels -q 0.5 -e 'AF = 0.5' -o ${OUT} ${1}.vcf.gz
    bcftools index ${OUT}
}

if [[ ! -f ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.major.vcf.gz ]]; then
    filter_major_snp ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites
fi
if [[ ! -f hg19.major.vcf.gz ]]; then
    bcftools annotate --rename-chrs hg19.chrom_map -O z -o  hg19.major.vcf.gz ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.major.vcf.gz
fi
if [[ ! -f hg19.major.vcf.gz.csi ]]; then
    bcftools index hg19.major.vcf.gz
fi

bcftools consensus -f hg19.fa -o hg19_1kgmaj.fa hg19.major.vcf.gz
