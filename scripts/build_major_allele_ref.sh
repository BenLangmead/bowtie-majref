#!/bin/bash

set -ex

if [ "$#" -ne "3" ]; then
    echo "Usage: $0 prefix vcf fasta" >&3
    exit 1
fi
msg=$(bcftools -h | wc -l)
if [ "$msg" -gt "1" ]; then
    echo "bcftools has already been loaded";
else
    echo "bcftools is not available. Exit";
    exit 1
fi
msg=$(samtools view | wc -l)
if [ "$msg" -gt "1" ]; then
    echo "samtools has already been loaded";
else
    echo "samtools is not available. Exit";
    exit 1
fi
PREFIX=$1
VCF="${PREFIX}.vcf"
VCF_R_GZ="${PREFIX}.recode.vcf.gz"
MA_REF="${PREFIX}.fa"

bcftools view -O z -m 2 -M 2 -v snps -q 0.5 -e 'AF = 0.5' $2 > $VCF_R_GZ
bcftools index $VCF_R_GZ
bcftools consensus -f ${3} $VCF_R_GZ > $MA_REF

