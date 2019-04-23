#!/bin/bash

set -ex

if [ "$#" -ne "3" ]; then
    echo "Usage: $0 prefix vcf fasta" >&3
else
    msg=$(bcftools -h | wc -l)
    if [ "$msg" -gt "1" ]; then
        echo "bcftools has already been loaded";
    else
        echo "loading bcftools...";
        module load bcftools;
    fi
    msg=$(vcftools -h | wc -l)
    if [ "$msg" -gt "1" ]; then
        echo "vcftools has already been loaded";
    else
        echo "loading vcftools...";
        module load vcftools;
    fi
    msg=$(samtools view | wc -l)
    if [ "$msg" -gt "1" ]; then
        echo "samtools has already been loaded";
    else
        echo "loading samtools...";
        module load samtools;
    fi
    PREFIX=$1
    VCF="${PREFIX}.vcf"
    VCF_R_GZ="${PREFIX}.recode.vcf.gz"
    MA_REF="${PREFIX}.fa"
    NUM_TH=16
    #: bcftools dont solve mnps
    #: vcftools need sample to calculate frequency (cannot only look at INFO:AF)
    
    bcftools view -O z -v snps --threads $NUM_TH -q 0.5 -e 'AF = 0.5' $2 | vcftools --gzvcf - --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout --remove-indels | bgzip -@ $NUM_TH > $VCF_R_GZ
    bcftools index $VCF_R_GZ
    bcftools consensus -f ${3} $VCF_R_GZ > $MA_REF

fi
