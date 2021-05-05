#!/bin/bash

set -ex

# 1. Download GRCh37 genome
mkdir -p hg19
REF="hg19.fa"
if [[ ! -f ${REF} ]] ; then
    cd hg19
    wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*.fa.gz
    gunzip -f *.fa.gz
    for f in `ls -1 | sort -k1,1 -V`; do
        cat $f >> ../${REF}
    done
    cd ..
fi

VCF="ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz"
# 2. Download WGS vcf from 1000 Genomes Project
if [[ ! -f ${VCF} ]] ; then
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
fi


# 3. Extract major alleles
filter_major_snv() {
    # Get major-alleles (SNVs) from a VCF.
    # $1: path to the input VCF file
    # $2: path to the output VCF file
    IN=$1
    OUT=$2
    bcftools view -O z -m 2 -M 2 -v snps -q 0.5 -e 'AF = 0.5' -o ${OUT} ${IN}
    bcftools index ${OUT}
}

filter_major_snv_indel() {
    # Get major-alleles (SNVs and indels) from a VCF.
    # $1: path to the input VCF file
    # $2: path to the output VCF file
    IN=$1
    OUT=$2
    bcftools view -O z -m 2 -M 2 -v snps,indels -q 0.5 -e 'AF = 0.5' -o ${OUT} ${IN}
    bcftools index ${OUT}
}


MAJOR_VCF_SNV_RAW="wgs.concat.major_snvs.vcf.gz"
MAJOR_VCF_SNV="wgs.concat.major_snvs.renamed.vcf.gz"
if [[ ! -f ${MAJOR_VCF_SNV_RAW} ]]; then
    filter_major_snv ${VCF} ${MAJOR_VCF_SNV_RAW}
fi
if [[ ! -f ${MAJOR_VCF_SNV}.csi ]]; then
    bcftools annotate --rename-chrs hg19.chrom_map -O z -o ${MAJOR_VCF_SNV} ${MAJOR_VCF_SNV_RAW}
    bcftools index ${MAJOR_VCF_SNV}
fi

MAJOR_VCF_SNV_INDEL_RAW="wgs.concat.major_snvindels.vcf.gz"
MAJOR_VCF_SNV_INDEL="wgs.concat.major_snvindels.renamed.vcf.gz"
if [[ ! -f ${MAJOR_VCF_SNV_INDEL_RAW} ]]; then
    filter_major_snv_indel ${VCF} ${MAJOR_VCF_SNV_INDEL_RAW}
fi
if [[ ! -f ${MAJOR_VCF_SNV_INDEL}.csi ]]; then
    bcftools annotate --rename-chrs hg19.chrom_map -O z -o ${MAJOR_VCF_SNV_INDEL} ${MAJOR_VCF_SNV_INDEL_RAW}
    bcftools index ${MAJOR_VCF_SNV_INDEL}
fi

# if [[ ! -f hg19.major.vcf.gz ]]; then
#     bcftools annotate --rename-chrs hg19.chrom_map -O z -o  hg19.major.vcf.gz ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.major.vcf.gz
# fi
# if [[ ! -f hg19.major.vcf.gz.csi ]]; then
#     bcftools index hg19.major.vcf.gz
# fi

if [[ ! -f hg19_1kgmaj_snvs.fa ]]; then
    bcftools consensus -f ${REF} -o hg19_1kgmaj_snvs.fa ${MAJOR_VCF_SNV}
fi
if [[ ! -f hg19_1kgmaj_snvindels.fa ]]; then
    bcftools consensus -f ${REF} -o hg19_1kgmaj_snvindels.fa ${MAJOR_VCF_SNV_INDEL}
fi

# Index REF
if [[ ! -f ${REF}.fai ]]; then
    samtools faidx ${REF}
fi
# Build levioSAM index
if [[ ! -f hg19_1kgmaj_snvindels.lft ]]; then
    cut -f 1-2 ${REF}.fai > hg19.length_map
    leviosam serialize -v ${MAJOR_VCF_SNV_INDEL} -k hg19.length_map -p hg19_1kgmaj_snvindels
fi
