#!/bin/bash

set -ex

REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# 1. Download hg38 genome
if [[ ! -f ${REF} ]]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    bgzip -d ${REF}.gz
fi

# GRCh38 p13
# if [[ -z GCF_000001405.39_GRCh38.p13_genomic.fna ]]; then
#     https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
#     bgzip -d GCF_000001405.39_GRCh38.p13_genomic.fna.gz
# fi

# 2. Download WGS VCFs from 1000 Genomes Project
if [[ ! -f ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz ]]; then
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz
fi
# Convert contig namings from "1" to "chr1"
if [[ ! -f ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.renamed.vcf.gz ]]; then 
    bcftools annotate --rename-chrs GRCh38.chrom_map -O z -o ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.renamed.vcf.gz ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz
fi

# We use the "lift-over" calls in the non-pseudo-autosomal regions since they are not included
# in the GRCh38-based 1KGP calls we use.
# 
# See:
# "Variants on the Y chromosome or in regions of the X chromosome outside the pseudo-autosomal 
#  regions were discarded due to the ploidy settings used in this work."
#       Lowy-Gallego E, Fairley S, Zheng-Bradley X et al. Variant calling on the GRCh38 assembly 
#       with the data from phase three of the 1000 Genomes Project. 
#       Wellcome Open Res 2019, 4:50 (https://doi.org/10.12688/wellcomeopenres.15126.2)

if [[ ! -f ALL.chrY_GRCh38_sites.20170504.vcf.gz ]]; then
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chrY_GRCh38_sites.20170504.vcf.gz
fi
if [[ ! -f ALL.chrY_GRCh38_sites.20170504.renamed.vcf.gz ]]; then 
    bcftools annotate --rename-chrs GRCh38.chrom_map -O z -o ALL.chrY_GRCh38_sites.20170504.renamed.vcf.gz ALL.chrY_GRCh38_sites.20170504.vcf.gz
fi
if [[ ! -f ALL.chrX_GRCh38_sites.20170504.vcf.gz ]]; then
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chrX_GRCh38_sites.20170504.vcf.gz
fi
if [[ ! -f ALL.chrX_GRCh38_sites.20170504.renamed.vcf.gz ]]; then 
    bcftools annotate --rename-chrs GRCh38.chrom_map -O z -o ALL.chrX_GRCh38_sites.20170504.renamed.vcf.gz ALL.chrX_GRCh38_sites.20170504.vcf.gz
fi

WG_VCF="wgs.concat.vcf.gz"
if [[ ! -f ${WG_VCF} ]]; then
    bcftools concat ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.renamed.vcf.gz ALL.chrX_GRCh38_sites.20170504.renamed.vcf.gz ALL.chrY_GRCh38_sites.20170504.renamed.vcf.gz | bcftools sort -O z -o ${WG_VCF}
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

MAJOR_VCF_SNV="wgs.concat.major_snvs.vcf.gz"
# Obtain major-alleles. Output: wgs.concat.major_snvs.vcf.gz or wgs.concat.major_snvindels.vcf.gz
if [[ ! -f ${MAJOR_VCF_SNV} ]]; then
    filter_major_snv ${WG_VCF} ${MAJOR_VCF_SNV}
fi
MAJOR_VCF_SNV_INDEL="wgs.concat.major_snvindels.vcf.gz"
if [[ ! -f ${MAJOR_VCF_SNV_INDEL} ]]; then
    filter_major_snv_indel ${WG_VCF} ${MAJOR_VCF_SNV_INDEL}
fi

# 4. Build the major-alelle reference
if [[ ! -f grch38_1kgmaj_snvs.fa ]]; then
    bcftools consensus -f ${REF} -o grch38_1kgmaj_snvs.fa ${MAJOR_VCF_SNV}
fi
if [[ ! -f grch38_1kgmaj_snvindels.fa ]]; then
    bcftools consensus -f ${REF} -o grch38_1kgmaj_snvindels.fa ${MAJOR_VCF_SNV_INDEL}
fi

# Index REF
if [[ ! -f ${REF}.fai ]]; then
    samtools faidx ${REF}
fi
# Build levioSAM index
if [[ ! -f grch38_1kgmaj_snvindels.lft ]]; then
    cut -f 1-2 ${REF}.fai > GRCh38.length_map
    leviosam serialize -v ${MAJOR_VCF_SNV_INDEL} -k GRCh38.length_map -p grch38_1kgmaj_snvindels
fi
