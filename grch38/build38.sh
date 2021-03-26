#!/bin/bash

set -ex

# 1. Download hg38 genome
if [[ ! -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ]]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    bgzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
fi

# GRCh38 p13
# if [[ -z GCF_000001405.39_GRCh38.p13_genomic.fna ]]; then
#     https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
#     bgzip -d GCF_000001405.39_GRCh38.p13_genomic.fna.gz
# fi

# 2. Download WGS vcf from 1000 Genomes Project
if [[ ! -f ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz ]]; then
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz
fi
# Convert contig namings from "1" to "chr1"
if [[ ! -f ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.renamed.vcf.gz ]]; then 
    bcftools annotate --rename-chrs GRCh38.chrom_map -O z -o ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.renamed.vcf.gz ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz
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

# We used the "lift-over" calls in the non-pseudo-autosomal regions since they were not included
# in the GRCh38-based 1KGP calls we used.
# 
# See:
#   Lowy-Gallego E, Fairley S, Zheng-Bradley X et al. Variant calling on the GRCh38 assembly 
#   with the data from phase three of the 1000 Genomes Project [version 2; peer review: 2 approved]. 
#   Wellcome Open Res 2019, 4:50 (https://doi.org/10.12688/wellcomeopenres.15126.2)
# "Variants on the Y chromosome or in regions of the X chromosome outside the pseudo-autosomal 
# regions were discarded due to the ploidy settings used in this work."

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

if [[ ! -f wgs.concat.vcf.gz ]]; then
    bcftools concat -O z -o wgs.concat.unsorted.vcf.gz ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.renamed.vcf.gz ALL.chrX_GRCh38_sites.20170504.renamed.vcf.gz ALL.chrY_GRCh38_sites.20170504.renamed.vcf.gz
    bcftools sort -O z -o wgs.concat.vcf.gz wgs.concat.unsorted.vcf.gz
fi

# Obtain major-alleles. Output: wgs.concat.major.vcf.gz
filter_major_snp wgs.concat

# Build the major-alelle reference
bcftools consensus -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o grch38_1kgmaj.fa wgs.concat.major.vcf.gz

