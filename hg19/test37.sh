if [[ ! -f ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.renamed.vcf.gz ]]; then
    bcftools annotate --rename-chrs hg19.chrom_map -O z -o  ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.renamed.vcf.gz ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
fi
bgzip -cd ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.renamed.vcf.gz |\
    python ../scripts/vcf_processing.py --out_var_loc --min_af 0.5 --max_allele_len 1 | \
    python ../scripts/test_major_allele_ref.py -m hg19_1kgmaj_snvs.fa \
    -r hg19.fa > test_h37.txt
