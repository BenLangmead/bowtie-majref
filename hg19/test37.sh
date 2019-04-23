bgzip -cd ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz |\
    python ../scripts/vcf_processing.py --out_var_loc 1 --min_af 0.5 --max_allele_len 1 | \
    python ../scripts/test_major_allele_ref.py -m h37_1kgmaj.fa \
    -r hg19.fa > test_h37.txt
