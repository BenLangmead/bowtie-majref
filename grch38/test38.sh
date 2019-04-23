bgzip -cd ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz | \
    python ../scripts/vcf_processing.py --out_var_loc 1 --min_af 0.5 --max_allele_len 1 | \
    python ../scripts/test_major_allele_ref.py -m h38_1kgmaj.fa \
    -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > test_h38.txt

bgzip -cd h38_liftover/ALL.chrY_GRCh38_sites.20170504.vcf.gz | \
    python ../scripts/vcf_processing.py --out_var_loc 1 --min_af 0.5 --max_allele_len 1 | \
    python ../scripts/test_major_allele_ref.py -m h38_1kgmaj.fa \
    -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -t chrY > test_h38_y.txt

bgzip -cd ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites_v42.vcf.gz | \
    grep 'chrX' | \
    python ../scripts/vcf_processing.py --out_var_loc 1 --min_af 0.5 --max_allele_len 1 > var_x_h38.txt
bgzip -cd h38_liftover/ALL.chrX_GRCh38_sites.20170504.vcf.gz | \
    python ../scripts/vcf_processing.py --out_var_loc 1 --min_af 0.5 --max_allele_len 1 >> var_x_h38.txt
python ../scripts/test_major_allele_ref.py -m h38_1kgmaj.fa \
    -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -t chrX -c var_x_h38.txt > test_h38_x.txt
