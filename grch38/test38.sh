set -ex

MAJ="grch38_1kgmaj_snvs.fa"
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

if [[ -f wgs.concat.vcf.gz ]]; then
    VCF="wgs.concat.vcf.gz"
    bgzip -cd ${VCF} | \
        python ../scripts/vcf_processing.py --out_var_loc --min_af 0.5 --max_allele_len 1 | \
        python ../scripts/test_major_allele_ref.py -m ${MAJ} \
        -r ${REF} > test_h38.txt
else
    VCF="ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.renamed.vcf.gz"
    VCF_X="ALL.chrX_GRCh38_sites.20170504.vcf.gz"
    VCF_Y="ALL.chrY_GRCh38_sites.20170504.vcf.gz"
    bgzip -cd ${VCF} | \
        python ../scripts/vcf_processing.py --out_var_loc --min_af 0.5 --max_allele_len 1 | \
        python ../scripts/test_major_allele_ref.py -m ${MAJ} -r ${REF} > test_h38.txt

    bgzip -cd ${VCF_Y} | \
        python ../scripts/vcf_processing.py --out_var_loc --min_af 0.5 --max_allele_len 1 | \
        python ../scripts/test_major_allele_ref.py -m ${MAJ} -r ${REF} -t chrY > test_h38_y.txt

    bgzip -cd ${VCF} | \
        grep 'chrX' | \
        python ../scripts/vcf_processing.py --out_var_loc --min_af 0.5 --max_allele_len 1 > var_x_h38.txt
    bgzip -cd ${VCF_X} | \
        python ../scripts/vcf_processing.py --out_var_loc --min_af 0.5 --max_allele_len 1 >> var_x_h38.txt
    python ../scripts/test_major_allele_ref.py -m ${MAJ} \
        -r ${REF} -t chrX -c var_x_h38.txt > test_h38_x.txt
fi
