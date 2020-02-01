#!/bin/bash

set -ex

#   # 1. Download hg38 genome
#   wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#   bgzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#   
#   # 2. Download WGS vcf from 1000 Genomes Project
#   wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz
#   
#   # 3. Convert the header from v4.3 to v4.2
#   sh ../scripts/convert_header_v43_to_v42.sh ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites
#   
# 4. Build major allele ref (no chrY)
sh ../scripts/build_major_allele_ref.sh h38_1kgmaj_no_y ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites_v42.vcf.gz GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# 5. Prepare chrY genome
mkdir h38_liftover
cd h38_liftover
samtools faidx ../GCA_000001405.15_GRCh38_no_alt_analysis_set.fna chrY > chrY.fa
# rename chrY to Y to fit the naming of the liftover chrY vcf
sed -i 's/chr//' chrY.fa
mv ../h38_1kgmaj_no_y.fa .
# rename chrX to X to fit the naming of the liftover chrX vcf
sed -i 's/chrX/X/' h38_1kgmaj_no_y.fa

# 6. Download chrY vcf from 1000 Genomes Project (GRCh38-liftover)
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chrY_GRCh38_sites.20170504.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chrX_GRCh38_sites.20170504.vcf.gz

# 7. build major allele ref
sh ../../scripts/build_major_allele_ref.sh h38_1kgmaj_y ALL.chrY_GRCh38_sites.20170504.vcf.gz chrY.fa
sh ../../scripts/build_major_allele_ref.sh h38_1kgmaj_no_y_update_x ALL.chrX_GRCh38_sites.20170504.vcf.gz h38_1kgmaj_no_y.fa

# rename X and Y to chrX and chrY to fit GRCh38 naming rules
sed -i 's/>Y/>chrY/' chrY.fa
sed -i 's/>X/>chrX/' h38_1kgmaj_no_y_update_x.fa

# 8. Merge two genomes
cd ..
mv h38_liftover/h38_1kgmaj_no_y_update_x.fa .
python ../scripts/replace_fa_seq.py -f h38_1kgmaj_no_y_update_x.fa -s chrY -n h38_liftover/h38_1kgmaj_y.fa > h38_1kgmaj.fa
# python ../scripts/replace_fa_seq.py -f h38_1kgmaj_no_y_update_x.fa -s Y -n h38_liftover/h38_1kgmaj_y.fa > h38_1kgmaj.fa
