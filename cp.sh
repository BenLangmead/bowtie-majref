#!/bin/bash

set -ex

BASE=/scratch/groups/blangme2/naechyun/major_allele

mkdir -p scripts
cp ${BASE}/scripts/build_major_allele_ref.sh scripts/
#cp ${BASE}/scripts/convert_header_v43_to_v42.sh scripts/
#cp ${BASE}/scripts/replace_fa_seq.py scripts/
cp ${BASE}/scripts/test_major_allele_ref.py scripts/
cp ${BASE}/scripts/test_major_allele_ref.sh scripts/
cp ${BASE}/scripts/vcf_processing.py scripts/

mkdir -p hg19
cp ${BASE}/grch37/build37.sh hg19/
cp ${BASE}/grch37/test37.sh hg19/

mkdir -p grch38
cp ${BASE}/grch38/build38.sh grch38/
cp ${BASE}/grch38/test38.sh grch38/

