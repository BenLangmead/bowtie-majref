#!/bin/sh

set -ex

AR=h37_1kgmaj.fa.gz
DST_AR=hg19_1kgmaj.fa.gz

for PTH in ftp.ccb/pub/data/bowtie_indexes ftp.ccb/pub/data/bowtie2_indexes ; do
    scp ${AR} gwln1:${PTH}/${DST_AR}
    ssh gwln1 chmod 664 ${PTH}/${DST_AR}
done
