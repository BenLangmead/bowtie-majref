#!/bin/sh

set -ex

AR=h38_1kgmaj.fa.gz
DST_AR=grch38_1kgmaj.fa.gz

for PTH in ftp.ccb/pub/data/bowtie_indexes ftp.ccb/pub/data/bowtie2_indexes ; do
    scp ${AR} gwln1:${PTH}/${DST_AR}
    ssh gwln1 chmod 664 ${PTH}/${DST_AR}
done
