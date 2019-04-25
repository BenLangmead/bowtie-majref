#!/bin/sh

AR=hg19_1kgmaj_bt2.zip
PTH=ftp.ccb/pub/data/bowtie2_indexes

scp ${AR} gwln1:${PTH}/
ssh gwln1 chmod 664 ${PTH}/${AR}
