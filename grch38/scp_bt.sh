#!/bin/sh

AR=grch38_1kgmaj_bt.zip
PTH=ftp.ccb/pub/data/bowtie_indexes

scp ${AR} gwln1:${PTH}/
ssh gwln1 chmod 664 ${PTH}/${AR}
