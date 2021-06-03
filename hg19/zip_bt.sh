#!/bin/sh

BASE="hg19_1kgmaj_snvs"
EXT=ebwt

zip ${BASE}_bt.zip \
    ${BASE}.1.${EXT} \
    ${BASE}.2.${EXT} \
    ${BASE}.3.${EXT} \
    ${BASE}.4.${EXT} \
    ${BASE}.rev.1.${EXT} \
    ${BASE}.rev.2.${EXT} \
    wgs.concat.major_snvs.vcf.gz \
    wgs.concat.major_snvs.vcf.gz.csi \
    README.md

BASE="hg19_1kgmaj_snvindels"
EXT=ebwt

zip ${BASE}_bt.zip \
    ${BASE}.1.${EXT} \
    ${BASE}.2.${EXT} \
    ${BASE}.3.${EXT} \
    ${BASE}.4.${EXT} \
    ${BASE}.rev.1.${EXT} \
    ${BASE}.rev.2.${EXT} \
    ${BASE}.lft \
    wgs.concat.major_snvindels.vcf.gz \
    wgs.concat.major_snvindels.vcf.gz.csi \
    README.md

