#!/bin/sh

BASE=hg19_1kgmaj
EXT=bt2

zip ${BASE}.zip \
    ${BASE}.1.${EXT} \
    ${BASE}.2.${EXT} \
    ${BASE}.3.${EXT} \
    ${BASE}.4.${EXT} \
    ${BASE}.rev.1.${EXT} \
    ${BASE}.rev.2.${EXT} \
    README.md

