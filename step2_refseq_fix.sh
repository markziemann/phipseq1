#!/bin/bash

# reference sequences were obtained from BTV_FMD_final_oligos.csv
cut -f-2 refseq.tsv  | sed 1d | sed 's@^@>@' | tr '\t' '\n' > refseq.fa
