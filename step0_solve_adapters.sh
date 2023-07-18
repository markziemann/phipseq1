#!/bin/bash

# Goal of this work is to identify the adapter sequences.
# Using BTV1-1 dataset initially,
# First step is to assemble forward and reverse reads
pear -j10  -f BTV1-1_S1_L001_R1_001.fastq.gz -r BTV1-1_S1_L001_R2_001.fastq.gz  -o  BTV1-1

# Next step is to see what part of the sequence is common, that is the adapter.
# identify the common 5 prime adapter
# modify x to get the common part of the sequence
x=49
sed -n '2~4p' BTV1-1.assembled.fastq | head -1000 | cut -c-$x | sort | uniq -c | sort -k1nr | head -50
# 5prime adapter = GGAGCTGTCGTATTCCAGTCAGGTGTGATGCTCGGGGATCCGAATTCA

# identify the 3' adapter
# modify y to get the common sequence
sed -n '2~4p' BTV1-1.assembled.fastq | head -1000 | rev | cut -c-67 | rev | sort | uniq -c | sort -k1nr | head -50
# 3prime adapter = AAGCTTGCGGCCGCACTCGAGTAACTAGTTAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTT
