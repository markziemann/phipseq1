#!/bin/bash

# Previously I used the solve_adapters.sh script to obtain the adapter sequences
# which need to be trimmed.


FILES=$(find . | grep _R1_ | grep gz$) ; echo $FILES

for FQZ1 in $FILES ; do

# step 1: merge fwd and reverse reads
FQZ2=$(echo $FQZ1 | sed 's/_R1_/_R2_/')
BASE=$(echo $FQZ1 | cut -d '_' -f1)
echo ; echo Running $FQZ1 and $FQZ2 ; echo
echo ; echo Merging read pairs ; echo
pear -j10  -f $FQZ1 -r $FQZ2  -o  $BASE

# step 2: remove the 3' adapter
FQA=$(echo $BASE.assembled.fastq)
MERGELEN=$(sed -n '2~4p' $FQA | awk '{print length}' | numaverage | numround -n 0.1)
echo ; echo Removing 3p adater ; echo
skewer -t 10 -x adapters.fa $FQA -l 10 \
  && rm $BASE.unassembled.forward.fastq $BASE.unassembled.reverse.fastq $BASE.assembled.fastq

# step 3: in preparation to remove the 5' adapter, must reverse complement the sequence
FQAT=$BASE.assembled-trimmed.fastq
FQATRC=$BASE.assembled-trimmed_rc.fastq
echo ; echo Reverse complementing the seq ; echo
seqtk seq -r $FQAT > $FQATRC \
  && rm $FQAT

# step 4: remove the 5' adapter in  reverse complement orientation
echo ; echo Removing 5p adapter ; echo
skewer -t 10 -x TGAATTCGGATCCCCGAGCATCACACCTGACTGGAATACGACAGCTCC $FQATRC

# confirm trimming observe the difference
FQATRCT=$BASE.assembled-trimmed_rc-trimmed.fastq
PRETRIM_LEN=$(sed -n '2~4p' $FQATRC | head -100 | awk '{print length}' | numaverage)
POSTTRIM_LEN=$(sed -n '2~4p' $FQATRCT | head -100 | awk '{print length}' | numaverage)
echo ; echo before trimming: $PRETRIM_LEN bp, after trimming: $POSTTRIM_LEN bp. ; echo

# step 5: flip the orientation back
echo ; echo Reverse complementing the seq back to normal orientation ; echo
FQCLEAN=$BASE.clean.fastq
seqtk seq -r $FQATRCT > $FQCLEAN \
  && rm $FQATRC $FQATRCT

# step 6: check the resulting sequence
STARTREADS=$(zcat $FQZ1 | sed -n '2~4p' | wc -l)
ENDREADS=$(sed -n '2~4p' $FQCLEAN | wc -l)
echo ; echo Cleaning complete. No. reads at start: $STARTREADS. No. reads at end: $ENDREADS ; echo
STARTLEN=$(zcat $FQZ1 | sed -n '2~4p' | awk '{print length}' | numaverage | numround -n 0.1)
AVELEN=$(sed -n '2~4p' $FQCLEAN	 | awk '{print length}' | numaverage | numround -n 0.1)
echo ; echo Cleaning complete. Read length at start: $STARTLEN x2. Merged read length: $MERGELEN. Final read length: $AVELEN ; echo

# Step 7: make a histogram if inser lengths
sed -n '2~4p' $FQCLEAN | awk '{print length}' > len.txt
Rscript -e 'x <- as.numeric(readLines("len.txt")) ; png("hist.png") ; hist(x,xlab="frag len (bp)") ; dev.off()'
mv len.txt $BASE.len.txt
mv hist.png $BASE.hist.png

done
