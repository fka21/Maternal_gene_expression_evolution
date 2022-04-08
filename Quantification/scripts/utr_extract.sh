#!/bin/bash

conda activate ortho

while read k 
do

cd "$k"

DIR=$(echo ${k})
CDNA=$(ls ${k}/*cdna.all.fa)
CDS=$(ls ${k}/*cds.all.fa)
SPECIES=$(echo ${k} | awk -F "/" '{print $8}')

perl ~/work/Ferenc/Tools/ExUTR/bin/3UTR_ext_20170816.pl \
-i1 $CDNA \
-i2 $CDS \
-d ~/work/Ferenc/Tools/ExUTR/3UTRef.Inv.fasta \
-a 3 \
-o ${DIR}/${SPECIES}_3UTR.fasta \
-x 2500 \
-m 200

done < /home/s9/tsai/work/Ferenc/Datasets/loop.txt #loop.txt contains filepaths for directories assigned for different species

