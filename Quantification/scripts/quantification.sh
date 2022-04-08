#!/bin/bash

conda activate salmon

for k in /home/s9/tsai/work/Ferenc/Datasets/*
do

cd "$k"

DIR=$(echo ${k})
GTF=$(ls ${k}/*gtf.gz)
GFF=$(ls ${k}/*gff3.gz)
GENOME=$(ls ${k}/*toplevel.fa.gz)
CDNA=$(ls ${k}/*cdna.all.fa.gz)

echo $DIR
 
grep "^>" <(gunzip -c $GENOME ) | cut -d " " -f1 > ${DIR}/decoys.txt
sed -i.bak -e 's/>//g' ${DIR}/decoys.txt
cat $CDNA $GENOME > ${DIR}/gentrome.fa.gz

   for i in ${DIR}/*_1_fastp.fastq.gz
   do
     SAMPLE1=$(echo ${i} | awk -F "/" '{print $9}' | sed "s/_1_\.fastp\.fastq\.gz//")
     SAMPLE2=$(echo ${i} | awk -F "/" '{print $9}' | sed "s/_1_fastp\.fastq\.gz/_2_fastp.fastq.gz/")
     salmon quant -i ${DIR}/salmon_index -l A -1 ${DIR}/${SAMPLE1} -2 ${DIR}/${SAMPLE2} -o ${DIR}/${SAMPLE1}_salmon --validateMappings -p 10 --seqBias --gcBias --posBias --numBootstraps 100
   done

   for i in ${DIR}/*_RNA-Seq.fastp.fastq.gz
   do
     SAMPLE=$(echo ${i} | awk -F "/" '{print $9}')
     echo $SAMPLE
     salmon quant -i ${DIR}/salmon_index -l A -r ${DIR}/${SAMPLE} -o ${DIR}/${SAMPLE}_salmon -p 2 --seqBias --posBias --validateMappings --numBootstraps 100
   done

 done


done

