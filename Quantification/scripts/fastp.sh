#!/bin/bash

conda activate quant

for k in /home/s9/tsai/work/Ferenc/Datasets/*
do

cd "$k"

DIR=$(echo ${k})

 for i in ${DIR}/*_1.fastq.gz
 do
   SAMPLE=$(echo ${i} | awk -F "/" '{print $9}' | sed "s/_1\.fastq\.gz//")
   
   echo ${SAMPLE}   

   fastp \
   -i ${DIR}/${SAMPLE}_1.fastq.gz \
   -I ${DIR}/${SAMPLE}_2.fastq.gz \
   -o ${DIR}/${SAMPLE}_1_fastp.fastq.gz \
   -O ${DIR}/${SAMPLE}_2_fastp.fastq.gz \
   -q 25 \
   -5 \
   --cut_front_window_size 5 \
   --cut_front_mean_quality 25 \
   -r \
   --cut_right_window_size 5 \
   --cut_right_mean_quality 25 \
   --correction \
   --low_complexity_filter \
   --detect_adapter_for_pe \
   -w 7 \
   -l 35
 done

  for i in ${DIR}/*.fastq
  do
   SAMPLE=$(echo ${i} | awk -F "/" '{print $9}' | sed "s/\.fastq//")

   echo ${SAMPLE}

   fastp \
   -i ${DIR}/${SAMPLE}.fastq \
   -o ${DIR}/${SAMPLE}_RNA-Seq.fastp.fastq \
   -q 25 \
   -5 \
   --cut_front_window_size 5 \
   --cut_front_mean_quality 25 \
   -r \
   --cut_right_window_size 5 \
   --cut_right_mean_quality 25 \
   --low_complexity_filter \
   -w 16 \
   -l 35
 done

done
