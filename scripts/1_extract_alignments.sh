#!/bin/bash

#Extract region(s) from the alignments
#edit the path if the tool is not in the PATH
samtools="samtools"

###example I:11953012-11962484; if more then 1 region, use samtools with "-L regions.bed"
region="I:11953012-11962484"


cd BAMs;
#create a new directory if it not exists
mkdir -p gene;


for file in $(ls *.bam);do

  $samtools view -b $file $region > gene/${file/bam/gene.bam}
  $samtools index gene/${file/bam/gene.bam}

done
