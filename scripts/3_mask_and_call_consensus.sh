#!/bin/bash

#tools: bedops, samtools, bedtools, bcftools (version 1.7)
samtools="samtools"
bcftools="scripts/bcftools-1.7/bcftools"
bedtools="bedtools"
bedops="bedops"


ref="ref_245/c_elegans.PRJNA13758.WS245.genomic.fa"
region="I:11953012-11962484"


cd BAMs/gene;

for file in $(ls *gene.bam);do

  #low-coverage sites
  samtools depth -a -r $region ${file/.bam/.re.bam}| awk '{ if ($3 < 3) { print $1 "\t" $2-1 "\t" $2 }}' - > ${file/.bam/.DP3.mask.bed}

  #extract high-quality indels
  zcat ${file/.bam/filt.soft}.vcf.gz |vcf2bed --deletions - > ${file/.bam/.indel.bed}
  zcat ${file/.bam/filt.soft}.vcf.gz |vcf2bed --insertions - >> ${file/.bam/.indel.bed}

  $bedtools sort -i ${file/.bam/.indel.bed} > ${file/.bam/.indel.s.bed}
  #remove these indels from the mask
  $bedtools merge -i ${file/.bam/.DP3.mask.bed} | $bedtools subtract -a - -b ${file/.bam/.indel.s.bed}  > ${file/.bam/.DP3.mask.noindel.bed}


#call consensus + IUPAC codes (works only with bcftools version 1.7)
name=${file/.gene.bam/}
$samtools faidx $ref $region | $BCFtools consensus -I -s $name -m ${file/.bam/.DP3.mask.noindel.bed} ${file/.bam/filt.soft}.vcf.gz > ${file/.bam/.0.fasta}

#edit names
sed -i "s/I/$name-gene-I/" ${file/.bam/.0.fasta}

