#!/bin/bash

#tools: samtools, bcftools, GATK, HTSlib
samtools="samtools"
bcftools="bcftools"
GATK="$EBROOTGATK/GenomeAnalysisTK.jar"
bgzip="HTSlib/1.5/bin/bgzip"
BCFtools="scripts/bcftools-1.7/bcftools"


#reference genome
ref="ref_245/c_elegans.PRJNA13758.WS245.genomic.fa"

############################### NOTE ################################
#!reference should be indexed by samtools and picard tools
#samtools faidx $ref
#picard="picard.jar"
#java -jar $picard CreateSequenceDictionary R=$ref O=${ref/.fa/.dict}
#####################################################################

cd BAMs/gene;

for file in $(ls *gene.bam);do
  echo "Processing $file";

  #realign reads with indels

  java -Xmx5g -jar $GATK -T RealignerTargetCreator -R $ref -I $file -o ${file/.bam/.intervals} ;
  java -Xmx5g -jar $GATK -T IndelRealigner -I $file -R $ref -targetIntervals ${file/.bam/.intervals} -o ${file/.bam/.re.bam};

  #call variants
  java -Xmx1g -jar $GATK -R $ref -T HaplotypeCaller -I ${file/.bam/.re.bam} -o ${file/.bam/.raw.g.vcf} -ploidy 2 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000;
  java -Xmx1g -jar $GATK -R $ref -T GenotypeGVCFs --variant ${file/.bam/.raw.g.vcf} -o ${file/.bam/.genotype.vcf} -ploidy 2 ;

  #filter
  java -Xmx1g -jar $GATK  -T VariantFiltration  -R $ref  -V ${file/.bam/.genotype.vcf} \
      --filterExpression "DP<3 || QD < 1.0 || FS > 200.0 || MQ < 20.0 || ReadPosRadPosRankSum < -20.0 || SOR > 10.0" \
      --filterName "SOFTFILT" \
      -o ${file/.bam/filt.soft}.vcf

      grep -P "#|PASS" ${file/.bam/filt.soft}.vcf > ${file/.bam/filt.soft}.pass.vcf


      $bgzip -f ${file/.bam/filt.soft}.vcf
      $BCFtools index -f ${file/.bam/filt.soft}.vcf.gz


done

echo "Done!"
