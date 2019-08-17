#!/bin/bash
#SBATCH -A phillipslab
#SBATCH --partition=short	### Partition
#SBATCH --job-name=cons_step3    ### Job Name
#SBATCH --time=1:00:00        ### WallTime
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --ntasks=1             ### Number of tasks per array job
#SBATCH --array=0-325           ### Array index
#SBATCH --mail-type=END,FAIL              ### Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=my@gmail.com  ### Where to send mail
#SBATCH --cpus-per-task=1            ### Number of CPU cores per task

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#load modules
#tools: bedops, samtools, bedtools, bcftools (version 1.7)
module load racs-eb BEDOPS/2.4.26 bedtools samtools bcftools

ref="/projects/phillipslab/ateterina/CeNDR/ref_245/c_elegans.PRJNA13758.WS245.genomic.fa"
region="I:11953012-11962484"
BCFtools="/projects/phillipslab/ateterina/scripts/bcftools-1.7/bcftools"

cd BAMs/hsf1;

list=(*.hsf1.bam);
file=${list[$SLURM_ARRAY_TASK_ID]};

#low-coverage sites
samtools depth -a -r $region ${file/.bam/.re.bam}| awk '{ if ($3 < 3) { print $1 "\t" $2-1 "\t" $2 }}' - > ${file/.bam/.DP3.mask.bed}

#extract high-quality indels
zcat ${file/.bam/filt.soft}.vcf.gz |vcf2bed --deletions - > ${file/.bam/.indel.bed}
zcat ${file/.bam/filt.soft}.vcf.gz |vcf2bed --insertions - >> ${file/.bam/.indel.bed}


#remove these indels from the mask
bedtools merge -i ${file/.bam/.DP3.mask.bed} | bedtools subtract -a - -b ${file/.bam/.indel.s.bed}  > ${file/.bam/.DP3.mask.noindel.bed}


#call consensus + IUPAC codes (works only with bcftools version 1.7)
name=${file/.hsf1.bam/}
samtools faidx $ref $region | $BCFtools consensus -I -s $name -m ${file/.bam/.DP3.mask.noindel.bed} ${file/.bam/filt.soft}.vcf.gz > ${file/.bam/.0.fasta}

#edit names
sed -i "s/I/$name-hsf1-I/" ${file/.bam/.0.fasta}


#clean up
#rm *.bed *.vcf *.gz *.re.* *.csi *.intervals

#combine all sequences
#for i in *0.fasta; do cat $i >> CeNDR_hsf-1_I_11953012-11962484_326_strains.fasta; done
#sed -i "s/:/-/g" CeNDR_hsf-1_I_11953012-11962484_326_strains.fasta
