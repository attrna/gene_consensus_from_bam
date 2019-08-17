#!/bin/bash

#clean up
rm *.bed *.vcf *.gz *.re.* *.csi *.intervals

#combine all sequences
for i in *0.fasta; do
  cat $i >> All_sequences.fasta;
done

sed -i "s/:/-/g" All_sequences.fasta
