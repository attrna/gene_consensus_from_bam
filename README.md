# **Generate a gene consensus sequence from an individual whole-genome sequencing**

Scripts to exctract gene/region sequences from the CeNDR database (https://www.elegansvariation.org, *Caenorhabditis elegans* Natural Diversity Resource)

## Required tools:
- samtools, bcftools (version 1.7 or more recent), htslib [download](http://www.htslib.org/download/)
- bedtools [download](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- GATK [download](https://software.broadinstitute.org/gatk/download/)
- Picard [download](https://broadinstitute.github.io/picard/)
- bedops [download](https://bedops.readthedocs.io/en/latest/)

## The pipeline:
* 0 step - load the data from CeNDR;
* 1 step - extract the region(s) of interest;
* 2 step - call and filter variants;
* 3 step - create a mask and call the consensus sequences;
* 4 step - clean up and combine the sequences.

Before use, edit the paths to the tools if necessary.
