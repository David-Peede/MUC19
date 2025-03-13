# `vntr`

This directory contains the code to replicate the VNTR data used in this study.

## Packages

All software is publicly available and their documentation can be viewed at the following places:

- [`InVNTR v1.4`](https://github.com/ValdmanisLab/InVNTR)
- [`RepeatLengthEstimatorFromShortRead`](https://github.com/ValdmanisLab/RepeatLengthEstimatorFromShortRead)


## Data

All data is publicly available and can be downloaded from the following locations:

- [1000 Genomes Project](http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/)
- [Human Pangenome Reference Consortium](https://projects.ensembl.org/hprc/)
- [Human Genome Structural Variant Consortium](https://www.internationalgenome.org/data-portal/data-collection/hgsvc2)

## Long Read Data

Long read data was processed with `InVNTR`, which extracted sequence information from haplotype assemblies corresponding with `Chr12:40482139-40491565` on hg38, and attempted to decompose into individual repeat motifs. We followed up by manually ensuring that the VNTR decomposed the motifs by correcting motifs that were not properly split into 30bp motifs in VNTR_frequency.csv. This was done by using find and replace in a text editor to add commas to a transposed VNTR.csv file output.

## Short Read Data

Short read data was processed with `RepeatLengthEstimatorFromShortRead`. This bash script follows the best practices outlined in [Course et. al., 2021](10.1101/gr.275560.121) to count and extract reads for the focal region in hg38 `Chr12:40482139-40491565` and two control regions `Chr7:5500000-5600000` and `Chr12:6490000-6590000`. We ran this script in a folder with short read whole genome sequences files. It counts the number of reads in each of the 3 regions and divides the number of reads by the number of base pairs in each region, averages the two control regions, and divides the focal region value by the control region value to come up with an enrichment value for each sample. This enrichment value was then multiplied by the number of base pairs in the focal region in hg38 as well as the number of repeat copies in hg38 in order to estimate the length of the repeat for each short read whole genome and the number of repeat copies. 
