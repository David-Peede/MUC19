# `vntr`

This directory contains the code to replicate the VNTR data used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`InVNTR v1.4`](https://github.com/ValdmanisLab/InVNTR)

## Data

All data is publicly available and can be downloaded from the following locations:

- [1000 Genomes Project](http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/)
- [Human Pangenome Reference Consortium](https://projects.ensembl.org/hprc/)
- [Human Genome Structural Variant Consortium](https://www.internationalgenome.org/data-portal/data-collection/hgsvc2)

## `InVNTR` Preprocessing

Following the best practices outlined in [Course et. al., 2021](10.1101/gr.275560.121) count and extract reads for the focal region `Chr12:40482139-40491565` and two control regions `Chr7:5500000-5600000` and `Chr12:6490000-6590000` using the `InVNTR` program.

## `InVNTR` Postprocessing

In order to ensure that the repeat has been properly split into individual motifs, especially in repeats with high internal variability like mucins, it is important to ensure that all of the most common motifs in the `VNTR_frequency.csv` are the expected length. If motifs are too long and they include multiple repeat units in a single motif, we use find and replace to add comma's in the appropriate location in `VNTR.csv`. If the frame of the motif was shifted, or if the motifs were smaller than expected due to the delimiter appearing within the motif, we used find and replace to remove the commas.

### Output Overview

- `VNTR.csv`
  - CSV with each allele named after the file, cut into simple motifs based on the length or delimiter provided from the beginning of the allele.
- `VNTR_allele_length.csv`
  - CSV with the length of each allele accross the dataset.
- `VNTR_frequency.csv`
  - CSV with the frequency of each repeat motif accross the entire dataset.
- `VNTR_alleles.txt `
  - Text file with each allele after the name of the file.
- `VNTR_errors.txt`
  - Text file with any errors (i.e., `filename has only reverse end` or `filename has no sign of tandem repeat in forward, reverse`).