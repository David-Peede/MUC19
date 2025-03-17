# `hmmix_tracts`

This directory contains the code to replicate the archaic tracts data used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`hmmix`](https://github.com/LauritsSkov/Introgression-detection)

## Data

All data is publicly available and can be downloaded from the following locations:

- [EPO Ancestral Allele Calls](http://ftp.ensembl.org/pub/release-74/fasta/ancestral_alleles/)
- [Altai Denisovan Data](http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Denisova/)
- [Altai Neanderthal Data](http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Altai/)
- [Chagyrskaya Neandertal Data](http://ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/)
- [Vindija Neanderthal Data](http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/)
- [Simons Genome Diversity Project Phased Dataset](https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/)
- [Phase 3 Release of the 1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)

## Code

The two shell scripts—`tgp_hmmix_tutorial.sh` and `papuan_hmmix_tutorial.sh`—provide comprehensive, step-by-step tutorials for applying `hmmix` to individuals from the 1000 Genomes Project and the Papuan population from the Simons Genome Diversity Project. The unfiltered decoded and annotated tracts for all autosomes can be found on Zenodo. Additionally, the pre-processed archaic tract information used as inputs for downstream analyses in the paper is also available on Zenodo:

- Chromosome 12's archaic tracts for the 1000 Genomes Project individuals at `hmmix_tracts/tgp_hmmix_haplotype_tracts_chr12.csv.gz`.
- Archaic tracts overlapping the _MUC19_ region for the 1000 Genomes Project individuals at `hmmix_tracts/tgp_hmmix_haplotype_tracts_muc19.csv.gz`.
- Autosomal archaic tracts for the Papuan individuals from the Simons Genome Diversity Project at `hmmix_tracts/papuan_hmmix_haplotype_tracts.csv.gz`.
- Coordinates for the putative Denisovan tracts for the Papuan individuals from the Simons Genome Diversity Project at `hmmix_tracts/den_tracts_in_pap.tar.gz`.

