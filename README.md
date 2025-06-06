# The *MUC19* Gene: An Evolutionary History of Recurrent Introgression and Natural Selection

Code associated with Villanea and Peede et. al. 202X.

## Directory Overview

```bash
├── amr_lai
│   ├── anc_props
│   └── region_beds
├── analyses_nbs
│   ├── dataframes
│   └── supp_tables
├── ancient_americans
├── annotations
├── arc_snp_density
├── figure_nbs
├── heterozygosity
├── hg38_data
├── hmmix_tracts
├── iHS
├── introgression
├── meta_data
├── pbs
├── pbs_sims
│   ├── data
│   └── scripts
├── psuedo_ancestry_painting
├── sequence_divergence
├── vcf_data
│   ├── ann_summary
│   ├── muc19
│   └── phasing
├── vntr
└── windowing
```

### Data Generation Directories

- `amr_lai`
- `ancient_americans`
- `annotations`
- `arc_snp_density`
- `heterozygosity`
- `hg38_data`
- `hmmix_tracts`
- `iHS`
- `introgression`
- `pbs`
- `pbs_sims`
- `psuedo_ancestry_painting`
- `sequence_divergence`
- `vcf_data`
- `vntr`
- `windowing`

### Data Analysis Directories

- `analyses_nbs`
- `figure_nbs`
- `meta_data`

### Notes

This repository contains all code, meta information, and final results. Key datasets of interest have been uploaded to Zenodo (see `./ZENODO_README.md`). All of the intermediary data used to generate the final set of results in `./analyses_nbs/dataframes` have been uploaded to DRYAD (see `./DRYAD_README.md`). As this code is associated with a manuscript that is currently going through the peer-review process, and due to the fact that GitHub only stores up to 100Mb per repository, some directories are compressed using `.tar.gz` or `.zip`. To extract the files from those directories, you will need to run either `tar -xf target_directory.tar.gz` or `unzip target_directory.zip`.