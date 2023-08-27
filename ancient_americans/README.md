# `ancient_americans`

This directory contains the code to replicate the ancient American results for ONLY our region of interest used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`Samtools v1.9`](http://www.htslib.org/doc/samtools.html)
- [`ANGSD v0.92`](http://www.popgen.dk/angsd/index.php/ANGSD)

## Data

All data is publicly available after signing data agreements and can be downloaded from the following publications:

- [Scheib et al., 2018](https://www.science.org/doi/10.1126/science.aar6851)
- [Lindo et al., 2018](https://www.science.org/doi/10.1126/sciadv.aau4921)
- [de la Fuente et al., 2018](https://doi.org/10.1073/pnas.1715688115)
- [Moreno-Mayar et al., 2018](https://doi.org/10.1038/nature25173)
- [Rasmussen et al., 2014](https://doi.org/10.1038/nature13025)
- [Villa-Islas et al., 2023](https://www.science.org/doi/10.1126/science.add6142)

## Code

__Extract the region of interest from the individual bam files.__

```bash
# Subset the focal region from the ancient American bam files.
cat ../meta_data/adna_samps.txt | while read ANC; do
samtools view -b -o ${ANC}_slc2a13_cntn1.bam ${ANC}_sorted_rmdup_autosomes.bam 12:40148827-41466217
done
```

__Index the individual bam files.__

```bash
# Index the ancient American bam files.
cat ../meta_data/adna_samps.txt | while read ANC; do
samtools index ${ANC}_slc2a13_cntn1.bam
done
```

__Compute the per-site allelic depth from high quality reads for each individual.__

```bash
# Compute depth of coverage per individual and per-site as well as depth per-base with a minimum quality score of 30 and a minimum of 1 individual per site.
cat ../meta_data/adna_samps.txt | while read ANC; do
angsd -i ${ANC}_slc2a13_cntn1.bam -doCounts 1 -minQ 30 -minInd 1 -dumpCounts 3 -doDepth 1 -out ${ANC}_ac
done
```

__Merge the outputs of `ANGSD` for analyses.__

```bash
# Combine the positions and allele information.
gunzip *.gz
cat ../meta_data/adna_samps.txt | while read ANC; do
paste -d'\t' ${ANC}_ac.pos ${ANC}_ac.counts > ${ANC}_ac.txt
done
```

