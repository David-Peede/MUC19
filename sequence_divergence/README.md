# `sequence_divergence`

This directory contains the code to replicate the sequence divergence results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute sequence divergence between haplotypes in the 1000 Genomes Project and the archaic diplotypes.__

```bash
# Compute sequence divergence between modern haplotypes and archaic diplotypes in non-overlapping windows.
for CHR in {1..22}; do
python tgp_haplotype_arc_diplotype_sequence_divergence_windows.py ${CHR} 748
done
```

