# `arc_snp_density`

This directory contains the code to replicate the archaic SNP density results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute archaic SNP density in non-overlapping windows.__

```bash
# Compute the archaic SNP density in non-overlapping 72kb windows.
for CHR in {1..22}; do
python tgp_arc_snp_density_windows.py ${CHR} 72
done
```

