# `heterozygosity`

This directory contains the code to replicate the heterozygosity results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute the number of heterozygous sites amongst the archaics.__

```bash
# Compute the number of heterozygous sites per archaic in non-overlapping windows.
for CHR in {1..22}; do
python arc_heterozygosity_windows.py ${CHR} 72
done
```

__Compute the number of heterozygous sites amongst individuals in the combined dataset.__

```bash
# Compute the number of heterozygous sites per individual in non-overlapping windows.
for CHR in {1..22}; do
python tgp_heterozygosity_windows.py ${CHR} 72
done
```

