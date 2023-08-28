# `pbs`

This directory contains the code to replicate the empirical PBS results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Compute PBS partitioned by ancestry components.__

```bash
# Compute PBS (MXL:CHB:CEU) for every SNP.
for CHR in {1..22}; do
python mxl_pbs_chromosome.py ${CHR}
done
```

__Compute the windowed average PBS value.__

```bash
# Compute PBS (MXL:CHB:CEU) in non-overlapping windows.
for CHR in {1..22}; do
python mxl_pbs_windows.py ${CHR} 748
done
```



__Compute PBS for Denisovan-specific sites.__

```bash
# Compute PBS (MXL:CHB:CEU) for every Denisovan-specific SNP.
for CHR in {1..22}; do
python den_sites_pbs_chromosome.py ${CHR}
done
```

