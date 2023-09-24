# `analyses_nbs`

This directory contains the code to replicate the analyses used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`matpltlib v3.5.2`](https://matplotlib.org/)
- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)
- [`SciPy v1.8.0`](https://docs.scipy.org/doc/scipy/)

## Code

__Create the `meta_data`, `dataframes`, `supp_figures`, and `supp_tables` directories since GitHub only allows for 100Mb of storage.__

```bash
# Extract the meta_data directory.
tar -xf meta_data.tar.gz
# Extract the dataframes directory.
tar -xf dataframes.tar.gz
# Extract the supp_figures directory.
tar -xf supp_figures.tar.gz
# Extract the supp_tables directory.
tar -xf supp_tables.tar.gz
```

__Run the code necessary to generate any other input files, which can be found in the `README.md` of each directory.__