# `annotations`

This directory contains the code to replicate the per-NCBI RefSeq Select genes information used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__QC NCBI RefSeq Select genes across autosomes for the TGP + Denisovan combined dataset.__
```bash
for CHR in {1..22}; do
python tgp_single_archaic_gene_qc_v_revisions.py ${CHR} tgp_den_masked_no_aa
done
```

**Consolidate the NCBI RefSeq Select gene information.**
```bash
python consolidate_tgp_single_archaic_genes_v_revisions.py tgp_den_masked_no_aa
```