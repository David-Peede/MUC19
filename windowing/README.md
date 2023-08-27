# `windowing`

This directory contains the code to replicate the non-overlapping windows used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__Create non-overlapping window quality control reports.__

```bash
# QC non-overlapping windows across the autosomes.
for CHR in {1..22}; do for WIND in 72 748; do
python arc_window_qc.py ${WIND} ${CHR}
python tgp_window_qc.py ${WIND} ${CHR}
done; done
```

__Consolidate the quality control reports for analyses.__

```bash
# Consolidate the non-overlapping window information.
for WIND in 72 748; do
python consolidate_nonoverlapping_windows.py ${WIND} arc
python consolidate_nonoverlapping_windows.py ${WIND} tgp
done
```

__Compute the effective sequence lengths for the two regions of interest.__

```bash
# Determine the effective sequence lengths for the focal 72kb region.
python compute_esl.py 12 40758000 40830000 72kb
python compute_sgdp_arc_72kb_esl.py
# Determine the effective sequence lengths for the 748kb longest inferred introgressed tract in MXL.
python compute_esl.py 12 40269000 41017000 748kb
```

