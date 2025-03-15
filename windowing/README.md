# `windowing`

This directory contains the code to replicate the non-overlapping windows used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Code

__QC non-overlapping windows across autosomes for the single archaic datasets.__
```bash
for CHR in {1..22}; do for PREFIX in den_masked_no_aa alt_masked_no_aa cha_masked_no_aa vin_masked_no_aa; do
python single_archaic_window_qc_v_revisions..py 72 ${CHR} ${PREFIX}
done; done
```


__QC non-overlapping windows across autosomes for the combined archaic datasets.__
```bash
for CHR in {1..22}; do for PREFIX in arcs_masked_aa arcs_masked_no_aa; do
python all_archaics_window_qc_v_revisions.py 72 ${CHR} ${PREFIX}
done; done
```


__QC non-overlapping windows across autosomes for 1000 Genomes Project datasets.__
```bash
for CHR in {1..22}; do for WIND in 72 742; do for PREFIX in tgp_mod_aa tgp_mod_no_aa; do
python tgp_window_qc_v_revisions.py ${WIND} ${CHR} ${PREFIX}
done; done; done
```


__QC non-overlapping windows across autosomes for the single archaic + modern human datasets.__
```bash
for CHR in {1..22}; do for WIND in 72 742; do for PREFIX in tgp_den_masked_aa tgp_den_masked_no_aa tgp_alt_masked_aa tgp_alt_masked_no_aa tgp_cha_masked_aa tgp_cha_masked_no_aa tgp_vin_masked_aa tgp_vin_masked_no_aa; do
python tgp_sgdp_single_archaic_window_qc_v_revisions.py ${WIND} ${CHR} ${PREFIX}
done; done; done
```


__QC non-overlapping windows across autosomes for the all archaics + modern human combined datasets.__
```bash
for CHR in {1..22}; do for WIND in 72 742 670; do for PREFIX in tgp_arcs_masked_no_aa tgp_arcs_masked_aa; do
python tgp_sgdp_all_archaics_window_qc_v_revisions.py ${WIND} ${CHR} ${PREFIX}
done; done; done
```


**Consolidate the non-overlapping window QC information.**
```bash
for PREFIX in den_masked_no_aa alt_masked_no_aa cha_masked_no_aa vin_masked_no_aa arcs_masked_aa arcs_masked_no_aa tgp_mod_aa tgp_mod_no_aa tgp_den_masked_aa tgp_den_masked_no_aa tgp_alt_masked_aa tgp_alt_masked_no_aa tgp_cha_masked_aa tgp_cha_masked_no_aa tgp_vin_masked_aa tgp_vin_masked_no_aa tgp_arcs_masked_no_aa tgp_arcs_masked_aa; do
python consolidate_all_archaics_tgp_sgdp_windows_v_revisions.py 72 ${PREFIX}
done
for PREFIX in tgp_mod_aa tgp_mod_no_aa tgp_den_masked_aa tgp_den_masked_no_aa tgp_alt_masked_aa tgp_alt_masked_no_aa tgp_cha_masked_aa tgp_cha_masked_no_aa tgp_vin_masked_aa tgp_vin_masked_no_aa tgp_arcs_masked_no_aa tgp_arcs_masked_aa; do
python consolidate_all_archaics_tgp_sgdp_windows_v_revisions.py 742 ${PREFIX}
done
for PREFIX in tgp_arcs_masked_no_aa tgp_arcs_masked_aa; do
python consolidate_all_archaics_tgp_sgdp_windows_v_revisions.py 670 ${PREFIX}
done
```


**‌Compute the effective sequence lengths for the focal regions.**
```bash
python compute_region_effective_sequence_lengths_v_revisions.py 12 40759001 40831000 72kb
python compute_region_effective_sequence_lengths_v_revisions.py 12 40272001 41014000 742kb
```


**‌QC the Denisovan introgressed tracts found in Papuans from the Simons Genome Diversity Project across the autosomes for all Papuans.**
```bash
for CHR in {1..22}; do for PAP in S_Papuan-1.DG S_Papuan-2.DG S_Papuan-3.DG S_Papuan-4.DG S_Papuan-5.DG S_Papuan-6.DG S_Papuan-7.DG S_Papuan-8.DG S_Papuan-9.DG S_Papuan-10.DG S_Papuan-11.DG S_Papuan-12.DG S_Papuan-13.DG S_Papuan-14.DG B_Papuan-15.DG; do
python sgdp_denisovan_intro_tracts_in_papuans_qc_v_revisions.py ${PAP} ${CHR} sgdp_den_masked_no_aa
done; done
```


**Consolidate the Denisovan introgressed tracts found in Papuans from the Simons Genome Diversity Project QC information.**
```bash
for PAP in S_Papuan-1.DG S_Papuan-2.DG S_Papuan-3.DG S_Papuan-4.DG S_Papuan-5.DG S_Papuan-6.DG S_Papuan-7.DG S_Papuan-8.DG S_Papuan-9.DG S_Papuan-10.DG S_Papuan-11.DG S_Papuan-12.DG S_Papuan-13.DG S_Papuan-14.DG B_Papuan-15.DG; do
python consolidate_sgdp_denisovan_intro_tracts_in_papuans_v_revisions.py ${PAP} sgdp_den_masked_no_aa
done
```

