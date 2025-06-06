# `pbs_sims`

This directory contains the code to replicate the PBS simulations.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`SLiM v4.1`](https://messerlab.org/slim/)
- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas v1.4.2`](https://pandas.pydata.org/docs/)
- [`scikit-allel v1.3.5`](https://scikit-allel.readthedocs.io/en/stable/index.html)

## Overview

- `data`
  - _MUC19_ annotation inputs for `SLiM`.
- `scripts`
  - `.slim` files correspond to the SLiMulation scripts.
	  - You can execute any of the `.slim` files like so: `slim -d nrun=${REP} {slim_script_prefix}.slim` where `${REP}` is an arbitrary replicate number, NOT the random seed.
  - `.py` files correspond to how we extracted our results from the SLiMulations.
	  - Note that based on our conditioning on introgressed alleles being present in the sampled individuals (see our supplementary material), I had to manually check the SLiMulation output to confirm that it was a valid replicate for our purposes, which is further exacerbated when also conditioning on the positively selected allele being of archaic origin. Unfortunately, I could not figure out the optimal way to do this that, but here is what I did at the end of the day:
		  - Step 1) After simulating approximately a thousand replicates, I generated QC reports to determine how many replicates I could currently retain.
			  - Neutral/Negative SLiMulations: `for MODEL in neutral negative; do python compile_slimulations_qc_summary_v_revisions.py ${MODEL} recomb_map VCF_s${MODEL}; done`
			  - Positive SLiMulations: `for SEL_COEFF in 1 01 0015; do python compile_positive_slimulations_qc_summary_v_revisions.py ${SEL_COEFF}; done`
		  - Step 2) Manually check the QC files to see if I have reached my target number of SLiMulation replicates per selection scenario. Once I reached my target number of SLiMulations, I would advance to Step 3. Otherwise, I would run another batch of SLiMulations based on the number of additional replicates needed to reach my target, and then redo Step 1.
			  - Neutral/Negative SLiMulations: 10k replicates.
			  - Positive SLiMulations: 1k replicates.
		  - Step 3) Extract the information needed from the QC'ed replicates for hypothesis testing.
			  - Neutral/Negative SLiMulations: `for MODEL in neutral negative; do python extract_neutral_negative_recomb_map_slimulated_results_v_revisions.py ${MODEL}; done`
			  - Positive SLiMulations: `for SEL_COEFF in 1 01 0015; do python extract_positive_slimulations_info_v_revisions.py ${SEL_COEFF}; done`
