# `sup_intro_sims`

This directory contains the code to replicate the super-archaic introgression simulation results used in this study.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`NumPy v1.22.3`](https://numpy.org/doc/stable/reference/index.html)
- [`msprime v1.1.1`](https://tskit.dev/msprime/docs/stable/intro.html#)

## Code

__Run 1,000 replicate simulations per pairwise combination of the time of introgression, time of super-archaic divergence, Denisovan effective population size, and Altai Neanderthal effective population size.__

```bash
# Run simulations based on the Argweaver-D model of human evolution.
for TGF in 2.5 3.0; do for DIV in 1.0 1.5 1.75 2.0 2.25 2.5 2.75 3.0; do for DEN in 2500 3400; do for NEA in 3400 2500; do
python sup_intro_model.py ${TGF} ${DIV} ${NEA} ${DEN}
done; done; done; done
```

