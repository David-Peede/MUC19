# Import packages.
import math
import msprime
import numpy as np
import sys

### sys.argv[1] = time of gene flow ###
### sys.argv[2] = divergence time for the super archaic population ###
### sys.argv[3] = Ne_{NEA} ###
### sys.argv[4] = Ne_{DEN} ###


# Define a model based on the Argweaver-D paper.
def hubisz_demo(f_sup, tgf, t_sup_div, ne_nea, ne_den):
    # Times are provided in years, so we convert into generations.
    generation_time = 29
    # Intialize sampling times for the Altai Nean. and the Denisovan based
    # on the ages from Prufer et al 2017.
    nea_samp_time = 122e3 / 29
    den_samp_time = 72e3 / 29
    # Intialize a demographic model.
    model = msprime.Demography()
    # AFR.
    model.add_population(name='AFR', initial_size=23700)
    # NEA.
    model.add_population(name='NEA', initial_size=ne_nea, default_sampling_time=nea_samp_time)
    # DEN.
    model.add_population(name='DEN', initial_size=ne_den, default_sampling_time=den_samp_time)
    # ARC.
    model.add_population(name='ARC', initial_size=7100)
    # SUP.
    model.add_population(name='SUP', initial_size=2000)
    # HOM.
    model.add_population(name='HOM', initial_size=18500)
    # ANC.
    model.add_population(name='ANC', initial_size=18500)
    # Introgression from SUP to DEN.
    model.add_mass_migration(time=(tgf / generation_time), source='DEN', dest='SUP', proportion=f_sup)
    # DEN and NEA originate from ARC.
    model.add_population_split(time=(415e3 / generation_time), derived=['NEA', 'DEN'], ancestral='ARC')
    # ARC and AFR originate from HOM.
    model.add_population_split(time=(575e3 / generation_time), derived=['ARC', 'AFR'], ancestral='HOM')
    # SUP and HOM originat from ANC.
    model.add_population_split(time=(t_sup_div / generation_time), derived=['HOM', 'SUP'], ancestral='ANC')
    return model

# Read in the sys args.
sup2den = 1.0
t_intro = float(sys.argv[1]) * 100_000
t_div = float(sys.argv[2]) * 1_000_000
nea_ne = int(sys.argv[3])
den_ne = int(sys.argv[4])
# Intialize the results matrix.
results_matrix = np.empty((1_000, 2))

# For 1000 simulation replicates...
for rep in range(1000):
    # Run an ancestry simulation.
    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(504, ploidy=2, population='AFR'),
            msprime.SampleSet(1, ploidy=2, population='NEA'),
            msprime.SampleSet(1, ploidy=2, population='DEN'),
        ],
        demography=hubisz_demo(sup2den, t_intro, t_div, nea_ne, den_ne),
        sequence_length=72_000,
        recombination_rate=6.986e-9,
        random_seed=rep+1,
    )
    # Overlay mutations.
    mts = msprime.sim_mutations(
        tree_sequence=ts, rate=1.429e-8,
        model='jc69', random_seed=rep+1,
        discrete_genome=False,
    )
    # Extract the genotype matrix.
    genotype_matrix = mts.genotype_matrix()
    # Extract the allele frequencies for each population.
    afr_freqs = np.sum(genotype_matrix[:, 0:(504*2)], axis=1) / genotype_matrix[:, 0:(504*2)].shape[1]
    nea_freqs = np.sum(genotype_matrix[:, (504*2):(504*2)+2], axis=1) / genotype_matrix[:, (504*2):(504*2)+2].shape[1]
    den_freqs = np.sum(genotype_matrix[:, (504*2)+2:], axis=1) / genotype_matrix[:, (504*2)+2:].shape[1]
    # Determine the number of shared derived and shared ancestral sites.
    der_div = np.where((afr_freqs >= 0.95) & (nea_freqs == 1.0) & (den_freqs == 0.0))[0].size
    anc_div = np.where((afr_freqs <= 0.05) & (nea_freqs == 0.0) & (den_freqs == 1.0))[0].size
    # Fill the results matrix.
    results_matrix[rep, :] = np.array([der_div, anc_div])

# Export the results matrix.
np.savetxt(
    f'./simulated_results/sup2den_f{sup2den}_tgf_{t_intro}kya_tdiv_{t_div}mya_neaNe_{nea_ne}_denNe_{den_ne}.csv.gz',
    results_matrix, fmt='%d', delimiter=',', newline='\n',
)
