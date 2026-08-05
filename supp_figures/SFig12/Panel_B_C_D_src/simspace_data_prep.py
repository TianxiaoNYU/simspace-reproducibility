import warnings

from SpatialSim.src.SpatialSim.offsets import generate_offsets
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

import simspace as ss

shape = (150, 150)
num_iteration = 4
n_iter = 6
custom_neighbor = ss.spatial.generate_offsets(3, 'manhattan')

n_niche = (2, 3)
n_state = (8, 9)
seeds = (0, 1)

for niches in n_niche:
    for states in n_state:
        for seed in seeds:
            param = ss.util.generate_random_parameters(n_group=niches, n_state=states, seed=seed)
            sim = ss.util.sim_from_params(
                param,
                shape=shape,
                num_iteration=num_iteration,
                n_iter=n_iter,
                custom_neighbor=custom_neighbor,
                seed=seed
            )

            sim.fit_scdesign(
                f'{script_dir}/Xenium_reference_count.csv',
                f'{script_dir}/Xenium_reference_metadata.csv',
                'Cluster',
                'x_centroid',
                'y_centroid',
                seed=seed,
            )

            # Uncomment to save simulated data
            # sim.meta.to_csv(f"{script_dir}/benchmark_datasets/sim_meta_niche{niches}_state{states}_seed{seed}.csv", index=False)
            # sim.omics.to_csv(f"{script_dir}/benchmark_datasets/sim_omics_niche{niches}_state{states}_seed{seed}.csv", index=False)