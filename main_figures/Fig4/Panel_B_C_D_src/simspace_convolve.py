import warnings
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

kernel_size = (5, 7, 10, 15)
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
                f'{script_dir}/../Panel_B_C_D_data/Xenium_reference_count.csv',
                f'{script_dir}/../Panel_B_C_D_data/Xenium_reference_metadata.csv',
                'Cluster',
                'x_centroid',
                'y_centroid',
                seed=seed,
            )
            for kernel in kernel_size:
                kernels = (kernel, kernel)
                spot_meta, spot_omics = ss.util.convolve(sim, kernel=kernels, scale=1, conv_type = 'sum')

                ## Uncomment to save the output files
                # spot_meta.to_csv(f"{script_dir}/spot_meta_kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv", index=False)
                # spot_omics.to_csv(f"{script_dir}/spot_omics_kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv", index=False)

kernel_size = (5, 7, 10, 15)
n_niche = (3, )
n_state = (11, )
seeds = (0, 1, 3, 4, )
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
                f'{script_dir}/../Panel_B_C_D_data/Xenium_reference_count_2.csv',
                f'{script_dir}/../Panel_B_C_D_data/Xenium_reference_metadata_2.csv',
                'Cluster',
                'x_centroid',
                'y_centroid',
                seed=seed,
            )
            for kernel in kernel_size:
                kernels = (kernel, kernel)
                spot_meta, spot_omics = ss.util.convolve(sim, kernel=kernels, scale=1, conv_type = 'sum')

                ## Uncomment to save the output files
                # spot_meta.to_csv(f"{script_dir}/spot_meta_kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv", index=False)
                # spot_omics.to_csv(f"{script_dir}/spot_omics_kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv", index=False)

kernel_size = (5, 7, 10, 15)
n_niche = (3, )
n_state = (14, )
seeds = (0, 1, 2, 3)
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
                f'{script_dir}/../Panel_B_C_D_data/Xenium_reference_count_3.csv',
                f'{script_dir}/../Panel_B_C_D_data/Xenium_reference_metadata_3.csv',
                'Cluster',
                'x_centroid',
                'y_centroid',
                seed=seed,
            )
            for kernel in kernel_size:
                kernels = (kernel, kernel)
                spot_meta, spot_omics = ss.util.convolve(sim, kernel=kernels, scale=1, conv_type = 'sum')

                ## Uncomment to save the output files
                # spot_meta.to_csv(f"{script_dir}/spot_meta_kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv", index=False)
                # spot_omics.to_csv(f"{script_dir}/spot_omics_kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv", index=False)
