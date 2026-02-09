import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

import simspace as ss

########### Panel A:
shape = (150, 150)
num_iteration = 4
n_iter = 6
custom_neighbor = ss.spatial.generate_offsets(3, 'manhattan')

# param = ss.util.generate_random_parameters(n_group=2, n_state=9, seed=0)
# sim = ss.util.sim_from_params(
#     param,
#     shape=shape,
#     num_iteration=num_iteration,
#     n_iter=n_iter,
#     custom_neighbor=custom_neighbor,
#     seed=0
# )

# sim.fit_scdesign(
#     f'{script_dir}/Panel_A_B_C_data/Xenium_reference_count.csv',
#     f'{script_dir}/Panel_A_B_C_data/Xenium_reference_metadata.csv',
#     'Cluster',
#     'x_centroid',
#     'y_centroid',
#     seed=0,
# )

# kernel_size = (5, 7, 10, 15)
# for i, kernel in enumerate(kernel_size):
#     spot_meta, _ = ss.util.convolve(sim, kernel=(kernel,kernel))

#     cmap = sns.color_palette('tab20', n_colors=sim.num_states)
#     state_names = spot_meta.columns[2:]
#     state_name_mapping = {i: name for i, name in enumerate(state_names)}
#     state_colors = {state_name_mapping[i]: cmap[i] for i in range(len(state_names))}

#     fig = plt.figure(figsize=(6, 4), dpi=300)
#     ax = fig.add_subplot()
#     ax.set_aspect('equal')
#     ax = sns.scatterplot(data=sim.meta, x='col', y='row', hue='fitted_celltype', s=8, palette=state_colors)
#     ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='fitted_celltype')
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     plt.tight_layout()
#     plt.savefig(f'{script_dir}/SFig5_panel_A1.png')

#     ss.plot.spatial_pie(sim, spot_meta, (kernel,kernel))
#     plt.savefig(f'{script_dir}/SFig5_panel_A{2+i}.png')
#     plt.close()


########### Panel B:
param = ss.util.generate_random_parameters(n_group=3, n_state=11, seed=1)
sim = ss.util.sim_from_params(
    param,
    shape=shape,
    num_iteration=num_iteration,
    n_iter=n_iter,
    custom_neighbor=custom_neighbor,
    seed=1
)

sim.fit_scdesign(
    f'{script_dir}/Panel_A_B_C_data/Xenium_reference_count_2.csv',
    f'{script_dir}/Panel_A_B_C_data/Xenium_reference_metadata_2.csv',
    'Cluster',
    'x_centroid',
    'y_centroid',
    seed=1,
)

kernel_size = (5, 7, 10, 15)
for i, kernel in enumerate(kernel_size):
    spot_meta, _ = ss.util.convolve(sim, kernel=(kernel,kernel))

    cmap = sns.color_palette('tab20', n_colors=sim.num_states)
    state_names = spot_meta.columns[2:]
    state_name_mapping = {i: name for i, name in enumerate(state_names)}
    state_colors = {state_name_mapping[i]: cmap[i] for i in range(len(state_names))}

    fig = plt.figure(figsize=(6, 4), dpi=300)
    ax = fig.add_subplot()
    ax.set_aspect('equal')
    ax = sns.scatterplot(data=sim.meta, x='col', y='row', hue='fitted_celltype', s=8, palette=state_colors)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='fitted_celltype')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.tight_layout()
    plt.savefig(f'{script_dir}/SFig5_panel_B1.png')

    ss.plot.spatial_pie(sim, spot_meta, (kernel,kernel))
    plt.savefig(f'{script_dir}/SFig5_panel_B{2+i}.png')
    plt.close()

########### Panel C:
param = ss.util.generate_random_parameters(n_group=3, n_state=14, seed=2)
sim = ss.util.sim_from_params(
    param,
    shape=shape,
    num_iteration=num_iteration,
    n_iter=n_iter,
    custom_neighbor=custom_neighbor,
    seed=2
)

sim.fit_scdesign(
    f'{script_dir}/Panel_A_B_C_data/Xenium_reference_count_3.csv',
    f'{script_dir}/Panel_A_B_C_data/Xenium_reference_metadata_3.csv',
    'Cluster',
    'x_centroid',
    'y_centroid',
    seed=2,
)

kernel_size = (5, 7, 10, 15)
for i, kernel in enumerate(kernel_size):
    spot_meta, _ = ss.util.convolve(sim, kernel=(kernel,kernel))

    cmap = sns.color_palette('tab20', n_colors=sim.num_states)
    state_names = spot_meta.columns[2:]
    state_name_mapping = {i: name for i, name in enumerate(state_names)}
    state_colors = {state_name_mapping[i]: cmap[i] for i in range(len(state_names))}

    fig = plt.figure(figsize=(6, 4), dpi=300)
    ax = fig.add_subplot()
    ax.set_aspect('equal')
    ax = sns.scatterplot(data=sim.meta, x='col', y='row', hue='fitted_celltype', s=8, palette=state_colors)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='fitted_celltype')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.tight_layout()
    plt.savefig(f'{script_dir}/SFig5_panel_C1.png')

    ss.plot.spatial_pie(sim, spot_meta, (kernel,kernel))
    plt.savefig(f'{script_dir}/SFig5_panel_C{2+i}.png')
    plt.close()