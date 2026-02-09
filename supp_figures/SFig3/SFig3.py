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
def ss_statistics_phi(
        shape = (100, 100), 
        n_group = 3, 
        n_state = 9, 
        num_iteration = 4, 
        n_iter = 6, 
        custom_neighbor = ss.spatial.generate_offsets(3, 'manhattan'),
        seed = 0,   
        phi = 5,
    ):
    param_1 = ss.util.generate_random_parameters(n_group=n_group, n_state=n_state, seed=seed)
    param_1['phi_replicates'] = phi

    sim1 = ss.util.sim_from_params(
        param_1,
        shape=shape,
        num_iteration=num_iteration,
        n_iter=n_iter,
        custom_neighbor=custom_neighbor,
        seed=seed
    )  

    mi_sim = ss.spatial.integrate_morans_I(sim1.meta['state'], sim1.meta[['row', 'col']], sim1.meta['state'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(sim1.meta['state'], sim1.meta[['row', 'col']], sim1.meta['state'].unique().tolist())
    ci_sim = ss.spatial.calculate_interaction_score(sim1.meta['state'], sim1.meta[['row', 'col']], sim1.meta['state'].unique().tolist())

    return [mi_sim, gc_sim, ci_sim]

phi_list = [1, 2, 3, 4, 4.5, 5, 6, 7, 8]
sim_phi_mi = [[] for _ in range(9)]
sim_phi_gc = [[] for _ in range(9)]
sim_phi_ci = [[] for _ in range(9)]

for i in range(9):
    for j in range(8):
        ## SimSpace
        sim_stats = ss_statistics_phi(
            shape=(100, 100), n_group=2, n_state=9, num_iteration=4, n_iter=6, 
            custom_neighbor=ss.spatial.generate_offsets(3, 'manhattan'), seed=j,
            phi=phi_list[i])
        sim_phi_mi[i].append(sim_stats[0])
        sim_phi_gc[i].append(sim_stats[1])
        sim_phi_ci[i].append(sim_stats[2].values.reshape(-1).tolist())

sim_phi_mi = [np.concatenate(item) for item in sim_phi_mi]
sim_phi_gc = [np.concatenate(item) for item in sim_phi_gc]
sim_phi_ci = [np.concatenate(item) for item in sim_phi_ci]

Xenium_mi_sim = []
Xenium_gc_sim = []
Xenium_ci_sim = []
merfish_mi_sim = []
merfish_gc_sim = []
merfish_ci_sim = []

for i in range(8):
    ## Xenium
    tile_path = f'{script_dir}/Panel_A_B_C_data/Xenium/tile_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    cells['x_centroid'] = (cells['x_centroid'] - cells['x_centroid'].min()) / (cells['x_centroid'].max() - cells['x_centroid'].min()) * 100
    cells['y_centroid'] = (cells['y_centroid'] - cells['y_centroid'].min()) / (cells['y_centroid'].max() - cells['y_centroid'].min()) * 100
    mi_sim = ss.spatial.integrate_morans_I(cells['Cluster'], cells[['x_centroid', 'y_centroid']], cells['Cluster'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(cells['Cluster'], cells[['x_centroid', 'y_centroid']], cells['Cluster'].unique().tolist())
    ci_sim = ss.spatial.calculate_interaction_score(cells['Cluster'], cells[['x_centroid', 'y_centroid']], cells['Cluster'].unique().tolist())
    Xenium_mi_sim.append(mi_sim)
    Xenium_gc_sim.append(gc_sim)
    Xenium_ci_sim.append(ci_sim.values.reshape(-1).tolist())
    ## MERFISH
    tile_path = f'{script_dir}/Panel_A_B_C_data/MERFISH/merfish_data_sample_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    cells['Centroid_X'] = (cells['Centroid_X'] - cells['Centroid_X'].min()) / (cells['Centroid_X'].max() - cells['Centroid_X'].min()) * 100
    cells['Centroid_Y'] = (cells['Centroid_Y'] - cells['Centroid_Y'].min()) / (cells['Centroid_Y'].max() - cells['Centroid_Y'].min()) * 100
    mi_sim = ss.spatial.integrate_morans_I(cells['Cell_class'], cells[['Centroid_X', 'Centroid_Y']], cells['Cell_class'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(cells['Cell_class'], cells[['Centroid_X', 'Centroid_Y']], cells['Cell_class'].unique().tolist())
    ci_sim = ss.spatial.calculate_interaction_score(cells['Cell_class'], cells[['Centroid_X', 'Centroid_Y']], cells['Cell_class'].unique().tolist())
    merfish_mi_sim.append(mi_sim)
    merfish_gc_sim.append(gc_sim)
    merfish_ci_sim.append(ci_sim.values.reshape(-1).tolist())

Xenium_mi_sim = np.concatenate(Xenium_mi_sim)
Xenium_gc_sim = np.concatenate(Xenium_gc_sim)
Xenium_ci_sim = np.concatenate(Xenium_ci_sim)
merfish_mi_sim = np.concatenate(merfish_mi_sim)
merfish_gc_sim = np.concatenate(merfish_gc_sim)
merfish_ci_sim = np.concatenate(merfish_ci_sim)

plt_data_mi = [
    Xenium_mi_sim, merfish_mi_sim, sim_phi_mi[0], sim_phi_mi[1], 
    sim_phi_mi[2], sim_phi_mi[3], sim_phi_mi[4], sim_phi_mi[5], 
    sim_phi_mi[6], sim_phi_mi[7], sim_phi_mi[8], 
]

plt_data_gc = [
    Xenium_gc_sim, merfish_gc_sim, sim_phi_gc[0], sim_phi_gc[1], 
    sim_phi_gc[2], sim_phi_gc[3], sim_phi_gc[4], sim_phi_gc[5], 
    sim_phi_gc[6], sim_phi_gc[7], sim_phi_gc[8], 
]

plt_data_ci = [
    Xenium_ci_sim, merfish_ci_sim, sim_phi_ci[0], sim_phi_ci[1], 
    sim_phi_ci[2], sim_phi_ci[3], sim_phi_ci[4], sim_phi_ci[5], 
    sim_phi_ci[6], sim_phi_ci[7], sim_phi_ci[8], 
]

fig, axs = plt.subplots(1, 3, figsize=(10, 3.5), dpi=300)
palette = sns.color_palette("muted", 11)
palette = list(palette)
first = np.array(palette[0])
# blend the first color toward a light blue and set as the second color
palette[1] = tuple(np.clip(first * 0.4 + np.array([0.7, 0.85, 1.0]) * 0.6, 0, 1))
palette[2:] =  sns.color_palette("viridis", 9)

# Prepare data for violin plots
labels = [
    "Xenium", "MERFISH", "1", "2", 
    "3", "4", "4.5", "5", 
    "6", "7", "8", 
]

# Panel 1: Moran's I
sns.violinplot(data=plt_data_mi, ax=axs[0], palette=palette)
axs[0].set_xticklabels(labels, rotation=30, ha='center')
axs[0].set_xlabel("Phi for Spatial Smoothness")
axs[0].set_ylabel("Moran's I")
axs[0].set_title("Moran's I Distribution")

# Panel 2: Geary's C
sns.violinplot(data=plt_data_gc, ax=axs[1], palette=palette)
axs[1].set_xticklabels(labels, rotation=30, ha='center')
axs[1].set_xlabel("Phi for Spatial Smoothness")
axs[1].set_ylabel("Geary's C")
axs[1].set_title("Geary's C Distribution")

# Panel 3: Cell Type Interaction
sns.violinplot(data=plt_data_ci, ax=axs[2], palette=palette)
axs[2].set_xticklabels(labels, rotation=30, ha='center')
axs[2].set_xlabel("Phi for Spatial Smoothness")
axs[2].set_ylabel("Cell Type Interaction")
axs[2].set_title("Cell Type Interaction Distribution")

plt.tight_layout()
plt.savefig(f'{script_dir}/SFig3_panel_A1.png')

fig, axes = plt.subplots(3, 3, figsize=(9, 9), dpi=300)
axes = axes.flatten()

for i, phi in enumerate(phi_list):
    param = ss.util.generate_random_parameters(n_group=2, n_state=9, seed=0)
    param['phi_replicates'] = phi

    sim = ss.util.sim_from_params(
        param,
        shape=(100, 100),
        num_iteration=4,
        n_iter=6,
        custom_neighbor=ss.spatial.generate_offsets(3, 'manhattan'),
        seed=0
    )

    ax = axes[i]
    sns.scatterplot(
        x='row', y='col', data=sim.meta,
        hue='state', palette='tab20', legend=False, s=12, ax=ax
    )
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_title(f'SimSpace: Phi = {phi}')
    ax.set_xticks([])
    ax.set_yticks([])

plt.tight_layout()
plt.savefig(f'{script_dir}/SFig3_panel_A2.png')


########### Panel B:
def ss_statistics(
        shape = (100, 100), 
        n_group = 3, 
        n_state = 9, 
        num_iteration = 4, 
        n_iter = 6, 
        custom_neighbor = ss.spatial.generate_offsets(3, 'manhattan'),
        seed = 0,  
    ):
    param_1 = ss.util.generate_random_parameters(n_group=n_group, n_state=n_state, seed=seed)

    sim1 = ss.util.sim_from_params(
        param_1,
        shape=shape,
        num_iteration=num_iteration,
        n_iter=n_iter,
        custom_neighbor=custom_neighbor,
        seed=seed
    )  

    mi_sim = ss.spatial.integrate_morans_I(sim1.meta['state'], sim1.meta[['row', 'col']], sim1.meta['state'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(sim1.meta['state'], sim1.meta[['row', 'col']], sim1.meta['state'].unique().tolist())
    ci_sim = ss.spatial.calculate_interaction_score(sim1.meta['state'], sim1.meta[['row', 'col']], sim1.meta['state'].unique().tolist())

    return [mi_sim, gc_sim, ci_sim]

neighborhood_sizes = [2, 3, 4, 5, 6, 7]
sim_neighbor_mi = [[] for _ in range(6)]
sim_neighbor_gc = [[] for _ in range(6)]
sim_neighbor_ci = [[] for _ in range(6)]

for i, neighborhood_size in enumerate(neighborhood_sizes):
    for j in range(8):
        ## SimSpace
        sim_stats = ss_statistics(
            shape=(100, 100), n_group=2, n_state=9, num_iteration=4, n_iter=6, 
            custom_neighbor=ss.spatial.generate_offsets(neighborhood_size, 'manhattan'), seed=j)
        sim_neighbor_mi[i].append(sim_stats[0])
        sim_neighbor_gc[i].append(sim_stats[1])
        sim_neighbor_ci[i].append(sim_stats[2].values.reshape(-1).tolist())
sim_neighbor_mi = [np.concatenate(item) for item in sim_neighbor_mi]
sim_neighbor_gc = [np.concatenate(item) for item in sim_neighbor_gc]
sim_neighbor_ci = [np.concatenate(item) for item in sim_neighbor_ci]

plt_data_mi = [
    Xenium_mi_sim, merfish_mi_sim, sim_neighbor_mi[0], sim_neighbor_mi[1], 
    sim_neighbor_mi[2], sim_neighbor_mi[3], sim_neighbor_mi[4], sim_neighbor_mi[5], 
]

plt_data_gc = [
    Xenium_gc_sim, merfish_gc_sim, sim_neighbor_gc[0], sim_neighbor_gc[1], 
    sim_neighbor_gc[2], sim_neighbor_gc[3], sim_neighbor_gc[4], sim_neighbor_gc[5], 
]

plt_data_ci = [
    Xenium_ci_sim, merfish_ci_sim, sim_neighbor_ci[0], sim_neighbor_ci[1], 
    sim_neighbor_ci[2], sim_neighbor_ci[3], sim_neighbor_ci[4], sim_neighbor_ci[5], 
]

# Prepare data for violin plots
labels = [
    "Xenium", "MERFISH", "2", 
    "3", "4", "5", "6", "7", 
]

fig, axs = plt.subplots(1, 3, figsize=(10, 3.5), dpi=300)
palette = sns.color_palette("muted", 11)
palette = list(palette)
first = np.array(palette[0])
# blend the first color toward a light blue and set as the second color
palette[1] = tuple(np.clip(first * 0.4 + np.array([0.7, 0.85, 1.0]) * 0.6, 0, 1))
palette[2:] =  sns.color_palette("viridis", 6)

# Panel 1: Moran's I
sns.violinplot(data=plt_data_mi, ax=axs[0], palette=palette)
axs[0].set_xticklabels(labels, rotation=30, ha='center')
axs[0].set_xlabel("Neighborhood Size (Manhattan Distance)")
axs[0].set_ylabel("Moran's I")
axs[0].set_title("Moran's I Distribution")

# Panel 2: Geary's C
sns.violinplot(data=plt_data_gc, ax=axs[1], palette=palette)
axs[1].set_xticklabels(labels, rotation=30, ha='center')
axs[1].set_xlabel("Neighborhood Size (Manhattan Distance)")
axs[1].set_ylabel("Geary's C")
axs[1].set_title("Geary's C Distribution")

# Panel 3: Cell Type Interaction
sns.violinplot(data=plt_data_ci, ax=axs[2], palette=palette)
axs[2].set_xticklabels(labels, rotation=30, ha='center')
axs[2].set_xlabel("Neighborhood Size (Manhattan Distance)")
axs[2].set_ylabel("Cell Type Interaction Score")
axs[2].set_title("Cell Type Interaction Score Distribution")

plt.tight_layout()
plt.savefig(f'{script_dir}/SFig3_panel_B1.png')

fig, axes = plt.subplots(3, 3, figsize=(9, 9), dpi=300)
axes = axes.flatten()

for i, neighborhood_size in enumerate(neighborhood_sizes):
    param = ss.util.generate_random_parameters(n_group=2, n_state=9, seed=0)
    sim = ss.util.sim_from_params(
        param,
        shape=(100, 100),
        num_iteration=4,
        n_iter=6,
        custom_neighbor=ss.spatial.generate_offsets(neighborhood_size, 'manhattan'),
        seed=0
    )

    ax = axes[i]
    sns.scatterplot(
        x='row', y='col', data=sim.meta,
        hue='state', palette='tab20', legend=False, s=12, ax=ax
    )
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_title(f'Neighborhood size = {neighborhood_size}')
    ax.set_xticks([])
    ax.set_yticks([])

plt.tight_layout()
plt.savefig(f'{script_dir}/SFig3_panel_B2.png')


########### Panel C:
def ss_statistics_s(
        shape = (100, 100), 
        n_group = 3, 
        n_state = 9, 
        num_iteration = 4, 
        n_iter = 6, 
        sigma = 0,
        seed = 0,  
    ):
    param_1 = ss.util.generate_random_parameters(n_group=n_group, n_state=n_state, seed=seed)

    n_group = param_1['n_group']
    n_state = param_1['n_state']
    
    niche_theta = np.zeros((n_group, n_group))
    niche_theta[np.triu_indices(n_group, 1)] = param_1['niche_theta']
    niche_theta = niche_theta + niche_theta.T - np.diag(niche_theta.diagonal())
    np.fill_diagonal(niche_theta, 1)

    theta_list = []
    for i in range(n_group):
        theta_tmp = np.zeros((n_state, n_state))
        theta_tmp[np.triu_indices(n_state, 1)] = param_1['theta_list'][i]
        theta_tmp = theta_tmp + theta_tmp.T - np.diag(theta_tmp.diagonal())
        np.fill_diagonal(theta_tmp, 1)
        theta_list.append(theta_tmp)
    
    density_replicates = np.array(param_1['density_replicates'])
    density_replicates[density_replicates < 0] = 0
    phi_replicates = param_1['phi_replicates']

    Sim = ss.SimSpace(
        shape = shape,
        num_states = n_state,
        num_iterations= num_iteration,
        theta=theta_list,
        phi=phi_replicates,
        neighborhood=ss.spatial.generate_offsets(3, 'manhattan'), 
        random_seed=seed,
        )
    Sim.initialize()   # Initialize the grid
    Sim.create_niche(num_niches=n_group, n_iter=n_iter, theta_niche=niche_theta)
    Sim.gibbs_sampler()    # Gibbs sampling
    Sim.density_sampler(density_replicates)  # Cell density of each niche
    Sim.perturbation(step = sigma)     # Perturbation

    mi_sim = ss.spatial.integrate_morans_I(Sim.meta['state'], Sim.meta[['row', 'col']], Sim.meta['state'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(Sim.meta['state'], Sim.meta[['row', 'col']], Sim.meta['state'].unique().tolist())
    ci_sim = ss.spatial.calculate_interaction_score(Sim.meta['state'], Sim.meta[['row', 'col']], Sim.meta['state'].unique().tolist())
    return [mi_sim, gc_sim, ci_sim]


sigma_list = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
sim_neighbor_mi = [[] for _ in range(6)]
sim_neighbor_gc = [[] for _ in range(6)]
sim_neighbor_ci = [[] for _ in range(6)]

for i, sigma in enumerate(sigma_list):
    for j in range(8):
        ## SimSpace
        sim_stats = ss_statistics_s(
            shape=(100, 100), n_group=2, n_state=9, num_iteration=4, n_iter=6, 
            sigma=sigma, seed=j)
        sim_neighbor_mi[i].append(sim_stats[0])
        sim_neighbor_gc[i].append(sim_stats[1])
        sim_neighbor_ci[i].append(sim_stats[2].values.reshape(-1).tolist())
sim_neighbor_mi = [np.concatenate(item) for item in sim_neighbor_mi]
sim_neighbor_gc = [np.concatenate(item) for item in sim_neighbor_gc]
sim_neighbor_ci = [np.concatenate(item) for item in sim_neighbor_ci]


plt_data_mi = [
    Xenium_mi_sim, merfish_mi_sim, sim_neighbor_mi[0], sim_neighbor_mi[1], 
    sim_neighbor_mi[2], sim_neighbor_mi[3], sim_neighbor_mi[4], sim_neighbor_mi[5], 
]

plt_data_gc = [
    Xenium_gc_sim, merfish_gc_sim, sim_neighbor_gc[0], sim_neighbor_gc[1], 
    sim_neighbor_gc[2], sim_neighbor_gc[3], sim_neighbor_gc[4], sim_neighbor_gc[5], 
]

plt_data_ci = [
    Xenium_ci_sim, merfish_ci_sim, sim_neighbor_ci[0], sim_neighbor_ci[1], 
    sim_neighbor_ci[2], sim_neighbor_ci[3], sim_neighbor_ci[4], sim_neighbor_ci[5], 
]

# Prepare data for violin plots
labels = [
    "Xenium", "MERFISH", "0", 
    "0.1", "0.2", "0.3", "0.4", "0.5", 
]

fig, axs = plt.subplots(1, 3, figsize=(10, 3.5), dpi=300)
palette = sns.color_palette("muted", 11)
palette = list(palette)
first = np.array(palette[0])
# blend the first color toward a light blue and set as the second color
palette[1] = tuple(np.clip(first * 0.4 + np.array([0.7, 0.85, 1.0]) * 0.6, 0, 1))
palette[2:] =  sns.color_palette("viridis", 6)

# Panel 1: Moran's I
sns.violinplot(data=plt_data_mi, ax=axs[0], palette=palette)
axs[0].set_xticklabels(labels, rotation=30, ha='center')
axs[0].set_xlabel("Perturbation Level (Sigma)")
axs[0].set_ylabel("Moran's I")
axs[0].set_title("Moran's I Distribution")

# Panel 2: Geary's C
sns.violinplot(data=plt_data_gc, ax=axs[1], palette=palette)
axs[1].set_xticklabels(labels, rotation=30, ha='center')
axs[1].set_xlabel("Perturbation Level (Sigma)")
axs[1].set_ylabel("Geary's C")
axs[1].set_title("Geary's C Distribution")

# Panel 3: Cell Type Interaction
sns.violinplot(data=plt_data_ci, ax=axs[2], palette=palette)
axs[2].set_xticklabels(labels, rotation=30, ha='center')
axs[2].set_xlabel("Perturbation Level (Sigma)")
axs[2].set_ylabel("Cell Type Interaction Score")
axs[2].set_title("Cell Type Interaction Score Distribution")

plt.tight_layout()
plt.savefig(f'{script_dir}/SFig3_panel_C1.png')

fig, axes = plt.subplots(3, 3, figsize=(9, 9), dpi=300)
axes = axes.flatten()

for j, sigma in enumerate(sigma_list):
    param_1 = ss.util.generate_random_parameters(n_group=2, n_state=9, seed=0)
    n_group = param_1['n_group']
    n_state = param_1['n_state'] 
    niche_theta = np.zeros((n_group, n_group))
    niche_theta[np.triu_indices(n_group, 1)] = param_1['niche_theta']
    niche_theta = niche_theta + niche_theta.T - np.diag(niche_theta.diagonal())
    np.fill_diagonal(niche_theta, 1)
    theta_list = []
    for i in range(n_group):
        theta_tmp = np.zeros((n_state, n_state))
        theta_tmp[np.triu_indices(n_state, 1)] = param_1['theta_list'][i]
        theta_tmp = theta_tmp + theta_tmp.T - np.diag(theta_tmp.diagonal())
        np.fill_diagonal(theta_tmp, 1)
        theta_list.append(theta_tmp)
    
    density_replicates = np.array(param_1['density_replicates'])
    density_replicates[density_replicates < 0] = 0
    phi_replicates = param_1['phi_replicates']

    Sim = ss.SimSpace(
        shape = (100, 100),
        num_states = n_state,
        num_iterations= 6,
        theta=theta_list,
        phi=phi_replicates,
        neighborhood=ss.spatial.generate_offsets(3, 'manhattan'), 
        random_seed=0,
        )
    Sim.initialize()   # Initialize the grid
    Sim.create_niche(num_niches=n_group, n_iter=6, theta_niche=niche_theta)
    Sim.gibbs_sampler()    # Gibbs sampling
    Sim.density_sampler(density_replicates)  # Cell density of each niche
    Sim.perturbation(step = sigma)     # Perturbation

    ax = axes[j]
    sns.scatterplot(
        x='row', y='col', data=Sim.meta,
        hue='state', palette='tab20', legend=False, s=12, ax=ax
    )
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_title(f'Perturbation Sigma = {sigma}')
    ax.set_xticks([])
    ax.set_yticks([])

plt.tight_layout()
plt.savefig(f'{script_dir}/SFig3_panel_C2.png')