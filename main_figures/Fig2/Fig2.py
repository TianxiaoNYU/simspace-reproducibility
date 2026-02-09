import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from math import pi
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc
import os

import simspace as ss


########## Panel B:
shape = (100, 100)
script_dir = os.path.dirname(os.path.abspath(__file__))
param_file = os.path.join(script_dir, 'fitted_params.json')
param_1 = ss.util.generate_random_parameters(n_group=3, n_state=8, seed=42)

sim1 = ss.util.sim_from_params(
    param_1,
    shape=shape,
    num_iteration=4,
    n_iter=6,
    custom_neighbor=ss.spatial.generate_offsets(3, 'manhattan'),
    seed=0
)

sim1.plot(figsize=(5, 5), dpi=300, size=14)
plt.savefig(f'{script_dir}/Fig2_panel_B1.png', bbox_inches='tight', dpi=300)
plt.close()  # Close the figure to prevent it from showing

shape = (100, 100)
param_1 = ss.util.generate_random_parameters(n_group=2, n_state=9, seed=42)

sim1 = ss.util.sim_from_params(
    param_1,
    shape=shape,
    num_iteration=4,
    n_iter=6,
    custom_neighbor=ss.spatial.generate_offsets(3, 'manhattan'),
    seed=0
)

sim1.plot(figsize=(5, 5), dpi=300, size=14)
plt.savefig(f'{script_dir}/Fig2_panel_B2.png', bbox_inches='tight', dpi=300)
plt.close()  # Close the figure to prevent it from showing


shape = (100, 100)
param_1 = ss.util.generate_random_parameters(n_group=3, n_state=9, seed=42)

n_group = param_1['n_group']
diag_length = len(param_1['theta_list'][0])
n_state = int((1 + np.sqrt(1 + 8 * diag_length)) / 2)
if n_state * (n_state - 1) != 2 * diag_length:
    raise ValueError("Invalid theta matrix size.")

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
    num_iterations= 4,
    theta=theta_list,
    phi=phi_replicates,
    neighborhood=ss.spatial.generate_offsets(3, 'manhattan'), 
    random_seed=42,
    )
Sim.initialize()   # Initialize the grid
x, y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
distance_from_center = np.sqrt((x - shape[1] // 2)**2 + (y - shape[0] // 2)**2)
Sim.grid[distance_from_center > 50] = -1
Sim.create_niche(num_niches=n_group, n_iter=6, theta_niche=niche_theta)
Sim.gibbs_sampler()    # Gibbs sampling
Sim.density_sampler(density_replicates)  # Cell density of each niche
Sim.perturbation(step = 0.2)
Sim.plot(figsize=(5, 5), dpi=300, size=14)
plt.savefig(f'{script_dir}/Fig2_panel_B3.png', bbox_inches='tight', dpi=300)
plt.close()  # Close the figure to prevent it from showing


shape = (100, 100)
param_1 = ss.util.generate_random_parameters(n_group=3, n_state=9, seed=1111)

n_group = param_1['n_group']
diag_length = len(param_1['theta_list'][0])
n_state = int((1 + np.sqrt(1 + 8 * diag_length)) / 2)
if n_state * (n_state - 1) != 2 * diag_length:
    raise ValueError("Invalid theta matrix size.")

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
    num_iterations= 4,
    theta=theta_list,
    phi=phi_replicates,
    neighborhood=ss.spatial.generate_offsets(3, 'manhattan'), 
    random_seed=42,
    )
Sim.initialize()   # Initialize the grid

x, y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
distance_from_center = np.sqrt((x - shape[1] // 2)**2 + (y - shape[0] // 2)**2)

# Introduce anisotropy by scaling x and y differently
anisotropy_factor_x = 1.2  # Stretch in x-direction
anisotropy_factor_y = 0.8  # Compress in y-direction
anisotropic_distance = np.sqrt(((x - shape[1] // 2) * anisotropy_factor_x)**2 + 
                                ((y - shape[0] // 2) * anisotropy_factor_y)**2)

random_noise = np.random.normal(0, 5, size=anisotropic_distance.shape)
irregular_shape = anisotropic_distance + random_noise
Sim.grid[irregular_shape > 50] = -1

Sim.create_niche(num_niches=n_group, n_iter=6, theta_niche=niche_theta)
Sim.gibbs_sampler()    # Gibbs sampling
Sim.density_sampler(density_replicates)  # Cell density of each niche
Sim.perturbation(step = 0.2)
Sim.plot(figsize=(5, 5), dpi=300, size=14)
plt.savefig(f'{script_dir}/Fig2_panel_B4.png', bbox_inches='tight', dpi=300)
plt.close()  # Close the figure to prevent it from showing



########## Panel C:
shape = (100, 100)
param_1 =ss.util.generate_random_parameters(n_group=3, n_state=9, seed=42)

n_group = param_1['n_group']
diag_length = len(param_1['theta_list'][0])
n_state = int((1 + np.sqrt(1 + 8 * diag_length)) / 2)
if n_state * (n_state - 1) != 2 * diag_length:
    raise ValueError("Invalid theta matrix size.")

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
    num_iterations= 4,
    theta=theta_list,
    phi=phi_replicates,
    neighborhood=ss.spatial.generate_offsets(3, 'manhattan'), 
    random_seed=42,
    )
Sim.initialize()   # Initialize the grid
x, y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
distance_from_center = np.sqrt((x - shape[1] // 2)**2 + (y - shape[0] // 2)**2)
Sim.grid[distance_from_center > 50] = -1
Sim.create_niche(num_niches=n_group, n_iter=6, theta_niche=niche_theta)
Sim.gibbs_sampler()    # Gibbs sampling
Sim.density_sampler(density_replicates)  # Cell density of each niche
Sim.perturbation(step = 0.2)
Sim.create_omics(n_genes=1000, bg_ratio=0.6, bg_param = (1, 0.5), marker_param = (3.5, 1.6), spatial=False)
ss.plot.plot_gene(
    Sim.meta, Sim.omics['Gene_2'], size=14, 
    figsize=(5,5), dpi=300, title='Marker Gene Expression',)
plt.savefig(f'{script_dir}/Fig2_panel_C.png', bbox_inches='tight', dpi=300)
plt.close()  # Close the figure to prevent it from showing

########## Panel D:
def ripley_K_simple(points, grid_size, radii):
    """
    Uncorrected Ripley's K for a rectangular window using cKDTree.
    
    points: (N,2) array of [x,y]
    window: (xmin, xmax, ymin, ymax)
    radii : 1D array of r values (ascending)
    
    Returns: K(r) array same shape as radii
    """
    points = np.asarray(points, float)
    radii  = np.asarray(radii, float)
    n = len(points)
    if n < 2:
        return np.zeros_like(radii)

    A = grid_size ** 2

    # build once and get all unordered pairs up to max r
    rmax = float(radii.max())
    tree = cKDTree(points)
    pairs = list(tree.query_pairs(rmax))  # list of (i,j) with i<j
    if not pairs:
        return np.zeros_like(radii)

    # pairwise distances for those pairs
    p = np.array(pairs, dtype=int)
    d = np.hypot(points[p[:,0],0] - points[p[:,1],0],
                 points[p[:,0],1] - points[p[:,1],1])

    # sort once → cumulative counts
    order = np.argsort(d)
    d_sorted = d[order]
    cum_pairs = np.arange(1, len(d_sorted)+1)  # 1..m

    # for each r, number of pairs with d <= r
    K = np.zeros_like(radii, dtype=float)
    factor = A / (n * (n - 1)) * 2.0  # ×2 because sum_{i≠j} counts both directions
    for k, r in enumerate(radii):
        idx = np.searchsorted(d_sorted, r, side='right') - 1
        if idx >= 0:
            K[k] = factor * cum_pairs[idx]
    return K

def ripley_L_simple(points, grid_size):
    radii = np.linspace(0, grid_size/4, 25)
    K = ripley_K_simple(points, grid_size, radii)
    L = np.sqrt(K / pi)  # often plot L(r) - r
    return L - radii

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
    lfunc_sim = ripley_L_simple(sim1.meta[['col', 'row']].values, shape[0])

    return [mi_sim, gc_sim, ci_sim, lfunc_sim]


sim_stats_list = []
Xenium_stats_list = []
merfish_stats_list = []
sccube_stats_list = []
SRTsim_stats_list = []
scMultiSim_stats_list = []
random_control_list = []
for i in range(8):
    ## SimSpace
    sim_stats = ss_statistics(
        shape=(100, 100), n_group=2, n_state=9, num_iteration=4, n_iter=6, 
        custom_neighbor=ss.spatial.generate_offsets(3, 'manhattan'), seed=i)
    sim_stats_list.append(sim_stats)
    ## Xenium
    tile_path = f'{script_dir}/Panel_D_data/Xenium/tile_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    mi_sim = ss.spatial.integrate_morans_I(cells['Cluster'], cells[['x_centroid', 'y_centroid']], cells['Cluster'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(cells['Cluster'], cells[['x_centroid', 'y_centroid']], cells['Cluster'].unique().tolist())
    cells['x_centroid'] = (cells['x_centroid'] - cells['x_centroid'].min()) / (cells['x_centroid'].max() - cells['x_centroid'].min()) * 100
    cells['y_centroid'] = (cells['y_centroid'] - cells['y_centroid'].min()) / (cells['y_centroid'].max() - cells['y_centroid'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['Cluster'], cells[['x_centroid', 'y_centroid']], cells['Cluster'].unique().tolist())
    L_sim = ripley_L_simple(cells[['x_centroid', 'y_centroid']].values, 100)
    Xenium_stats_list.append([mi_sim, gc_sim, ci_sim, L_sim])
    ## scCube
    tile_path = f'{script_dir}/Panel_D_data/sccube/sc_meta_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    mi_sim = ss.spatial.integrate_morans_I(cells['Cell_type'], cells[['point_x', 'point_y']], cells['Cell_type'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(cells['Cell_type'], cells[['point_x', 'point_y']], cells['Cell_type'].unique().tolist())
    cells['point_x'] = (cells['point_x'] - cells['point_x'].min()) / (cells['point_x'].max() - cells['point_x'].min()) * 100
    cells['point_y'] = (cells['point_y'] - cells['point_y'].min()) / (cells['point_y'].max() - cells['point_y'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['Cell_type'], cells[['point_x', 'point_y']], cells['Cell_type'].unique().tolist())
    L_sim = ripley_L_simple(cells[['point_x', 'point_y']].values, 100)
    sccube_stats_list.append([mi_sim, gc_sim, ci_sim, L_sim])
    ## SRTsim   
    tile_path = f'{script_dir}/Panel_D_data/SRTsim/ref_free_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    mi_sim = ss.spatial.integrate_morans_I(cells['group'], cells[['x', 'y']], cells['group'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(cells['group'], cells[['x', 'y']], cells['group'].unique().tolist())
    cells['x'] = (cells['x'] - cells['x'].min()) / (cells['x'].max() - cells['x'].min()) * 100
    cells['y'] = (cells['y'] - cells['y'].min()) / (cells['y'].max() - cells['y'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['group'], cells[['x', 'y']], cells['group'].unique().tolist())
    L_sim = ripley_L_simple(cells[['x', 'y']].values, 100)
    SRTsim_stats_list.append([mi_sim, gc_sim, ci_sim, L_sim])
    ## scMultiSim
    tile_path = f'{script_dir}/Panel_D_data/scMultiSim/spatial_sim_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    mi_sim = ss.spatial.integrate_morans_I(cells['cell.type'], cells[['x', 'y']], cells['cell.type'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(cells['cell.type'], cells[['x', 'y']], cells['cell.type'].unique().tolist())
    cells['x'] = (cells['x'] - cells['x'].min()) / (cells['x'].max() - cells['x'].min()) * 100
    cells['y'] = (cells['y'] - cells['y'].min()) / (cells['y'].max() - cells['y'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['cell.type'], cells[['x', 'y']], cells['cell.type'].unique().tolist())
    L_sim = ripley_L_simple(cells[['x', 'y']].values, 100)
    scMultiSim_stats_list.append([mi_sim, gc_sim, ci_sim, L_sim])
    ## MERFISH
    tile_path = f'{script_dir}/Panel_D_data/MERFISH/merfish_data_sample_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    mi_sim = ss.spatial.integrate_morans_I(cells['Cell_class'], cells[['Centroid_X', 'Centroid_Y']], cells['Cell_class'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(cells['Cell_class'], cells[['Centroid_X', 'Centroid_Y']], cells['Cell_class'].unique().tolist())
    cells['Centroid_X'] = (cells['Centroid_X'] - cells['Centroid_X'].min()) / (cells['Centroid_X'].max() - cells['Centroid_X'].min()) * 100
    cells['Centroid_Y'] = (cells['Centroid_Y'] - cells['Centroid_Y'].min()) / (cells['Centroid_Y'].max() - cells['Centroid_Y'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['Cell_class'], cells[['Centroid_X', 'Centroid_Y']], cells['Cell_class'].unique().tolist())
    L_sim = ripley_L_simple(cells[['Centroid_X', 'Centroid_Y']].values, 100)
    merfish_stats_list.append([mi_sim, gc_sim, ci_sim, L_sim])

## MI
sim_mi = np.sort(np.concatenate([item[0] for item in sim_stats_list]))
xenium_mi = np.sort(np.concatenate([item[0] for item in Xenium_stats_list]))
sccube_mi = np.sort(np.concatenate([item[0] for item in sccube_stats_list]))
srt_mi = np.sort(np.concatenate([item[0] for item in SRTsim_stats_list]))
scMultiSim_mi = np.sort(np.concatenate([item[0] for item in scMultiSim_stats_list]))
merfish_mi = np.sort(np.concatenate([item[0] for item in merfish_stats_list]))
## CDF
xenium_cdf = np.arange(1, len(xenium_mi) + 1) / len(xenium_mi)
sim_cdf = np.arange(1, len(sim_mi) + 1) / len(sim_mi)
sccube_cdf = np.arange(1, len(sccube_mi) + 1) / len(sccube_mi)
srt_cdf = np.arange(1, len(srt_mi) + 1) / len(srt_mi)
scMultiSim_cdf = np.arange(1, len(scMultiSim_mi) + 1) / len(scMultiSim_mi)
merfish_cdf = np.arange(1, len(merfish_mi) + 1) / len(merfish_mi)

## GC
sim_gc = np.sort(np.concatenate([item[1] for item in sim_stats_list]))
xenium_gc = np.sort(np.concatenate([item[1] for item in Xenium_stats_list]))
sccube_gc = np.sort(np.concatenate([item[1] for item in sccube_stats_list]))
srt_gc = np.sort(np.concatenate([item[1] for item in SRTsim_stats_list]))
scMultiSim_gc = np.sort(np.concatenate([item[1] for item in scMultiSim_stats_list]))
merfish_gc = np.sort(np.concatenate([item[1] for item in merfish_stats_list]))
## CDF
gc_xenium_cdf = np.arange(1, len(xenium_gc) + 1) / len(xenium_gc)
gc_sim_cdf = np.arange(1, len(sim_gc) + 1) / len(sim_gc)
gc_sccube_cdf = np.arange(1, len(sccube_gc) + 1) / len(sccube_gc)
gc_srt_cdf = np.arange(1, len(srt_gc) + 1) / len(srt_gc)
gc_scMultiSim_cdf = np.arange(1, len(scMultiSim_gc) + 1) / len(scMultiSim_gc)
gc_merfish_cdf = np.arange(1, len(merfish_gc) + 1) / len(merfish_gc)

## L function
sim_L_arr = np.array([item[3] for item in sim_stats_list])
xenium_L_arr = np.array([item[3] for item in Xenium_stats_list])
sccube_L_arr = np.array([item[3] for item in sccube_stats_list])
srt_L_arr = np.array([item[3] for item in SRTsim_stats_list])
scMultiSim_L_arr = np.array([item[3] for item in scMultiSim_stats_list])
merfish_L_arr = np.array([item[3] for item in merfish_stats_list])
sim_L_arr = np.mean(sim_L_arr, axis=0)
xenium_L_arr = np.mean(xenium_L_arr, axis=0)
sccube_L_arr = np.mean(sccube_L_arr, axis=0)
srt_L_arr = np.mean(srt_L_arr, axis=0)
scMultiSim_L_arr = np.mean(scMultiSim_L_arr, axis=0)
merfish_L_arr = np.mean(merfish_L_arr, axis=0)


fig, axs = plt.subplots(1, 3, figsize=(10, 3.5), dpi=300)
palette = sns.color_palette("muted", 6)
palette = list(palette)
first = np.array(palette[0])
# blend the first color toward a light blue and set as the second color
palette[1] = tuple(np.clip(first * 0.4 + np.array([0.7, 0.85, 1.0]) * 0.6, 0, 1))

# Prepare data for violin plots
mi_data = [
    xenium_mi, merfish_mi, sim_mi, sccube_mi, srt_mi, scMultiSim_mi, 
]
mi_labels = ['Xenium', 'MERFISH', 'SimSpace', 'scCube', 'SRTsim', 'scMultiSim']

gc_data = [
    xenium_gc, merfish_gc, sim_gc, sccube_gc, srt_gc, scMultiSim_gc
]
gc_labels = mi_labels

L_data = [
    xenium_L_arr, merfish_L_arr, sim_L_arr, sccube_L_arr, srt_L_arr, scMultiSim_L_arr
]
L_labels = mi_labels

# Panel 1: Moran's I
sns.violinplot(data=mi_data, ax=axs[0], palette=palette)
axs[0].set_xticklabels(mi_labels, rotation=30, ha='center')
axs[0].set_ylabel("Moran's I")
axs[0].set_title("Moran's I Distribution")

# Panel 2: Geary's C
sns.violinplot(data=gc_data, ax=axs[1], palette=palette)
axs[1].set_xticklabels(gc_labels, rotation=30, ha='center')
axs[1].set_ylabel("Geary's C")
axs[1].set_title("Geary's C Distribution")

# Panel 3: Ripley's L (show as distribution at r=12, or mean across r)
L_r_idx = len(sim_L_arr) // 2  # middle radius
L_violin = [arr[L_r_idx] if arr.ndim == 1 else np.mean(arr) for arr in L_data]
sns.violinplot(data=[arr for arr in L_data], ax=axs[2], palette=palette)
axs[2].set_xticklabels(L_labels, rotation=30, ha='center')
axs[2].set_ylabel("Ripley's L (L(r)-r)")
axs[2].set_title("Ripley's L Function Distribution")

plt.tight_layout()
plt.savefig(f'{script_dir}/Fig2_panel_D.png', bbox_inches='tight', dpi=300)
plt.close()  

########## Panel E:
shape = (100, 100, 40)
num_iteration = 5
custom_neighbor = ss.spatial.generate_offsets3D(4, 'manhattan')

# parameters = load_parameters_from_json('/Users/zhaotianxiao/Library/CloudStorage/Dropbox/FenyoLab/Project/Spatialsim/output/best_solution_breast_23_0331_sparse.json')
parameters = ss.util.generate_random_parameters(n_group=3, n_state=8, seed=42)
n_state = parameters['n_state']
n_group = parameters['n_group']

niche_theta = np.zeros((n_group, n_group))
niche_theta[np.triu_indices(n_group, 1)] = parameters['niche_theta']
niche_theta = niche_theta + niche_theta.T - np.diag(niche_theta.diagonal())
np.fill_diagonal(niche_theta, 1)

theta_list = []
for i in range(n_group):
    theta_tmp = np.zeros((n_state, n_state))
    theta_tmp[np.triu_indices(n_state, 1)] = parameters['theta_list'][i]
    theta_tmp = theta_tmp + theta_tmp.T - np.diag(theta_tmp.diagonal())
    np.fill_diagonal(theta_tmp, 1)
    theta_list.append(theta_tmp)

density_replicates = np.array(parameters['density_replicates'])
density_replicates[density_replicates < 0] = 0
phi_replicates = parameters['phi_replicates']

Sim1 = ss.SimSpace(
    shape = shape,
    num_states = n_state,
    num_iterations= num_iteration,
    theta=theta_list,
    phi=phi_replicates,
    neighborhood=custom_neighbor, 
    random_seed=0)
Sim1.initialize3D()   # Initialize the grid
Sim1.create_niche3D(num_niches=n_group, n_iter=4, theta_niche=niche_theta)
Sim1.gibbs_sampler3D()    # Gibbs sampling
Sim1.density_sampler(density_replicates)  # Cell density of each niche
Sim1.perturbation3D(step = 0.2)     # Perturbation

from matplotlib.lines import Line2D
# Create a DataFrame for easier pslotting
df = pd.DataFrame(Sim1.meta, columns=['x', 'y', 'z', 'state'])
cmap = sns.color_palette(cc.glasbey, n_colors=8)
# cmap = sns.color_palette('tab20', n_colors=8)

# Plot the 3D scatter plot
fig = plt.figure(figsize=(6, 4), dpi=300)
ax = fig.add_subplot(111, projection='3d')
df['state'] = df['state'].astype('category')
scatter = ax.scatter(df['x'], df['y'], df['z'], c=df['state'].cat.codes, cmap=plt.cm.colors.ListedColormap(cmap), alpha=0.5, linewidths=0, s=3)
# Add legend for states
state_labels = [f'Type {cat + 1}' for cat in df['state'].cat.categories]
legend_elements = [Line2D([0], [0], marker='o', color='w', label=label,
                          markerfacecolor=cmap[i], markersize=6)
                   for i, label in enumerate(state_labels)]
ax.legend(handles=legend_elements, title='State', bbox_to_anchor=(1.25, 1), loc='upper left')
ax.set_aspect('equal')
ax.set_zticks([0, 20, 40])

# Set labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_title('3D State Distribution')

plt.savefig(f'{script_dir}/Fig2_panel_E.png', bbox_inches='tight', dpi=300)
plt.close()  

########## Panel F:
Sim1.create_omics(n_genes=100, bg_ratio=0, bg_param = (1, 0.5), marker_param = (5, 2), spatial=False)
df = pd.DataFrame(Sim1.meta, columns=['x', 'y', 'z', 'state'])
df['gene'] = Sim1.omics['Gene_48']
cmap = sns.color_palette('rocket_r', as_cmap=True)
df['alpha'] = df['gene']/df['gene'].max()
df['alpha'][df['alpha'] < 0.4] = 0.2
df['alpha'][df['alpha'] > 0.4] = 0.7

# Plot the 3D scatter plot
fig = plt.figure(figsize=(6, 3), dpi=300)
ax = fig.add_subplot(111, projection='3d')
df['state'] = df['state'].astype('category')
scatter = ax.scatter(df['x'], df['y'], df['z'], c=df['gene'], alpha=df['alpha'], linewidths=0, s=3, cmap=cmap)
# Add legend
cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.14)
# cbar.set_label(feature.name)
ax.set_aspect('equal')
ax.set_zticks([0, 20, 40])

# Set labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_title('Marker Gene Expression')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig2_panel_F.png', bbox_inches='tight', dpi=300)
plt.close()
