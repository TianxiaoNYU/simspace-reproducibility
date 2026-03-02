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
import os

import simspace as ss

script_dir_global = os.path.dirname(os.path.abspath(__file__))
script_dir = os.path.join(script_dir_global, '../../../')

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
SRTsim_stats_list = []
SRTsim_perturbed1_list = []
SRTsim_perturbed2_list = []

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
    ## SRTsim_perturbed_1
    tile_path = f'{script_dir}/Panel_D_data/SRTsim/supp_foldchange/perturbed_1/ref_free_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    mi_sim = ss.spatial.integrate_morans_I(cells['group'], cells[['x', 'y']], cells['group'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(cells['group'], cells[['x', 'y']], cells['group'].unique().tolist())
    cells['x'] = (cells['x'] - cells['x'].min()) / (cells['x'].max() - cells['x'].min()) * 100
    cells['y'] = (cells['y'] - cells['y'].min()) / (cells['y'].max() - cells['y'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['group'], cells[['x', 'y']], cells['group'].unique().tolist())
    L_sim = ripley_L_simple(cells[['x', 'y']].values, 100)
    SRTsim_perturbed1_list.append([mi_sim, gc_sim, ci_sim, L_sim])
    ## SRTsim_perturbed_2
    tile_path = f'{script_dir}/Panel_D_data/SRTsim/supp_foldchange/perturbed_2/ref_free_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    mi_sim = ss.spatial.integrate_morans_I(cells['group'], cells[['x', 'y']], cells['group'].unique().tolist())
    gc_sim = ss.spatial.integrate_gearys_C(cells['group'], cells[['x', 'y']], cells['group'].unique().tolist())
    cells['x'] = (cells['x'] - cells['x'].min()) / (cells['x'].max() - cells['x'].min()) * 100
    cells['y'] = (cells['y'] - cells['y'].min()) / (cells['y'].max() - cells['y'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['group'], cells[['x', 'y']], cells['group'].unique().tolist())
    L_sim = ripley_L_simple(cells[['x', 'y']].values, 100)
    SRTsim_perturbed2_list.append([mi_sim, gc_sim, ci_sim, L_sim])
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
srt_mi = np.sort(np.concatenate([item[0] for item in SRTsim_stats_list]))
srt_p1_mi = np.sort(np.concatenate([item[0] for item in SRTsim_perturbed1_list]))
srt_p2_mi = np.sort(np.concatenate([item[0] for item in SRTsim_perturbed2_list]))
merfish_mi = np.sort(np.concatenate([item[0] for item in merfish_stats_list]))
## CDF
xenium_cdf = np.arange(1, len(xenium_mi) + 1) / len(xenium_mi)
sim_cdf = np.arange(1, len(sim_mi) + 1) / len(sim_mi)
srt_cdf = np.arange(1, len(srt_mi) + 1) / len(srt_mi)
srt_p1_cdf = np.arange(1, len(srt_p1_mi) + 1) / len(srt_p1_mi)
srt_p2_cdf = np.arange(1, len(srt_p2_mi) + 1) / len(srt_p2_mi)
merfish_cdf = np.arange(1, len(merfish_mi) + 1) / len(merfish_mi)

## GC
sim_gc = np.sort(np.concatenate([item[1] for item in sim_stats_list]))
xenium_gc = np.sort(np.concatenate([item[1] for item in Xenium_stats_list]))
srt_gc = np.sort(np.concatenate([item[1] for item in SRTsim_stats_list]))
srt_p1_gc = np.sort(np.concatenate([item[1] for item in SRTsim_perturbed1_list]))
srt_p2_gc = np.sort(np.concatenate([item[1] for item in SRTsim_perturbed2_list]))
merfish_gc = np.sort(np.concatenate([item[1] for item in merfish_stats_list]))
## CDF
gc_xenium_cdf = np.arange(1, len(xenium_gc) + 1) / len(xenium_gc)
gc_sim_cdf = np.arange(1, len(sim_gc) + 1) / len(sim_gc)
gc_srt_cdf = np.arange(1, len(srt_gc) + 1) / len(srt_gc)
gc_srt_p1_cdf = np.arange(1, len(srt_p1_gc) + 1) / len(srt_p1_gc)
gc_srt_p2_cdf = np.arange(1, len(srt_p2_gc) + 1) / len(srt_p2_gc)
gc_merfish_cdf = np.arange(1, len(merfish_gc) + 1) / len(merfish_gc)

## L function
sim_L_arr = np.array([item[3] for item in sim_stats_list])
xenium_L_arr = np.array([item[3] for item in Xenium_stats_list])
srt_L_arr = np.array([item[3] for item in SRTsim_stats_list])
srt_p1_L_arr = np.array([item[3] for item in SRTsim_perturbed1_list])
srt_p2_L_arr = np.array([item[3] for item in SRTsim_perturbed2_list])
merfish_L_arr = np.array([item[3] for item in merfish_stats_list])

sim_L_arr = np.mean(sim_L_arr, axis=0)
xenium_L_arr = np.mean(xenium_L_arr, axis=0)
srt_L_arr = np.mean(srt_L_arr, axis=0)
srt_p1_L_arr = np.mean(srt_p1_L_arr, axis=0)
srt_p2_L_arr = np.mean(srt_p2_L_arr, axis=0)
merfish_L_arr = np.mean(merfish_L_arr, axis=0)


fig, axs = plt.subplots(1, 3, figsize=(10, 3.5), dpi=300)
palette = sns.color_palette("muted", 6)
palette = list(palette)
first = np.array(palette[0])
# blend the first color toward a light blue and set as the second color
palette[1] = tuple(np.clip(first * 0.4 + np.array([0.7, 0.85, 1.0]) * 0.6, 0, 1))

# Prepare data for violin plots
mi_data = [
    xenium_mi, merfish_mi, sim_mi, srt_mi, srt_p1_mi, srt_p2_mi
]
mi_labels = ['Xenium', 'MERFISH', 'SimSpace', 'SRTsim', 'Replicate 1', 'Replicate 2']

gc_data = [
    xenium_gc, merfish_gc, sim_gc, srt_gc, srt_p1_gc, srt_p2_gc
]
gc_labels = mi_labels

L_data = [
    xenium_L_arr, merfish_L_arr, sim_L_arr, srt_L_arr, srt_p1_L_arr, srt_p2_L_arr
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
plt.savefig(f'{script_dir}/Panel_D_data/SRTsim/supp_foldchange/Fig2_panel_D.png', bbox_inches='tight', dpi=300)
plt.close()  