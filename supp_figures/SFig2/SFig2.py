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
shape = (100, 100)
param_1 = ss.util.generate_random_parameters(n_group=2, n_state=7, seed=1)

n_group = param_1['n_group']
diag_length = len(param_1['theta_list'][0])
n_state = int((1 + np.sqrt(1 + 8 * diag_length)) / 2)
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
    theta=theta_list, # type: ignore
    phi=phi_replicates,
    neighborhood=ss.spatial.generate_offsets(3, 'manhattan'), 
    random_seed=42,
    )
Sim.initialize()   # Initialize the grid
## Create niches and run simulation
Sim.create_niche(num_niches=n_group, n_iter=6, theta_niche=niche_theta)
## Create a trapezoid boundary in the grid, which takes half of the area
x, y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
boundary_mask = (y < ((shape[0] / shape[1]) * x * 2 - 100))  # Trapezoid boundary
Sim.grid[boundary_mask] = -1

Sim.gibbs_sampler()    # Gibbs sampling
Sim.density_sampler(density_replicates) # Cell density of each niche # type: ignore 
Sim.perturbation(step = 0.2)
Sim.plot(figsize=(5, 5), dpi=300, size=14)
plt.savefig(os.path.join(script_dir, 'SFig2_panel_A1.png'), bbox_inches='tight', dpi=300)
Sim.plot_niche(figsize=(5, 5), dpi=300)
plt.savefig(os.path.join(script_dir, 'SFig2_panel_A2.png'), bbox_inches='tight', dpi=300)

shape = (100, 100)
param_1 = ss.util.generate_random_parameters(n_group=4, n_state=10, seed=42)

n_group = param_1['n_group']
diag_length = len(param_1['theta_list'][0])
n_state = int((1 + np.sqrt(1 + 8 * diag_length)) / 2)
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
    theta=theta_list, # type: ignore
    phi=phi_replicates,
    neighborhood=ss.spatial.generate_offsets(3, 'manhattan'), 
    random_seed=42,
    )
Sim.initialize()   # Initialize the grid
## Create niches and run simulation
Sim.manual_niche(
    pattern = {'domain_1': ['rectangle', [50, 12, 100, 25]],
               'domain_2': ['rectangle', [50, 37, 100, 25]],
               'domain_3': ['rectangle', [50, 62, 100, 25]]},
)

Sim.gibbs_sampler()    # Gibbs sampling
Sim.density_sampler(density_replicates)  # Cell density of each niche # type: ignore
Sim.perturbation(step = 0.2)
Sim.plot(figsize=(5, 5), dpi=250, size=14)
plt.savefig(os.path.join(script_dir, 'SFig2_panel_A3.png'), bbox_inches='tight', dpi=300)
Sim.plot_niche(figsize=(5, 5), dpi=300)
plt.savefig(os.path.join(script_dir, 'SFig2_panel_A4.png'), bbox_inches='tight', dpi=300)

shape = (100, 100)
param_1 = ss.util.generate_random_parameters(n_group=4, n_state=12, seed=42)

n_group = param_1['n_group']
diag_length = len(param_1['theta_list'][0])
n_state = int((1 + np.sqrt(1 + 8 * diag_length)) / 2)
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
    theta=theta_list, # type: ignore
    phi=phi_replicates,
    neighborhood=ss.spatial.generate_offsets(3, 'manhattan'), 
    random_seed=42,
    )
Sim.initialize() 

## Create a trapezoid boundary in the grid, which takes half of the area
x, y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
# create an ellipse mask for the boundary, which covers (0, 0) to (90, 90)
random_noise = np.random.normal(0, 5, size=shape)
boundary_mask = ((y+x+random_noise) / 200) ** 2 + ((y-x+random_noise)/60) ** 2 > 1
Sim.grid[boundary_mask] = -1

## Create niches and run simulation
Sim.manual_niche(
    pattern = {'domain_1': ['ellipse', [0, 0, 50, 50, 0]],
               'domain_2': ['ellipse', [0, 0, 75, 75, 0]],
               'domain_3': ['ellipse', [0, 0, 90, 90, 0]],
               },
)

Sim.gibbs_sampler()    # Gibbs sampling
Sim.density_sampler(density_replicates)  # Cell density of each niche # type: ignore
Sim.perturbation(step = 0.2)
Sim.plot(figsize=(5, 5), dpi=250, size=14)
plt.savefig(os.path.join(script_dir, 'SFig2_panel_A5.png'), bbox_inches='tight', dpi=300)
Sim.plot_niche(figsize=(5, 5), dpi=300)
plt.savefig(os.path.join(script_dir, 'SFig2_panel_A6.png'), bbox_inches='tight', dpi=300)

########### Panel B:
Sim.create_omics(n_genes=1000, bg_ratio=0.2, bg_param = (1, 0.5), marker_param = (3.5, 1.6), spatial=False)

# Create a dictionary for more efficient bulk updates
gene_meta_updates = {
    'Type_8': 7,
    'Type_5': 5.2,
    'Type_0': 5,
    'Type_7': 3.5,
    'Type_9': 3.5,
    'Type_10': 2.5
}
# Apply updates using .loc for proper DataFrame modification
for gene_type, value in gene_meta_updates.items():
    if gene_type in Sim.gene_meta.columns:
        Sim.gene_meta.loc[0, gene_type] = value

Sim.omics = ss.omics.simOmics(
    omics_meta=Sim.gene_meta, 
    meta=Sim.meta, 
    seed=Sim.seed,
    )

ss.plot.plot_gene(
    Sim.meta, Sim.omics['Gene_0'], size=12, 
    figsize=(4, 4), dpi=300, title='Expression with Gradient (Gene 0)',
    cmap='viridis')
plt.savefig(os.path.join(script_dir, 'SFig2_panel_B.png'), bbox_inches='tight', dpi=300)

########### Panel C:
plt_data = pd.DataFrame({
    'x': Sim.meta['row'],
    'y': Sim.meta['col'],
    'type': Sim.meta['state'],
    'Gene_0': Sim.omics['Gene_0'],
})
# distance along the y=x axis from the origin (non-negative)
plt_data['axis_dist'] = np.abs((plt_data['x'] + plt_data['y']) / np.sqrt(2))
plt.figure(figsize=(5.5,3), dpi=300)
plt.hist2d(
    data=plt_data, 
    x='axis_dist', 
    y='Gene_0', 
    bins=[20, 10]
)
plt.xlabel('Distance along y=x axis')
plt.ylabel('Gene 0 Expression')
plt.title('Gene 0 Expression vs Distance along y=x axis')
plt.savefig(os.path.join(script_dir, 'SFig2_panel_C.png'), bbox_inches='tight', dpi=300)


########### Panel D:
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
    ci_sim = ss.spatial.calculate_interaction_score(sim1.meta['state'], sim1.meta[['row', 'col']], sim1.meta['state'].unique().tolist())

    return [ci_sim]

sim_stats_list = []
Xenium_stats_list = []
merfish_stats_list = []
sccube_stats_list = []
SRTsim_stats_list = []
scMultiSim_stats_list = []
for i in range(8):
    ## SimSpace
    sim_stats = ss_statistics(
        shape=(100, 100), n_group=2, n_state=9, num_iteration=4, n_iter=6, 
        custom_neighbor=ss.spatial.generate_offsets(3, 'manhattan'), seed=i)
    sim_stats_list.append(sim_stats)
    ## Xenium
    tile_path = f'{script_dir}/Panel_D_data/Xenium/tile_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    cells['x_centroid'] = (cells['x_centroid'] - cells['x_centroid'].min()) / (cells['x_centroid'].max() - cells['x_centroid'].min()) * 100
    cells['y_centroid'] = (cells['y_centroid'] - cells['y_centroid'].min()) / (cells['y_centroid'].max() - cells['y_centroid'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['Cluster'], cells[['x_centroid', 'y_centroid']], cells['Cluster'].unique().tolist())
    Xenium_stats_list.append([ci_sim])
    ## scCube
    tile_path = f'{script_dir}/Panel_D_data/sccube/sc_meta_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    cells['point_x'] = (cells['point_x'] - cells['point_x'].min()) / (cells['point_x'].max() - cells['point_x'].min()) * 100
    cells['point_y'] = (cells['point_y'] - cells['point_y'].min()) / (cells['point_y'].max() - cells['point_y'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['Cell_type'], cells[['point_x', 'point_y']], cells['Cell_type'].unique().tolist())
    sccube_stats_list.append([ci_sim])
    ## SRTsim   
    tile_path = f'{script_dir}/Panel_D_data/SRTsim/ref_free_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    cells['x'] = (cells['x'] - cells['x'].min()) / (cells['x'].max() - cells['x'].min()) * 100
    cells['y'] = (cells['y'] - cells['y'].min()) / (cells['y'].max() - cells['y'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['group'], cells[['x', 'y']], cells['group'].unique().tolist())
    SRTsim_stats_list.append([ci_sim])
    ## scMultiSim
    tile_path = f'{script_dir}/Panel_D_data/scMultiSim/spatial_sim_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    cells['x'] = (cells['x'] - cells['x'].min()) / (cells['x'].max() - cells['x'].min()) * 100
    cells['y'] = (cells['y'] - cells['y'].min()) / (cells['y'].max() - cells['y'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['cell.type'], cells[['x', 'y']], cells['cell.type'].unique().tolist())
    scMultiSim_stats_list.append([ci_sim])
    ## MERFISH
    tile_path = f'{script_dir}/Panel_D_data/MERFISH/merfish_data_sample_{i}.csv'
    cells = pd.read_csv(tile_path, index_col=None)
    cells['Centroid_X'] = (cells['Centroid_X'] - cells['Centroid_X'].min()) / (cells['Centroid_X'].max() - cells['Centroid_X'].min()) * 100
    cells['Centroid_Y'] = (cells['Centroid_Y'] - cells['Centroid_Y'].min()) / (cells['Centroid_Y'].max() - cells['Centroid_Y'].min()) * 100
    ci_sim = ss.spatial.calculate_interaction_score(cells['Cell_class'], cells[['Centroid_X', 'Centroid_Y']], cells['Cell_class'].unique().tolist())
    merfish_stats_list.append([ci_sim])

## Cell Type Interaction
sim_ci = np.sort(np.concatenate([item[0].values.reshape(-1).tolist() for item in sim_stats_list]))
xenium_ci = np.sort(np.concatenate([item[0].values.reshape(-1).tolist() for item in Xenium_stats_list]))
sccube_ci = np.sort(np.concatenate([item[0].values.reshape(-1).tolist() for item in sccube_stats_list]))
srt_ci = np.sort(np.concatenate([item[0].values.reshape(-1).tolist() for item in SRTsim_stats_list]))
scMultiSim_ci = np.sort(np.concatenate([item[0].values.reshape(-1).tolist() for item in scMultiSim_stats_list]))
merfish_ci = np.sort(np.concatenate([item[0].values.reshape(-1).tolist() for item in merfish_stats_list]))

# Prepare data for violin plots
ci_labels = ['Xenium', 'MERFISH', 'SimSpace', 'scCube', 'SRTsim', 'scMultiSim']
ci_data = [
    xenium_ci, merfish_ci, sim_ci, sccube_ci, srt_ci, scMultiSim_ci
]
ci_labels = ci_labels

palette = sns.color_palette("muted", 6)
palette = list(palette)
first = np.array(palette[0])
# blend the first color toward a light blue and set as the second color
palette[1] = tuple(np.clip(first * 0.4 + np.array([0.7, 0.85, 1.0]) * 0.6, 0, 1))

plt.figure(figsize=(3.5, 3.5), dpi=300)
sns.violinplot(data=ci_data, palette=palette)
plt.xticks(ticks=np.arange(len(ci_labels)), labels=ci_labels, rotation=30, ha='center')
plt.ylabel('Cell Type Interaction Score')
plt.title('Cell Type Interaction Score Distribution')
plt.tight_layout()
plt.savefig(os.path.join(script_dir, 'SFig2_panel_D.png'), bbox_inches='tight', dpi=300)

########### Panel E:
cell_annotation = pd.read_csv(f'{script_dir}/Panel_E_F_data/cell_annotation.csv', index_col=None)
cell_annotation = cell_annotation.rename(columns={'Barcode': 'cell_id'})
cells = pd.read_csv(f'{script_dir}/Panel_E_F_data/cells.csv', index_col=None)
cells = cells.merge(cell_annotation, on='cell_id', how='left')
cells = cells[cells['Cluster'] != 'Unlabeled']
cells = cells[cells['transcript_counts'] > 0]
cells.reset_index(drop=True, inplace=True)
plt.figure(figsize=(9, 5), dpi=300)
sns.scatterplot(x='x_centroid', y='y_centroid', data=cells, hue='Cluster', s=0.6, palette='tab20', alpha=0.8, edgecolor=None)
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Xenium Breast Tumor Sample')
plt.gca().set_aspect('equal')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=5)
plt.savefig(os.path.join(script_dir, 'SFig2_panel_E.png'), bbox_inches='tight', dpi=300)

########### Panel F:
merfish_data = pd.read_csv(f'{script_dir}/Panel_E_F_data/MERFISH_meta.csv', index_col=None)
plt.figure(figsize=(6, 4), dpi=250)
sns.scatterplot(data=merfish_data, x='Centroid_X', y='Centroid_Y', hue='Cell_class', s=9)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(script_dir, 'SFig2_panel_F.png'), bbox_inches='tight', dpi=300)