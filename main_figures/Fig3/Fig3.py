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

########## Panel B:
cells = pd.read_csv(f"{script_dir}/Panel_B_C_D_data/Xenium_reference_metadata.csv")
cell_counts = cells['Cluster'].value_counts()
cells = cells[cells['Cluster'].isin(cell_counts[cell_counts >= 10].index)]
ranked_cell_types = {cell_type: rank for rank, cell_type in enumerate(cells['Cluster'].value_counts().index, 1)}
cells['celltype_rank'] = cells['Cluster'].map(ranked_cell_types)

cells['x_centroid'] = 100 * (cells['x_centroid'] - cells['x_centroid'].min()) / (cells['x_centroid'].max() - cells['x_centroid'].min())
cells['y_centroid'] = 100 * (cells['y_centroid'] - cells['y_centroid'].min()) / (cells['y_centroid'].max() - cells['y_centroid'].min())

plt.figure(figsize=(5, 3.5), dpi = 250)
sns.scatterplot(x='x_centroid', y='y_centroid', data=cells, hue='celltype_rank', palette='tab20', legend=True, s=10)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
# Map celltype_rank to Cluster labels for the legend
handles, labels = plt.gca().get_legend_handles_labels()
cluster_labels = {str(rank): cluster for cluster, rank in ranked_cell_types.items()}
new_labels = [cluster_labels.get(label, label) for label in labels]
# plt.legend(handles, new_labels, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend([], [], frameon=False)
plt.title('Xenium Breast Tumor Sample')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig3_panel_B1.png', bbox_inches='tight', dpi=300)


sccube_cells = pd.read_csv(f"{script_dir}/Panel_B_C_D_data/Xenium_Breast_scCube.csv",
                    index_col=0)
sccube_cell_counts = sccube_cells['Cell_type'].value_counts()
sccube_cells = sccube_cells[sccube_cells['Cell_type'].isin(sccube_cell_counts[sccube_cell_counts >= 10].index)]
ranked_cell_types = {cell_type: rank for rank, cell_type in enumerate(sccube_cells['Cell_type'].value_counts().index, 1)}
sccube_cells['celltype_rank'] = sccube_cells['Cell_type'].map(ranked_cell_types)

sccube_cells['point_x'] = 100 * (sccube_cells['point_x'] - sccube_cells['point_x'].min()) / (sccube_cells['point_x'].max() - sccube_cells['point_x'].min())
sccube_cells['point_y'] = 100 * (sccube_cells['point_y'] - sccube_cells['point_y'].min()) / (sccube_cells['point_y'].max() - sccube_cells['point_y'].min())

plt.figure(figsize=(5, 3.5), dpi = 250)
sns.scatterplot(x='point_x', y='point_y', data=sccube_cells, hue='celltype_rank', palette='tab20', legend=True, s=10)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend([], [], frameon=False)
plt.title('scCube Simulation')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig3_panel_B2.png', bbox_inches='tight', dpi=300)

shape = (100, 100)
sim = ss.util.sim_from_json(
    input_file=f'{script_dir}/Panel_B_C_D_data/simspace_fitted_params.json',
    shape=shape,
    num_iteration=4,
    n_iter=6,
)
simulation = sim.meta.copy()
simulation.columns = ['celltype', 'x', 'y', 'spatial_domain', 'celltype_rank']
plt.figure(figsize=(5, 3.5), dpi = 250)
sns.scatterplot(x='x', y='y', data=simulation, hue='celltype_rank', palette='tab20', legend=True, s=10)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend([], [], frameon=False)
plt.title('SimSpace Simulation')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig3_panel_B3.png', bbox_inches='tight', dpi=300)

########## Panel D:
shape = (100, 100)
sim = ss.util.sim_from_json(
    input_file=f'{script_dir}/Panel_B_C_D_data/simspace_fitted_params.json',
    shape=shape,
    num_iteration=4,
    n_iter=6,
)

simulation = sim.meta.copy()
simulation.columns = ['celltype', 'x', 'y', 'spatial_domain', 'celltype_rank']
plt.figure(figsize=(5, 3.5), dpi = 200)
sns.scatterplot(x='x', y='y', data=simulation, hue='spatial_domain', palette='muted', legend=True, s=10)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('SimSpace Simulation')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig3_panel_D1.png', bbox_inches='tight', dpi=300)

bank_ss_domain = pd.read_csv(f'{script_dir}/Panel_B_C_D_data/BANKSY_simspace_domain.csv')
palette = sns.color_palette("muted", 10)
palette[0] = palette[3]
palette[1] = palette[2]

plt.figure(figsize=(5, 3.5), dpi = 200)
sns.scatterplot(x='row', y='col', data=bank_ss_domain, hue='Tumor_label', palette=palette, legend=True, s=10)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
plt.legend([], [], frameon=False)
plt.title('SimSpace BANKSY Tumor Domain')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig3_panel_D2.png', bbox_inches='tight', dpi=300)


plt.figure(figsize=(5, 3.5), dpi = 200)
sns.scatterplot(x='row', y='col', data=bank_ss_domain, hue='UMAP1_cluster', palette='tab20', legend=True, s=10)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
plt.legend([], [], frameon=False)
plt.title('SimSpace BANKSY Clustering')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig3_panel_D3.png', bbox_inches='tight', dpi=300)

########## Panel E:
bank_xenium_domain = pd.read_csv(f'{script_dir}/Panel_B_C_D_data/BANKSY_xenium_domain.csv')

# Compute domain-wise cell type composition for SimSpace (ss) and Xenium
ct_order = pd.concat([bank_ss_domain['fitted_celltype'], bank_xenium_domain['fitted_celltype']]) \
             .value_counts().index.tolist()

# Cross-tab and normalize by domain (rows sum to 1)
ct_ss = pd.crosstab(bank_ss_domain['Tumor_label'], bank_ss_domain['fitted_celltype'])
comp_ss = ct_ss.reindex(columns=ct_order, fill_value=0).div(ct_ss.sum(axis=1), axis=0)

ct_xen = pd.crosstab(bank_xenium_domain['Tumor_label'], bank_xenium_domain['fitted_celltype'])
comp_xen = ct_xen.reindex(columns=ct_order, fill_value=0).div(ct_xen.sum(axis=1), axis=0)

# Common scale for colorbar
vmax = max(comp_ss.max().max(), comp_xen.max().max())

# Create percent annotation matrices
annot_ss  = comp_ss.map(lambda v: f"{v*100:.1f}%").to_numpy()
annot_xen = comp_xen.map(lambda v: f"{v*100:.1f}%").to_numpy()

# Plot side-by-side heatmaps
fig, axes = plt.subplots(1, 2, figsize=(6, 6), dpi=250)
sns.heatmap(comp_ss.T, ax=axes[0], cmap='rocket_r', vmin=0, vmax=vmax,
            annot=annot_ss.T, fmt='', cbar=False, cbar_kws={'label': 'Fraction'})
axes[0].set_title('SimSpace')
axes[0].set_xlabel('Tumor Domain')
axes[0].set_ylabel('Cell Type')
axes[0].set_aspect('equal')

sns.heatmap(comp_xen.T, ax=axes[1], cmap='rocket_r', vmin=0, vmax=vmax,
            annot=annot_xen.T, fmt='', cbar=True, cbar_kws={'label': 'Fraction'})
axes[1].set_title('Xenium')
axes[1].set_xlabel('Tumor Domain')
axes[1].set_ylabel(None)
axes[1].set_yticklabels([])
axes[1].set_aspect('equal')

plt.tight_layout()
plt.savefig(f'{script_dir}/Fig3_panel_E.png', bbox_inches='tight', dpi=300)