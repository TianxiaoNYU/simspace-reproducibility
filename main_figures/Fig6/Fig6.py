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

Marker_list = [
    'DAPI', 'MPO', 'Ecadherin', 'PDL1', 'CD163', 'PD1', 'CD47', 
    'GAL3', 'PARP1', 'LAG3', 'CD4', 'PI3KCA', 'TIM3', 'CD68', 
    'ER', 'PR', 'MSH2', 'CD8', 'MSH6', 'bCatenin1', 'HLAABC', 
    'MLH1', 'Ki67', 'CD20', 'ARID1A', 'IFNG', 'CD31', 'PMS', 
    'CD44', 'PanCytokeratin', 'CD3e']

########## Panel A:
sim = ss.util.sim_from_json(
    input_file=f'{script_dir}/Panel_A_B_C_data/simspace_fitted.json',
    shape=(50, 50),
    num_iteration=4,
    n_iter=6,
)
sim.update_seed(seed=1)
sim.fit_scdesign(
    f'{script_dir}/Panel_A_B_C_data/CODEX_expr_mat.csv',
    f'{script_dir}/Panel_A_B_C_data/CODEX_meta.csv',
    'phenotype',
    'X_centroid',
    'Y_centroid',
    seed=0,
)

cells = pd.read_csv(f"{script_dir}/Panel_A_B_C_data/CODEX_ref.csv")
cell_counts = cells['cell_type'].value_counts()
ranked_cell_types = {cell_type: rank for rank, cell_type in enumerate(cells['cell_type'].value_counts().index, 1)}
cells['celltype_rank'] = cells['cell_type'].map(ranked_cell_types)
cells['x_centroid'] = 50 * (cells['centroid_x'] - cells['centroid_x'].min()) / (cells['centroid_x'].max() - cells['centroid_x'].min())
cells['y_centroid'] = 50 * (cells['centroid_y'] - cells['centroid_y'].min()) / (cells['centroid_y'].max() - cells['centroid_y'].min())
cells['cell_type'] = cells['cell_type'].replace('CD4+ T cells', 'Helper T cells')
# Get unique cell types from both sim.meta['fitted_celltype'] and cells['phenotype']
sim_celltypes = sim.meta['fitted_celltype'].unique()
# cells_celltypes = cells['phenotype'].unique()
cells_celltypes = cells['cell_type'].unique()
all_celltypes = pd.unique(np.concatenate([sim_celltypes, cells_celltypes]))

# Create a consistent palette mapping
palette = dict(zip(
    all_celltypes,
    sns.color_palette('tab20', n_colors=len(all_celltypes))
))

plt.figure(figsize=(6, 3.5), dpi=250)
sns.scatterplot(
    x='row', y='col', data=sim.meta,
    hue='fitted_celltype', palette=palette, legend=True, s=15
)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('SimSpace Simulation')
plt.tight_layout()
plt.show()

plt.figure(figsize=(6, 3.5), dpi=250)
sns.scatterplot(
    x='x_centroid', y='y_centroid', data=cells,
    hue='cell_type', palette=palette, legend=True, s=15
)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('CODEX Tumor Sample')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig6_Panel_A1.png')


marker_expr = pd.read_csv(f"{script_dir}/Panel_A_B_C_data/CODEX_expr_mat.csv", index_col=0)
expr_cells = cells.copy()
expr_cells['PD1'] = marker_expr['PD1'].values
cmap = sns.color_palette('flare', as_cmap=True)
plt.figure(figsize=(6, 3.5), dpi=250)
sns.scatterplot(
    x='x_centroid', y='y_centroid', data=expr_cells,
    hue='PD1', legend=True, s=15, palette=cmap
)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('PD1 Marker Expression')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig6_Panel_A2.png')

ss.plot.plot_gene(
    sim.meta, sim.omics['PD1'], size=14, 
    figsize=(5,5), dpi=300, title='PD1 Marker Expression',)
plt.savefig(f'{script_dir}/Fig6_Panel_A3.png', bbox_inches='tight', dpi=300)


########## Panel B:
state_to_celltype = dict(zip(sim.meta.state.unique(), sim.meta.fitted_celltype.unique()))

hm_data = sim.theta[0].copy()
hm_data = pd.DataFrame(hm_data, index=range(hm_data.shape[0]), columns=range(hm_data.shape[1]))
hm_data.index = hm_data.index.map(state_to_celltype)
hm_data.columns = hm_data.columns.map(state_to_celltype)

hm_data = hm_data.reindex(sorted(hm_data.columns), axis=1)
hm_data = hm_data.reindex(sorted(hm_data.index), axis=0)

# make all diagonal elements NA for hm_data
hm_arr = hm_data.to_numpy(copy=True)
np.fill_diagonal(hm_arr, np.nan)

# plot na as gray
plt.figure(figsize=(7,7), dpi=250)
sns.heatmap(hm_arr, annot=True, cmap='viridis', vmin=0, fmt='.2g', 
            mask=np.eye(hm_arr.shape[0], dtype=bool), 
            cbar_kws={'label': 'Interaction Strength'}, square=True, 
            linewidths=0, linecolor='black')
plt.xlabel('Cell Type')
plt.ylabel('Cell Type')
plt.gca().collections[0].colorbar.remove()  # type: ignore
plt.gca().set_aspect('equal')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig6_Panel_B1.png')

# Get unique cell types from both sim.meta['fitted_celltype'] and cells['phenotype']
sim_celltypes = sim.meta['fitted_celltype'].unique()
# cells_celltypes = cells['phenotype'].unique()
cells_celltypes = cells['cell_type'].unique()
all_celltypes = pd.unique(np.concatenate([sim_celltypes, cells_celltypes]))

# Create a consistent palette mapping
palette = dict(zip(
    all_celltypes,
    sns.color_palette('tab20', n_colors=len(all_celltypes))
))

# Only show Cytotoxic T cells and Helper T cells in color, others in transparent gray
highlight_types = ['Cytotoxic T cells', 'Helper T cells']
highlight_colors = {k: palette[k] for k in highlight_types}
other_color = (0.5, 0.5, 0.5, 0.1)  # transparent gray

# Build new palette: keep original color for highlights, gray for others
plot_palette = {k: highlight_colors[k] if k in highlight_types else other_color for k in palette}

plt.figure(figsize=(5, 3.5), dpi=250)
sns.scatterplot(
    x='row', y='col', data=sim.meta,
    hue='fitted_celltype', palette=plot_palette, legend=True, s=15
)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend([], [], frameon=False)
plt.title('SimSpace Simulation')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig6_Panel_B2.png')

plt.figure(figsize=(5, 3.5), dpi=250)
sns.scatterplot(
    x='x_centroid', y='y_centroid', data=cells,
    hue='cell_type', palette=plot_palette, legend=True, s=15
)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
plt.legend([], [], frameon=False)
plt.title('CODEX Tumor Sample')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig6_Panel_B3.png')

# Only show Endothelial cells and Stromal cells in color, others in transparent gray
highlight_types = ['Endothelial cells', 'Stromal cells (undefined)']
highlight_colors = {k: palette[k] for k in highlight_types}
other_color = (0.5, 0.5, 0.5, 0.1)  # transparent gray

# Build new palette: keep original color for highlights, gray for others
plot_palette = {k: highlight_colors[k] if k in highlight_types else other_color for k in palette}

plt.figure(figsize=(5, 3.5), dpi=250)
sns.scatterplot(
    x='row', y='col', data=sim.meta,
    hue='fitted_celltype', palette=plot_palette, legend=True, s=15
)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend([], [], frameon=False)
plt.title('SimSpace Simulation')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig6_Panel_B4.png')

plt.figure(figsize=(5, 3.5), dpi=250)
sns.scatterplot(
    x='x_centroid', y='y_centroid', data=cells,
    hue='cell_type', palette=plot_palette, legend=True, s=15
)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
plt.legend([], [], frameon=False)
plt.title('CODEX Tumor Sample')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig6_Panel_B5.png')

########## Panel C:
highlight_types = ['Cytotoxic T cells', 'Helper T cells']
highlight_colors = {k: palette[k] for k in highlight_types}
other_color = (0.5, 0.5, 0.5, 0.1)  # transparent gray

# Build new palette: keep original color for highlights, gray for others
plot_palette = {k: highlight_colors[k] if k in highlight_types else other_color for k in palette}
fig, axes = plt.subplots(1, 7, figsize=(20, 4), dpi=250)
axes = axes.flatten()

for i in range(7):
    sim_p = ss.util.sim_from_json(
        input_file=f'{script_dir}/Panel_A_B_C_data/perturb/ss_perturb_{i+1}.json',
        shape=(50, 50),
        num_iteration=4,
        n_iter=6,
    )
    state_to_fitted_celltype = sim.meta[['state', 'fitted_celltype']].drop_duplicates().set_index('state')['fitted_celltype'].to_dict()
    sim_p.meta['fitted_celltype'] = sim_p.meta['state'].map(state_to_fitted_celltype)

    ax = axes[i]
    sns.scatterplot(
        x='row', y='col', data=sim_p.meta,
        hue='fitted_celltype', palette=plot_palette, legend=False, s=15, ax=ax
    )
    ax.set_xlabel('X')
    if i == 0:
        ax.set_ylabel('Y')
    else:
        ax.set_ylabel('')
        ax.set_yticklabels([])
    ax.set_aspect('equal')
    ax.set_title(f'SimSpace Simulation {i+1}')
    ax.set_xticks([])
    ax.set_yticks([])

plt.tight_layout()
plt.savefig(f'{script_dir}/Fig6_Panel_C.png')