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
import subprocess

######## Panel A:
shape = (100, 100)
# Generate a simulation space
sim1 = ss.util.sim_from_json(
    input_file=f'{script_dir}/Panel_A_data/simspace_fitted_params.json',
    shape=shape,
    num_iteration=4,
    n_iter=6,
)
sim1.update_seed(seed=1)
simulation = sim1.meta.copy()
simulation.columns = ['celltype', 'x', 'y', 'spatial_domain', 'celltype_rank']
tmp = simulation[simulation['celltype'] == -1].index
simulation = simulation.drop(tmp)
simulation = simulation.reset_index(drop=True)
simulation['celltype'] = simulation['celltype'].astype(int)
cell_counts = simulation['celltype'].value_counts()

plt.figure(figsize=(5, 5), dpi = 250)
sns.scatterplot(x='y', y='x', data=simulation, hue='celltype_rank', palette='tab20', legend=True, s=10)
plt.gca().set_aspect('equal')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig4_panel_A1.png', dpi=250)

sim1.update_seed(seed=1)
sim1.create_omics(n_genes=1000, bg_ratio=0.2)
cmap = sns.color_palette('tab20', n_colors=sim1.num_states)

kernel = (3, 3)
spot_meta, spot_omics = ss.util.convolve(sim1, kernel=kernel, scale=2)
# print(state_colors[0])
fig, ax = plt.subplots()
fig.set_size_inches(4, 3)
fig.set_dpi(150)
ax.set_aspect('equal')
for i in range(len(spot_meta)):
    centroid_x = spot_meta.iloc[i]['col']
    centroid_y = spot_meta.iloc[i]['row']
    state_proportions = spot_meta.iloc[i][2:]
    state_proportions = state_proportions[state_proportions > 0]
    state_proportions = state_proportions / state_proportions.sum()
    wedges, texts = ax.pie(state_proportions, 
                            colors=[cmap[state-1] for state in state_proportions.index], 
                            startangle=90, 
                            radius=kernel[0]/3, 
                            center=(centroid_x, centroid_y), 
                            frame=True,
                            )
plt.savefig(f'{script_dir}/Fig4_panel_A2.png', dpi=250, bbox_inches='tight')

kernel = (5, 5)
spot_meta, spot_omics = ss.util.convolve(sim1, kernel=kernel, scale=2)
# print(state_colors[0])
fig, ax = plt.subplots()
fig.set_size_inches(4, 3)
fig.set_dpi(150)
ax.set_aspect('equal')
for i in range(len(spot_meta)):
    centroid_x = spot_meta.iloc[i]['col']
    centroid_y = spot_meta.iloc[i]['row']
    state_proportions = spot_meta.iloc[i][2:]
    state_proportions = state_proportions[state_proportions > 0]
    state_proportions = state_proportions / state_proportions.sum()
    wedges, texts = ax.pie(state_proportions, 
                            colors=[cmap[state-1] for state in state_proportions.index], 
                            startangle=90, 
                            radius=kernel[0]/3, 
                            center=(centroid_x, centroid_y), 
                            frame=True,
                            )
plt.savefig(f'{script_dir}/Fig4_panel_A3.png', dpi=250, bbox_inches='tight')

kernel = (7, 7)
spot_meta, spot_omics = ss.util.convolve(sim1, kernel=kernel, scale=2)
# print(state_colors[0])
fig, ax = plt.subplots()
fig.set_size_inches(4, 3)
fig.set_dpi(150)
ax.set_aspect('equal')
for i in range(len(spot_meta)):
    centroid_x = spot_meta.iloc[i]['col']
    centroid_y = spot_meta.iloc[i]['row']
    state_proportions = spot_meta.iloc[i][2:]
    state_proportions = state_proportions[state_proportions > 0]
    state_proportions = state_proportions / state_proportions.sum()
    wedges, texts = ax.pie(state_proportions, 
                            colors=[cmap[state-1] for state in state_proportions.index], 
                            startangle=90, 
                            radius=kernel[0]/3, 
                            center=(centroid_x, centroid_y), 
                            frame=True,
                            )
plt.savefig(f'{script_dir}/Fig4_panel_A4.png', dpi=250, bbox_inches='tight')

kernel = (10, 10)
spot_meta, spot_omics = ss.util.convolve(sim1, kernel=kernel, scale=2)
# print(state_colors[0])
fig, ax = plt.subplots()
fig.set_size_inches(4, 3)
fig.set_dpi(150)
ax.set_aspect('equal')
for i in range(len(spot_meta)):
    centroid_x = spot_meta.iloc[i]['col']
    centroid_y = spot_meta.iloc[i]['row']
    state_proportions = spot_meta.iloc[i][2:]
    state_proportions = state_proportions[state_proportions > 0]
    state_proportions = state_proportions / state_proportions.sum()
    wedges, texts = ax.pie(state_proportions, 
                            colors=[cmap[state-1] for state in state_proportions.index], 
                            startangle=90, 
                            radius=kernel[0]/3, 
                            center=(centroid_x, centroid_y), 
                            frame=True,
                            )
plt.savefig(f'{script_dir}/Fig4_panel_A5.png', dpi=250, bbox_inches='tight')


######## Panel B, C:
tab10_palette = sns.color_palette("tab10", 10)
# Define the variable order
method_names = ['c2l', 'RCTD_full', 'RCTD_multi', 'CARD', 'Seurat', 'spatialDWLS']
# Create a mapping from method name to color
method_color_map = {name: tab10_palette[i] for i, name in enumerate(method_names)}

# PCC, RMSE, JSD results aggregation
pcc_res_df = pd.DataFrame()
rmse_res_df = pd.DataFrame()
jsd_res_df = pd.DataFrame()

kernel_size = (5, 7, 10, 15)
n_niche = (2, 3)
n_state = (8, 9)
seeds = (0, 1)

for kernel in kernel_size:
    for niches in n_niche:
        for states in n_state:
            for seed in seeds:
                pcc_res = pd.read_csv(f"{script_dir}/Panel_B_C_D_data/pcc/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv")
                pcc_res_long = pcc_res.melt(var_name='Metric', value_name='Value')
                pcc_res_long['Kernel'] = kernel
                pcc_res_long['Niche'] = niches
                pcc_res_long['State'] = states
                pcc_res_long['Seed'] = seed
                pcc_res_df = pd.concat([pcc_res_df, pcc_res_long], ignore_index=True)

                rmse_res = pd.read_csv(f"{script_dir}/Panel_B_C_D_data/rmse/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv")
                rmse_res_long = rmse_res.melt(var_name='Metric', value_name='Value')
                rmse_res_long['Kernel'] = kernel
                rmse_res_long['Niche'] = niches
                rmse_res_long['State'] = states
                rmse_res_long['Seed'] = seed
                rmse_res_df = pd.concat([rmse_res_df, rmse_res_long], ignore_index=True)

                jsd_res = pd.read_csv(f"{script_dir}/Panel_B_C_D_data/jaccard/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv")
                jsd_res_long = jsd_res.melt(var_name='Metric', value_name='Value')
                jsd_res_long['Kernel'] = kernel
                jsd_res_long['Niche'] = niches
                jsd_res_long['State'] = states
                jsd_res_long['Seed'] = seed
                jsd_res_df = pd.concat([jsd_res_df, jsd_res_long], ignore_index=True)

pcc_res_average = pcc_res_df.groupby(['Kernel', 'Niche', 'State', 'Metric', 'Seed']).agg({'Value': ['mean', 'std']}).reset_index()
rmse_res_average = rmse_res_df.groupby(['Kernel', 'Niche', 'State', 'Metric', 'Seed']).agg({'Value': ['mean', 'std']}).reset_index()
jsd_res_average = jsd_res_df.groupby(['Kernel', 'Niche', 'State', 'Metric', 'Seed']).agg({'Value': ['mean', 'std']}).reset_index()
# Ensure columns are accessible (flatten MultiIndex if needed)
df_avg = pcc_res_average.copy()
df_avg.columns = ['Kernel', 'Niche', 'State', 'Metric', 'Seed', 'Mean', 'Std']

plt.figure(figsize=(5, 3), dpi=300)
sns.lineplot(
    data=df_avg,
    x='Kernel',
    y='Mean',
    hue='Metric',
    style='Metric',
    palette=method_color_map,
    markers=True,
    dashes=True
)
plt.ylabel('PCC')
plt.xlabel('Kernel Size')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title ='Metric')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig4_panel_B1.png', dpi=300)

plt.figure(figsize=(6.5, 5), dpi=300)
sns.boxplot(
    data=pcc_res_df[pcc_res_df['Kernel'] == 10],
    y='Metric',
    x='Value',
    hue='Metric',
    fliersize=0,
    palette=method_color_map,
)
plt.xlim(-0.05, 1.03)
plt.title(f'PCC Distribution by Metric (Kernel=10)')
plt.ylabel('PCC Value')
plt.xlabel('Metric')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title ='Metric')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig4_panel_C1.png', dpi=300)

df_avg = rmse_res_average.copy()
df_avg.columns = ['Kernel', 'Niche', 'State', 'Metric', 'Seed', 'Mean', 'Std']
plt.figure(figsize=(5, 3), dpi=300)
sns.lineplot(
    data=df_avg,
    x='Kernel',
    y='Mean',
    hue='Metric',
    style='Metric',
    palette=method_color_map,
    markers=True,
    dashes=True
)
plt.ylabel('RMSE')
plt.xlabel('Kernel Size')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title ='Metric')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig4_panel_B2.png', dpi=300)

plt.figure(figsize=(6.5, 5), dpi=300)
sns.boxplot(
    data=rmse_res_df[rmse_res_df['Kernel'] == 10],
    y='Metric',
    x='Value',
    hue='Metric',
    fliersize=0,
    palette=method_color_map,
)
plt.xlim(-0.05, 0.4)
plt.title(f'RMSE Distribution by Metric (Kernel=10)')
plt.ylabel('RMSE Value')
plt.xlabel('Metric')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title ='Metric')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig4_panel_C2.png', dpi=300)

df_avg = jsd_res_average.copy()
df_avg.columns = ['Kernel', 'Niche', 'State', 'Metric', 'Seed', 'Mean', 'Std']
plt.figure(figsize=(5, 3), dpi=300)
sns.lineplot(
    data=df_avg,
    x='Kernel',
    y='Mean',
    hue='Metric',
    style='Metric',
    palette=method_color_map,
    markers=True,
    dashes=True
)
plt.ylabel('JSD')
plt.xlabel('Kernel Size')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title ='Metric')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig4_panel_B3.png', dpi=300)

plt.figure(figsize=(6.5, 5), dpi=300)
sns.boxplot(
    data=jsd_res_df[jsd_res_df['Kernel'] == 10],
    y='Metric',
    x='Value',
    hue='Metric',
    fliersize=0,
    palette=method_color_map,
)
plt.xlim(-0.05, 1.03)
plt.title(f'JSD Distribution by Metric (Kernel=10)')
plt.ylabel('JSD Value')
plt.xlabel('Metric')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title ='Metric')
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig4_panel_C3.png', dpi=300)


######## Panel D:
# PCC
pcc_res_df = pd.DataFrame()

kernel_size = (5, 7, 10, 15)
n_niche = (2, 3)
n_state = (8, 9)
seeds = (0, 1)

for kernel in kernel_size:
    for niches in n_niche:
        for states in n_state:
            for seed in seeds:
                pcc_res = pd.read_csv(f"{script_dir}/Panel_B_C_D_data/pcc/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv")
                pcc_res_long = pcc_res.melt(var_name='Metric', value_name='Value')
                pcc_res_long['Kernel'] = kernel
                pcc_res_long['Niche'] = niches
                pcc_res_long['State'] = states
                pcc_res_long['Seed'] = seed
                pcc_res_df = pd.concat([pcc_res_df, pcc_res_long], ignore_index=True)

n_niche = (3, )
n_state = (11, )
seeds = (0, 1, 3, 4)

for kernel in kernel_size:
    for niches in n_niche:
        for states in n_state:
            for seed in seeds:
                pcc_res = pd.read_csv(f"{script_dir}/Panel_B_C_D_data/pcc/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv")
                pcc_res_long = pcc_res.melt(var_name='Metric', value_name='Value')
                pcc_res_long['Kernel'] = kernel
                pcc_res_long['Niche'] = niches
                pcc_res_long['State'] = states
                pcc_res_long['Seed'] = seed
                pcc_res_df = pd.concat([pcc_res_df, pcc_res_long], ignore_index=True)

n_state = (14, )
seeds = (0, 1, 2, 3)

for kernel in kernel_size:
    for niches in n_niche:
        for states in n_state:
            for seed in seeds:
                pcc_res = pd.read_csv(f"{script_dir}/Panel_B_C_D_data/pcc/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv")
                pcc_res_long = pcc_res.melt(var_name='Metric', value_name='Value')
                pcc_res_long['Kernel'] = kernel
                pcc_res_long['Niche'] = niches
                pcc_res_long['State'] = states
                pcc_res_long['Seed'] = seed
                pcc_res_df = pd.concat([pcc_res_df, pcc_res_long], ignore_index=True)

pcc_res_average = pcc_res_df.groupby(['Kernel', 'Niche', 'State', 'Metric', 'Seed']).agg({'Value': ['mean', 'std']}).reset_index()
# Ensure columns are accessible (flatten MultiIndex if needed)
df_avg = pcc_res_average.copy()
df_avg.columns = ['Kernel', 'Niche', 'State', 'Metric', 'Seed', 'Mean', 'Std']


kernel_sizes = (5, 7, 10, 15)
for i, kernel in enumerate(kernel_sizes):
    plt.figure(figsize=(5, 3.2), dpi=300)
    sns.lineplot(
        data=df_avg[df_avg['Kernel'] == kernel],
        x='State',
        y='Mean',
        hue='Metric',
        style='Metric',
        palette=method_color_map,
        markers=True,
        dashes=True
    )
    plt.ylabel('PCC')
    plt.xlabel('# Simulated Cell Types')
    plt.title('Kernel Size: {}'.format(kernel))
    plt.ylim(0.3, 1)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title ='Metric')
    plt.tight_layout()
    plt.savefig(f'{script_dir}/Fig4_panel_D{i+1}.png', dpi=300)
