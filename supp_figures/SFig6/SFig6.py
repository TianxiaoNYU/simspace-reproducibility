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


plt_kernel_size = (5, 7, 15)
######## Panel A:
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
                pcc_res = pd.read_csv(f"{script_dir}/Panel_A_B_C_data/pcc/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv")
                pcc_res_long = pcc_res.melt(var_name='Metric', value_name='Value')
                pcc_res_long['Kernel'] = kernel
                pcc_res_long['Niche'] = niches
                pcc_res_long['State'] = states
                pcc_res_long['Seed'] = seed
                pcc_res_df = pd.concat([pcc_res_df, pcc_res_long], ignore_index=True)

                rmse_res = pd.read_csv(f"{script_dir}/Panel_A_B_C_data/rmse/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv")
                rmse_res_long = rmse_res.melt(var_name='Metric', value_name='Value')
                rmse_res_long['Kernel'] = kernel
                rmse_res_long['Niche'] = niches
                rmse_res_long['State'] = states
                rmse_res_long['Seed'] = seed
                rmse_res_df = pd.concat([rmse_res_df, rmse_res_long], ignore_index=True)

                jsd_res = pd.read_csv(f"{script_dir}/Panel_A_B_C_data/jaccard/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv")
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

for i, plt_kernel in enumerate(plt_kernel_size):
    plt.figure(figsize=(6.5, 5), dpi=300)
    sns.boxplot(
        data=pcc_res_df[pcc_res_df['Kernel'] == plt_kernel],
        y='Metric',
        x='Value',
        hue='Metric',
        fliersize=0,
        palette=method_color_map,
    )
    plt.xlim(-0.05, 1.03)
    plt.title(f'PCC Distribution by Metric (Kernel={plt_kernel})')
    plt.ylabel('PCC Value')
    plt.xlabel('Metric')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title ='Metric')
    plt.tight_layout()
    plt.savefig(f'{script_dir}/Fig4_panel_A{i+1}.png', dpi=300)

    plt.figure(figsize=(6.5, 5), dpi=300)
    sns.boxplot(
        data=rmse_res_df[rmse_res_df['Kernel'] == plt_kernel],
        y='Metric',
        x='Value',
        hue='Metric',
        fliersize=0,
        palette=method_color_map,
    )
    plt.xlim(-0.05, 0.4)
    plt.title(f'RMSE Distribution by Metric (Kernel={plt_kernel})')
    plt.ylabel('RMSE Value')
    plt.xlabel('Metric')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title ='Metric')
    plt.tight_layout()
    plt.savefig(f'{script_dir}/Fig4_panel_B{i+1}.png', dpi=300)
    plt.figure(figsize=(6.5, 5), dpi=300)
    sns.boxplot(
        data=jsd_res_df[jsd_res_df['Kernel'] == plt_kernel],
        y='Metric',
        x='Value',
        hue='Metric',
        fliersize=0,
        palette=method_color_map,
    )
    plt.xlim(-0.05, 1.03)
    plt.title(f'JSD Distribution by Metric (Kernel={plt_kernel})')
    plt.ylabel('JSD Value')
    plt.xlabel('Metric')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title ='Metric')
    plt.tight_layout()
    plt.savefig(f'{script_dir}/Fig4_panel_C{i+1}.png', dpi=300)