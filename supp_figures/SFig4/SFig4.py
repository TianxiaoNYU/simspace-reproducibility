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

########### Panel C:
bank_xenium_domain = pd.read_csv(f'{script_dir}/Panel_C_D_data/BANKSY_xenium_domain.csv')
plt.figure(figsize=(5, 3.5), dpi = 200)
sns.scatterplot(x='row', y='col', data=bank_xenium_domain, hue='UMAP1_cluster', palette='tab20', legend=True, s=10)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend([], [], frameon=False)
plt.title('Xenium BANKSY Clustering')
plt.tight_layout()
plt.savefig(f'{script_dir}/SFig4_panel_C1.png')

palette = sns.color_palette("muted", 10)
palette[0] = palette[2]
palette[1] = palette[3]
plt.figure(figsize=(5, 3.5), dpi = 200)
sns.scatterplot(x='row', y='col', data=bank_xenium_domain, hue='Tumor_label', palette=palette, legend=True, s=10)
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend([], [], frameon=False)
plt.title('Xenium BANKSY Clustering by Tumor')
plt.tight_layout()
plt.savefig(f'{script_dir}/SFig4_panel_C2.png')

########### Panel D:
ss_count = pd.read_csv(f'{script_dir}/Panel_C_D_data/simspace_fitted_count.csv', index_col=0)
sccube_count = pd.read_csv(f'{script_dir}/Panel_C_D_data/sccube_fitted_count.csv', index_col=0)
sccube_count = sccube_count.T
xenium_count = pd.read_csv(f'{script_dir}/Panel_C_D_data/Xenium_reference_count.csv', index_col=0)
xenium_count = xenium_count.T

# Compute fraction zeros, mean and variance per gene for the three datasets and plot relationships

# Use the intersection of genes present in all three datasets
common_genes = list(set(ss_count.columns) & set(sccube_count.columns) & set(xenium_count.columns))
ss_sub = ss_count[common_genes]
sccube_sub = sccube_count[common_genes]
xenium_sub = xenium_count[common_genes]

def gene_metrics(df, label):
    # df: genes in columns, cells/samples in rows
    means = df.mean(axis=0)
    variances = df.var(axis=0, ddof=1)
    frac_zero = (df == 0).sum(axis=0) / df.shape[0]
    out = pd.DataFrame({
        'gene': means.index,
        'mean': means.values,
        'variance': variances.values,
        'frac_zero': frac_zero.values,
        'dataset': label
    }).set_index('gene')
    return out

metrics_ss = gene_metrics(ss_sub, 'SimSpace')
metrics_sccube = gene_metrics(sccube_sub, 'scCube')
metrics_xenium = gene_metrics(xenium_sub, 'Xenium')

metrics = pd.concat([metrics_ss, metrics_sccube, metrics_xenium])

# Summary correlations (mean vs variance, mean vs frac_zero) per dataset
corrs = {}
for name, grp in metrics.groupby('dataset'):
    # handle potential constant series by checking std
    mean_var_corr = np.nan
    mean_fzero_corr = np.nan
    if grp['mean'].std() > 0 and grp['variance'].std() > 0:
        mean_var_corr = np.corrcoef(grp['mean'], grp['variance'])[0,1]
    if grp['mean'].std() > 0 and grp['frac_zero'].std() > 0:
        mean_fzero_corr = np.corrcoef(grp['mean'], grp['frac_zero'])[0,1]
    corrs[name] = {'mean_vs_variance': mean_var_corr, 'mean_vs_frac_zero': mean_fzero_corr}

print("Correlations (mean vs variance, mean vs fraction zeros):")
for k,v in corrs.items():
    print(f"  {k}: mean-variance {v['mean_vs_variance']:.4f}, mean-frac_zero {v['mean_vs_frac_zero']:.4f}")

# Plots: mean vs variance (log-log) and fraction zeros vs mean (log-x)
fig, axs = plt.subplots(1, 2, figsize=(8,4), dpi=300)

# Add a small offset to avoid log(0)
eps = 1e-6
sns.scatterplot(data=metrics.reset_index(), x='mean', y='variance', hue='dataset', ax=axs[0], s=30, alpha=0.6)
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_xlabel('Mean expression (log scale)')
axs[0].set_ylabel('Variance (log scale)')
axs[0].set_title('Mean vs Variance per gene')

sns.scatterplot(data=metrics.reset_index(), x='mean', y='frac_zero', hue='dataset', ax=axs[1], s=30, alpha=0.6)
axs[1].set_xscale('log')
axs[1].set_xlabel('Mean expression (log scale)')
axs[1].set_ylabel('Fraction zeros')
axs[1].set_title('Fraction zeros vs Mean expression per gene')
plt.tight_layout()
plt.savefig(f'{script_dir}/SFig4_panel_D.png')