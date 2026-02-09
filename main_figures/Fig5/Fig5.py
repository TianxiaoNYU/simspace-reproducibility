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

########## Panel B:
BSP_dir = f"{script_dir}/Panel_B_C_D_data/scBSP"
hotspot_dir = f"{script_dir}/Panel_B_C_D_data/hotspot"
giotto_rank_dir = f"{script_dir}/Panel_B_C_D_data/giotto_rank"
giotto_kmeans_dir = f"{script_dir}/Panel_B_C_D_data/giotto_kmeans"
MI_dir = f"{script_dir}/Panel_B_C_D_data/moransi"
SpatialDE_dir = f"{script_dir}/Panel_B_C_D_data/spatialde"
spark_dir = f"{script_dir}/Panel_B_C_D_data/spark"

BSP_fitted = pd.read_csv(f"{BSP_dir}/svg_scBSP_fitted.csv")
Hotspot_fitted = pd.read_csv(f"{hotspot_dir}/svg_hotspot_fitted.csv")
Giotto_fitted_1 = pd.read_csv(f"{giotto_kmeans_dir}/svg_giotto_kmeans_xenium.csv")
Giotto_fitted_2 = pd.read_csv(f"{giotto_rank_dir}/svg_giotto_rank_xenium.csv")
MI_fitted = pd.read_csv(f"{MI_dir}/svg_moransi_xenium.csv")
SpatialDE_fitted = pd.read_csv(f"{SpatialDE_dir}/svg_spatialde_xenium.csv")
Spark_fitted = pd.read_csv(f"{spark_dir}/svg_spark_fitted.csv", index_col=0)

BSP_ref = pd.read_csv(f"{BSP_dir}/svg_scBSP_reference.csv")
Hotspot_ref = pd.read_csv(f"{hotspot_dir}/svg_hotspot_reference.csv")
Giotto_ref_1 = pd.read_csv(f"{giotto_kmeans_dir}/svg_giotto_kmeans_ref.csv")
Giotto_ref_2 = pd.read_csv(f"{giotto_rank_dir}/svg_giotto_rank_ref.csv")
MI_ref = pd.read_csv(f"{MI_dir}/svg_moransi_ref.csv")
SpatialDE_ref = pd.read_csv(f"{SpatialDE_dir}/svg_spatialde_ref.csv")
Spark_ref = pd.read_csv(f"{spark_dir}/svg_spark_reference.csv", index_col=0)

Hotspot_fitted = Hotspot_fitted.sort_values(by='Gene')
Hotspot_fitted.reset_index(drop=True, inplace=True)

Hotspot_ref = Hotspot_ref.sort_values(by='Gene')
Hotspot_ref.reset_index(drop=True, inplace=True)

Giotto_fitted_1 = Giotto_fitted_1.sort_values(by='genes')
Giotto_fitted_1.reset_index(drop=True, inplace=True)

Giotto_ref_1 = Giotto_ref_1.sort_values(by='genes')
Giotto_ref_1.reset_index(drop=True, inplace=True)

Giotto_fitted_2 = Giotto_fitted_2.sort_values(by='genes')
Giotto_fitted_2.reset_index(drop=True, inplace=True)

Giotto_ref_2 = Giotto_ref_2.sort_values(by='genes')
Giotto_ref_2.reset_index(drop=True, inplace=True)

MI_fitted = MI_fitted.sort_values(by='genes')
MI_fitted.reset_index(drop=True, inplace=True)
MI_ref = MI_ref.sort_values(by='genes')
MI_ref.reset_index(drop=True, inplace=True)

Spark_fitted.index = Hotspot_fitted.index
Spark_ref.index = Hotspot_ref.index


combined = pd.DataFrame({
    'Gene': Hotspot_fitted['Gene'],
    'Hotspot_fitted': Hotspot_fitted['FDR'],
    'Hotspot_reference': Hotspot_ref['FDR'],
    'BSP_fitted': BSP_fitted['P_values'],
    'BSP_reference': BSP_ref['P_values'],
    'Giotto_kmeans_fitted': Giotto_fitted_1['adj.p.value'],
    'Giotto_kmeans_reference': Giotto_ref_1['adj.p.value'],
    'Giotto_rank_fitted': Giotto_fitted_2['adj.p.value'],
    'Giotto_rank_reference': Giotto_ref_2['adj.p.value'],
    'MI_fitted': MI_fitted['p.value'],
    'MI_reference': MI_ref['p.value'],
    'SpatialDE_fitted': SpatialDE_fitted['qval'],
    'SpatialDE_reference': SpatialDE_ref['qval'],
    'Spark_fitted': Spark_fitted['adjusted_pvalue'],
    'Spark_reference': Spark_ref['adjusted_pvalue']
})

# Convert p-values to ranks for each column in combined (smaller p-value = higher rank)
combined_ranks = combined.copy()
for col in combined.columns[1:]:  # skip 'Gene'
    combined_ranks[col] = combined[col].rank(method='min', ascending=True)

# Prepare fitted p-value data for heatmap
fitted_pvals = pd.DataFrame({
    'BSP': BSP_fitted['P_values'].values,
    'Hotspot': Hotspot_fitted['FDR'].values,
    'Giotto_kmeans': Giotto_fitted_1['adj.p.value'].values,
    'Giotto_rank': Giotto_fitted_2['adj.p.value'].values,
    "Moran's I": MI_fitted['p.value'].values,
    'SpatialDE': SpatialDE_fitted['qval'].values,
    'SPARK': Spark_fitted['adjusted_pvalue'].values
}, index=BSP_fitted['GeneNames'])

fitted_pvals.to_csv(f'{script_dir}/Panel_B_C_D_data/fitted_pvals.csv')
# Run the R script for panel B
result = subprocess.run(['Rscript', f'{script_dir}/Panel_B_C_D_src/panel_B.R'], 
                       capture_output=True, text=True)
if result.returncode != 0:
    print(f"Error running R script: {result.stderr}")
else:
    print("R script executed successfully")

########## Panel C:
from scipy import stats
from scipy.stats import linregress
# from scipy.stats import fisher_exact
from sklearn.metrics import f1_score

# Compute R^2 for each method
r_results = {}

# BSP
slope, intercept, r_value, p_value, std_err = linregress(
    combined_ranks['BSP_fitted'], combined_ranks['BSP_reference']
)
r_results['scBSP'] = r_value
# Hotspot (nonzero only)
slope, intercept, r_value, p_value, std_err = linregress(
    combined_ranks['Hotspot_fitted'], combined_ranks['Hotspot_reference']
)
r_results['Hotspot'] = r_value
# Giotto kmeans
slope, intercept, r_value, p_value, std_err = linregress(
    combined_ranks['Giotto_kmeans_fitted'],
    combined_ranks['Giotto_kmeans_reference']
)
r_results['Giotto(kmeans)'] = r_value
# Giotto rank
slope, intercept, r_value, p_value, std_err = linregress(
    combined_ranks['Giotto_rank_fitted'],
    combined_ranks['Giotto_rank_reference']
)
r_results['Giotto(rank)'] = r_value
# MI
slope, intercept, r_value, p_value, std_err = linregress(
    combined_ranks['MI_fitted'],
    combined_ranks['MI_reference']
)
r_results['Moran\'s I'] = r_value
# SPARK
slope, intercept, r_value, p_value, std_err = linregress(
    combined_ranks['Spark_fitted'],
    combined_ranks['Spark_reference']
)
r_results['SPARK'] = r_value
# SpatialDE (nonzero only)
slope, intercept, r_value, p_value, std_err = linregress(
    combined_ranks['SpatialDE_fitted'],
    combined_ranks['SpatialDE_reference']
)
r_results['SpatialDE'] = r_value

# Summarize
r_summary = pd.Series(r_results, name='r')

# Define threshold
threshold = 5e-2
# Prepare results dict
f1_ratios = {}

# Helper function to compute odds ratio
def compute_f1(fitted, reference):
    # True = significant (< threshold), False = not significant
    fitted_sig = fitted < threshold
    reference_sig = reference < threshold
    # Compute F1 score (treat reference as ground truth)
    return f1_score(reference_sig, fitted_sig)

# BSP
f1_ratios['scBSP'] = compute_f1(combined['BSP_fitted'], combined['BSP_reference'])
# Hotspot
f1_ratios['Hotspot'] = compute_f1(combined['Hotspot_fitted'], combined['Hotspot_reference'])
# Giotto kmeans
f1_ratios['Giotto(kmeans)'] = compute_f1(combined['Giotto_kmeans_fitted'], combined['Giotto_kmeans_reference'])
# Giotto rank
f1_ratios['Giotto(rank)'] = compute_f1(combined['Giotto_rank_fitted'], combined['Giotto_rank_reference'])
# MI
f1_ratios['Moran\'s I'] = compute_f1(combined['MI_fitted'], combined['MI_reference'])
# SpatialDE
f1_ratios['SpatialDE'] = compute_f1(combined['SpatialDE_fitted'], combined['SpatialDE_reference'])
# SPARK
f1_ratios['SPARK'] = compute_f1(combined['Spark_fitted'], combined['Spark_reference'])

# Summarize
f1_summary = pd.Series(f1_ratios, name='F1_score')


# Prepare results dict
kendall_taus = {}
def compute_kendall_tau(fitted, reference):
    # Compute Kendall's tau (treat reference as ground truth)
    return stats.kendalltau(reference, fitted).correlation

# BSP
kendall_taus['scBSP'] = compute_kendall_tau(combined_ranks['BSP_fitted'], combined_ranks['BSP_reference'])
# Hotspot
kendall_taus['Hotspot'] = compute_kendall_tau(combined_ranks['Hotspot_fitted'], combined_ranks['Hotspot_reference'])
# Giotto kmeans
kendall_taus['Giotto(kmeans)'] = compute_kendall_tau(combined_ranks['Giotto_kmeans_fitted'], combined_ranks['Giotto_kmeans_reference'])
# Giotto rank
kendall_taus['Giotto(rank)'] = compute_kendall_tau(combined_ranks['Giotto_rank_fitted'], combined_ranks['Giotto_rank_reference'])
# MI
kendall_taus['Moran\'s I'] = compute_kendall_tau(combined_ranks['MI_fitted'], combined_ranks['MI_reference'])
# SpatialDE
kendall_taus['SpatialDE'] = compute_kendall_tau(combined_ranks['SpatialDE_fitted'], combined_ranks['SpatialDE_reference'])
# SPARK
kendall_taus['SPARK'] = compute_kendall_tau(combined_ranks['Spark_fitted'], combined_ranks['Spark_reference'])

# Summarize
kendall_tau_summary = pd.Series(kendall_taus, name='Kendall\'s Tau')

# Sort r_summary by value in descending order
plt.figure(figsize=(3, 4), dpi=250)
r_sorted = r_summary.sort_values(ascending=False)
sns.barplot(x=r_sorted.index, y=r_sorted.values, palette='tab20')
plt.ylabel('Pearson Correlation r')
# plt.title('Pearson r for Each Method')
plt.xticks(rotation=55, ha='right', rotation_mode='anchor')
plt.ylim(0, 1)
plt.xlabel(None)
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig5_panel_C1.png')

plt.figure(figsize=(3, 4), dpi=250)
kendall_tau_sorted = kendall_tau_summary[r_sorted.index] # type: ignore
sns.barplot(x=kendall_tau_sorted.index, y=kendall_tau_sorted.values, palette='tab20')
plt.ylabel('Kendall\'s Tau Coefficient')
# plt.title('F1 Score of Spatial Methods')
plt.xticks(rotation=55, ha='right', rotation_mode='anchor')
plt.xlabel(None)
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig5_panel_C2.png')

plt.figure(figsize=(3, 4), dpi=250)
f1_sorted = f1_summary[r_sorted.index] # type: ignore
sns.barplot(x=f1_sorted.index, y=f1_sorted.values, palette='tab20')
plt.ylabel('F1 Score')
# plt.title('F1 Score of Spatial Methods')
plt.xticks(rotation=55, ha='right', rotation_mode='anchor')
plt.ylim(0, 1)
plt.xlabel(None)
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig5_panel_C3.png')


########## Panel D:
method_list = [
    'scBSP',
    'hotspot',
    'giotto_rank',
    'giotto_kmeans',
    'moransi',
    'spatialde',
    'spark'
]

replicate_index = [280, 281, 290, 291, 380, 381, 390, 391]

from scipy.stats import chi2 as chi2_dist
from scipy.stats import spearmanr
import itertools

def kendalls_w(matrix):
    """
    matrix: rows are items (genes), columns are rankings (from replicates)
    """
    k = matrix.shape[1]  # number of rankings
    n = matrix.shape[0]  # number of items

    # Compute row means
    row_means = np.mean(matrix, axis=1)
    S = np.sum((np.sum(matrix, axis=1) - k * row_means.mean())**2)
    W = 12 * S / (k**2 * (n**3 - n))

    # Significance test (Kendall & Gibbons, 1990)
    # Under H0, the test statistic is: chi2 = k*(n-1)*W, df = n-1
    chi2 = k * (n - 1) * W
    p_value = 1 - chi2_dist.cdf(chi2, df=n-1)

    return W, chi2, p_value

from scipy.stats import spearmanr
import itertools

def avg_spearman(rank_matrix):
    pairs = list(itertools.combinations(range(rank_matrix.shape[1]), 2))
    corrs = [spearmanr(rank_matrix[:, i], rank_matrix[:, j])[0] for i, j in pairs]
    return np.mean(corrs)

all_replicates = pd.DataFrame()
for replicate in replicate_index:
    replicate_file = f"{script_dir}/Panel_B_C_D_data/scBSP/svg_scBSP_{replicate}.csv"
    replicate_data = pd.read_csv(replicate_file, index_col=0)
    replicate_data['replicate'] = replicate
    all_replicates[f'{replicate}'] = replicate_data['P_values']

# Convert each column of all_replicates to ranks (1 = lowest value)
all_replicates_ranked = all_replicates.rank(axis=0, method='average', ascending=True)
scBSP_spearman = avg_spearman(all_replicates_ranked.values)
scBSP_kw = kendalls_w(all_replicates_ranked.values)

all_replicates = pd.DataFrame()
for replicate in replicate_index:
    replicate_file = f"{script_dir}/Panel_B_C_D_data/giotto_kmeans/svg_giotto_kmeans_{replicate}.csv"
    replicate_data = pd.read_csv(replicate_file, index_col=0)
    replicate_data['replicate'] = replicate
    all_replicates[f'{replicate}'] = -np.log(replicate_data['adj.p.value'])

# Convert each column of all_replicates to ranks (1 = lowest value)
all_replicates_ranked = all_replicates.rank(axis=0, method='average', ascending=True)
giotto_kmeans_spearman = avg_spearman(all_replicates_ranked.values)
giotto_kmeans_kw = kendalls_w(all_replicates_ranked.values)

all_replicates = pd.DataFrame()
for replicate in replicate_index:
    replicate_file = f"{script_dir}/Panel_B_C_D_data/giotto_rank/svg_giotto_rank_{replicate}.csv"
    replicate_data = pd.read_csv(replicate_file, index_col=0)
    replicate_data['replicate'] = replicate
    all_replicates[f'{replicate}'] = -np.log(replicate_data['adj.p.value'])

# Convert each column of all_replicates to ranks (1 = lowest value)
all_replicates_ranked = all_replicates.rank(axis=0, method='average', ascending=True)
giotto_rank_spearman = avg_spearman(all_replicates_ranked.values)
giotto_rank_kw = kendalls_w(all_replicates_ranked.values)

all_replicates = pd.DataFrame()
for replicate in replicate_index:
    replicate_file = f"{script_dir}/Panel_B_C_D_data/hotspot/svg_hotspot_{replicate}.csv"
    replicate_data = pd.read_csv(replicate_file, index_col=0)
    replicate_data['replicate'] = replicate
    all_replicates[f'{replicate}'] = -np.log(replicate_data['FDR'])

# Convert each column of all_replicates to ranks (1 = lowest value)
all_replicates_ranked = all_replicates.rank(axis=0, method='average', ascending=True)
hotspot_spearman = avg_spearman(all_replicates_ranked.values)
hotspot_kw = kendalls_w(all_replicates_ranked.values)

all_replicates = pd.DataFrame()
for replicate in replicate_index:
    replicate_file = f"{script_dir}/Panel_B_C_D_data/moransi/svg_moransi_{replicate}.csv"
    replicate_data = pd.read_csv(replicate_file, index_col=0)
    replicate_data['replicate'] = replicate
    all_replicates[f'{replicate}'] = replicate_data['p.value']

# Convert each column of all_replicates to ranks (1 = lowest value)
all_replicates_ranked = all_replicates.rank(axis=0, method='average', ascending=True)
moransi_spearman = avg_spearman(all_replicates_ranked.values)
MI_kw = kendalls_w(all_replicates_ranked.values)

all_replicates = pd.DataFrame()
for replicate in replicate_index:
    replicate_file = f"{script_dir}/Panel_B_C_D_data/spatialde/svg_spatialde_sim{replicate}.csv"
    replicate_data = pd.read_csv(replicate_file, index_col=0)
    replicate_data['replicate'] = replicate
    all_replicates[f'{replicate}'] = replicate_data['pval']

# Convert each column of all_replicates to ranks (1 = lowest value)
all_replicates_ranked = all_replicates.rank(axis=0, method='average', ascending=True)
spatialde_spearman = avg_spearman(all_replicates_ranked.values)
spatialde_kw = kendalls_w(all_replicates_ranked.values)

### spatialde
all_replicates = pd.DataFrame()
for replicate in replicate_index:
    replicate_file = f"{script_dir}/Panel_B_C_D_data/spark/svg_spark_{replicate}.csv"
    replicate_data = pd.read_csv(replicate_file, index_col=0)
    replicate_data['replicate'] = replicate
    all_replicates[f'{replicate}'] = replicate_data['adjusted_pvalue']

# Convert each column of all_replicates to ranks (1 = lowest value)
all_replicates_ranked = all_replicates.rank(axis=0, method='average', ascending=True)
spark_spearman = avg_spearman(all_replicates_ranked.values)
spark_kw = kendalls_w(all_replicates_ranked.values)

kw_summary = pd.Series({
    'scBSP': scBSP_kw[0],
    'Hotspot': hotspot_kw[0],
    'Giotto(kmeans)': giotto_kmeans_kw[0],
    'Giotto(rank)': giotto_rank_kw[0],
    'Moran\'s I': MI_kw[0],
    'SpatialDE': spatialde_kw[0],
    'SPARK': spark_kw[0]
}, name='Kendall\'s W')

plt.figure(figsize=(3, 4), dpi=250)
# kw_sorted = kw_summary.sort_values(ascending=False)
kw_sorted = kw_summary[r_sorted.index]
sns.barplot(x=kw_sorted.index, y=kw_sorted.values, palette='tab20')
plt.ylabel('Kendall\'s W Coefficient')
# plt.title('Pearson r for Each Method')
plt.xticks(rotation=55, ha='right', rotation_mode='anchor')
plt.ylim(0, 1)
plt.xlabel(None)
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig5_panel_D1.png')

spearman_summary = pd.Series({
    'scBSP': scBSP_spearman,
    'Hotspot': hotspot_spearman,
    'Giotto(kmeans)': giotto_kmeans_spearman,
    'Giotto(rank)': giotto_rank_spearman,
    'Moran\'s I': moransi_spearman,
    'SpatialDE': spatialde_spearman,
    'SPARK': spark_spearman,
}, name='Spearman Correlation')

plt.figure(figsize=(3, 4), dpi=250)
spearman_sorted = spearman_summary[r_sorted.index]
sns.barplot(x=spearman_sorted.index, y=spearman_sorted.values, palette='tab20')
plt.ylabel('Average Spearman Correlation')
# plt.title('Pearson r for Each Method')
plt.xticks(rotation=55, ha='right', rotation_mode='anchor')
plt.ylim(0, 1)
plt.xlabel(None)
plt.tight_layout()
plt.savefig(f'{script_dir}/Fig5_panel_D2.png')