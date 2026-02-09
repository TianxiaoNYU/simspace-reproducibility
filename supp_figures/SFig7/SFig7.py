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
from scipy.stats import linregress

########## Panel B:
BSP_dir = f"{script_dir}/data/scBSP"
hotspot_dir = f"{script_dir}/data/hotspot"
giotto_rank_dir = f"{script_dir}/data/giotto_rank"
giotto_kmeans_dir = f"{script_dir}/data/giotto_kmeans"
MI_dir = f"{script_dir}/data/moransi"
SpatialDE_dir = f"{script_dir}/data/spatialde"
spark_dir = f"{script_dir}/data/spark"

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

########## Panel scBSP:
slope, intercept, r_value, p_value, std_err = linregress((combined_ranks['BSP_fitted']), (combined_ranks['BSP_reference']))

plt.figure(figsize=(4, 4), dpi=250)
plt.scatter(combined_ranks['BSP_fitted'], combined_ranks['BSP_reference'], alpha=0.7)
plt.plot([combined_ranks['BSP_fitted'].min(), combined_ranks['BSP_fitted'].max()], 
         [slope * combined_ranks['BSP_reference'].min() + intercept, slope * combined_ranks['BSP_reference'].max() + intercept],
         color='red', linestyle='--', label=f"r: {r_value:.4f}")
plt.legend()
plt.xlabel('SimSpace scBSP Ranks')
plt.ylabel('Reference scBSP Ranks')
plt.title('scBSP Simulated vs Reference Ranks')
plt.savefig(f"{script_dir}/SFig7_scBSP.png", bbox_inches='tight')

########## Panel SPARK:
slope, intercept, r_value, p_value, std_err = linregress((combined_ranks['Spark_fitted']), (combined_ranks['Spark_reference']))

plt.figure(figsize=(4, 4), dpi=250)
plt.scatter(combined_ranks['Spark_fitted'], combined_ranks['Spark_reference'], alpha=0.7)
plt.plot([combined_ranks['Spark_fitted'].min(), combined_ranks['Spark_fitted'].max()], 
         [slope * combined_ranks['Spark_reference'].min() + intercept, slope * combined_ranks['Spark_reference'].max() + intercept],
         color='red', linestyle='--', label=f"r: {r_value:.4f}")
plt.legend()
plt.xlabel('SimSpace SPARK Ranks')
plt.ylabel('Reference SPARK Ranks')
plt.title('SPARK Simulated vs Reference Ranks')
plt.savefig(f"{script_dir}/SFig7_SPARK.png", bbox_inches='tight')

########### Panel Hotspot:
slope, intercept, r_value, p_value, std_err = linregress((combined_ranks['Hotspot_fitted']), (combined_ranks['Hotspot_reference']))

plt.figure(figsize=(4, 4), dpi=250)
plt.scatter(combined_ranks['Hotspot_fitted'], combined_ranks['Hotspot_reference'], alpha=0.7)
plt.plot([combined_ranks['Hotspot_fitted'].min(), combined_ranks['Hotspot_fitted'].max()], 
         [slope * combined_ranks['Hotspot_reference'].min() + intercept, slope * combined_ranks['Hotspot_reference'].max() + intercept],
         color='red', linestyle='--', label=f"r: {r_value:.4f}")
plt.legend()
plt.xlabel('SimSpace Hotspot Ranks')
plt.ylabel('Reference Hotspot Ranks')
plt.title('Hotspot Simulated vs Reference Ranks')
plt.savefig(f"{script_dir}/SFig7_Hotspot.png", bbox_inches='tight')

########### Panel Morans Index:
slope, intercept, r_value, p_value, std_err = linregress((combined_ranks['MI_fitted']), (combined_ranks['MI_reference']))

plt.figure(figsize=(4, 4), dpi=250)
plt.scatter(combined_ranks['MI_fitted'], combined_ranks['MI_reference'], alpha=0.7)
plt.plot([combined_ranks['MI_fitted'].min(), combined_ranks['MI_fitted'].max()], 
         [slope * combined_ranks['MI_reference'].min() + intercept, slope * combined_ranks['MI_reference'].max() + intercept],
         color='red', linestyle='--', label=f"r: {r_value:.4f}")
plt.legend()
plt.xlabel('SimSpace MI Ranks')
plt.ylabel('Reference MI Ranks')
plt.title('MI Simulated vs Reference Ranks')
plt.savefig(f"{script_dir}/SFig7_MI.png", bbox_inches='tight')

########### Panel Giotto kmeans:
slope, intercept, r_value, p_value, std_err = linregress((combined_ranks['Giotto_kmeans_fitted']), (combined_ranks['Giotto_kmeans_reference']))

plt.figure(figsize=(4, 4), dpi=250)
plt.scatter(combined_ranks['Giotto_kmeans_fitted'], combined_ranks['Giotto_kmeans_reference'], alpha=0.7)
plt.plot([combined_ranks['Giotto_kmeans_fitted'].min(), combined_ranks['Giotto_kmeans_fitted'].max()], 
         [slope * combined_ranks['Giotto_kmeans_reference'].min() + intercept, slope * combined_ranks['Giotto_kmeans_reference'].max() + intercept],
         color='red', linestyle='--', label=f"r: {r_value:.4f}")
plt.legend()
plt.xlabel('SimSpace Giotto Ranks')
plt.ylabel('Reference Giotto Ranks')
plt.title('Giotto (kmeans) Simulated vs Reference Ranks')
plt.savefig(f"{script_dir}/SFig7_Giotto_kmeans.png", bbox_inches='tight')

########### Panel Giotto rank:
slope, intercept, r_value, p_value, std_err = linregress((combined_ranks['Giotto_rank_fitted']), (combined_ranks['Giotto_rank_reference']))

plt.figure(figsize=(4, 4), dpi=250)
plt.scatter(combined_ranks['Giotto_rank_fitted'], combined_ranks['Giotto_rank_reference'], alpha=0.7)
plt.plot([combined_ranks['Giotto_rank_fitted'].min(), combined_ranks['Giotto_rank_fitted'].max()], 
         [slope * combined_ranks['Giotto_rank_reference'].min() + intercept, slope * combined_ranks['Giotto_rank_reference'].max() + intercept],
         color='red', linestyle='--', label=f"r: {r_value:.4f}")
plt.legend()
plt.xlabel('SimSpace Giotto Ranks')
plt.ylabel('Reference Giotto Ranks')
plt.title('Giotto (rank) Simulated vs Reference Ranks')
plt.savefig(f"{script_dir}/SFig7_Giotto_rank.png", bbox_inches='tight')

########### Panel SpatialDE:
slope, intercept, r_value, p_value, std_err = linregress((combined_ranks['SpatialDE_fitted']), (combined_ranks['SpatialDE_reference']))

plt.figure(figsize=(4, 4), dpi=250)
plt.scatter(combined_ranks['SpatialDE_fitted'], combined_ranks['SpatialDE_reference'], alpha=0.7)
plt.plot([combined_ranks['SpatialDE_fitted'].min(), combined_ranks['SpatialDE_fitted'].max()], 
         [slope * combined_ranks['SpatialDE_reference'].min() + intercept, slope * combined_ranks['SpatialDE_reference'].max() + intercept],
         color='red', linestyle='--', label=f"r: {r_value:.4f}")
plt.legend()
plt.xlabel('SimSpace SpatialDE Ranks')
plt.ylabel('Reference SpatialDE Ranks')
plt.title('SpatialDE Simulated vs Reference Ranks')
plt.savefig(f"{script_dir}/SFig7_SpatialDE.png", bbox_inches='tight')