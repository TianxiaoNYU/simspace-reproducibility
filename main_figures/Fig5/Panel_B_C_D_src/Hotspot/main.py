import hotspot
import anndata
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

load_path = os.path.join(script_dir, '../benchmark_datasets/')
save_path = os.path.join(script_dir)
niche_list = [2,3]
state_list = [8,9]
seed_list = [0,1]

for niche in niche_list:
    for state in state_list:
        for seed in seed_list:
            print(f'Processing niche {niche}, state {state}, seed {seed}')
            meta = pd.read_csv(f'{load_path}/sim_meta_niche{niche}_state{state}_seed{seed}.csv')
            counts = pd.read_csv(f'{load_path}/sim_omics_niche{niche}_state{state}_seed{seed}.csv', index_col=0)

            adata = anndata.AnnData(X=counts.values, obs=meta, var=pd.DataFrame(index=counts.columns))
            adata.layers['counts'] = adata.X.copy()
            adata.obsm['spatial'] = adata.obs[['row', 'col']].values
            adata.obs['total_counts'] = adata.X.sum(axis=1)

            hs = hotspot.Hotspot(
                adata,
                layer_key="counts",
                model='danb',
                latent_obsm_key="spatial",
                umi_counts_obs_key="total_counts"
            )
            hs.create_knn_graph(
                weighted_graph=False, n_neighbors=30,
            )
            hs_results = hs.compute_autocorrelations(jobs=4)
            hs_results.to_csv(f'{save_path}/svg_hotspot_{niche}{state}{seed}.csv')


meta = pd.read_csv(f'{load_path}/Xenium_ref/tile_23_meta.csv', index_col=0)
counts = pd.read_csv(f'{load_path}/Xenium_ref/tile_23_count.csv', index_col=0)
counts = counts.T
adata = anndata.AnnData(X=counts.values, obs=meta, var=pd.DataFrame(index=counts.columns))
adata.layers['counts'] = adata.X.copy()
adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values
adata.obs['total_counts'] = adata.X.sum(axis=1)

hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model='danb',
    latent_obsm_key="spatial",
    umi_counts_obs_key="total_counts"
)
hs.create_knn_graph(
    weighted_graph=False, n_neighbors=30,
)
hs_results = hs.compute_autocorrelations(jobs=4)
hs_results.to_csv(f'{save_path}/svg_hotspot_reference.csv')

meta = pd.read_csv(f'{load_path}/Xenium_fitted/tile_23_fitted_meta.csv', index_col=0)
counts = pd.read_csv(f'{load_path}/Xenium_fitted/tile_23_fitted_count.csv', index_col=0)
adata = anndata.AnnData(X=counts.values, obs=meta, var=pd.DataFrame(index=counts.columns))
adata.layers['counts'] = adata.X.copy()
adata.obsm['spatial'] = adata.obs[['row', 'col']].values
adata.obs['total_counts'] = adata.X.sum(axis=1)

hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model='danb',
    latent_obsm_key="spatial",
    umi_counts_obs_key="total_counts"
)
hs.create_knn_graph(
    weighted_graph=False, n_neighbors=30,
)
hs_results = hs.compute_autocorrelations(jobs=4)
hs_results.to_csv(f'{save_path}/svg_hotspot_fitted.csv')