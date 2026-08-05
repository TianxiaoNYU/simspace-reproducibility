import scanpy as sc
import pandas as pd
from sklearn import metrics
import torch

import matplotlib.pyplot as plt
import seaborn as sns

import os
script_dir = os.path.dirname(os.path.abspath(__file__))
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import SEDR
random_seed = 2023
SEDR.fix_seed(random_seed)

# gpu
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'

Xenium_count = pd.read_csv(f'{script_dir}/../../Panel_B_C_D_data/Xenium_reference_count.csv', index_col=0)
Xenium_meta = pd.read_csv(f'{script_dir}/../../Panel_B_C_D_data/Xenium_reference_metadata.csv', index_col=0)
cell_count = Xenium_meta['Cluster'].value_counts()
Xenium_meta['state_rank'] = 0
for i, cell_type in enumerate(cell_count.index):
    Xenium_meta.loc[Xenium_meta['Cluster'] == cell_type, 'state_rank'] = i + 1
Xenium_meta['x_centroid'] = 100 * (Xenium_meta['x_centroid'] - Xenium_meta['x_centroid'].min()) / (Xenium_meta['x_centroid'].max() - Xenium_meta['x_centroid'].min())
Xenium_meta['y_centroid'] = 100 * (Xenium_meta['y_centroid'] - Xenium_meta['y_centroid'].min()) / (Xenium_meta['y_centroid'].max() - Xenium_meta['y_centroid'].min())

simspace_count = pd.read_csv(f'{script_dir}/../../Panel_B_C_D_data/simspace_fitted_count.csv', index_col=0)
simspace_meta = pd.read_csv(f'{script_dir}/../../Panel_B_C_D_data/simspace_fitted_metadata.csv', index_col=0)
simspace_count = simspace_count.T
simspace_meta['row'] = 100 * (simspace_meta['row'] - simspace_meta['row'].min()) / (simspace_meta['row'].max() - simspace_meta['row'].min())
simspace_meta['col'] = 100 * (simspace_meta['col'] - simspace_meta['col'].min()) / (simspace_meta['col'].max() - simspace_meta['col'].min())
simspace_meta['row'] = simspace_meta['row'] + 200

sccube_count = pd.read_csv(f'{script_dir}/../../Panel_B_C_D_data/scCube_fitted_count.csv', index_col=0)
sccube_meta = pd.read_csv(f'{script_dir}/../../Panel_B_C_D_data/scCube_fitted_metadata.csv', index_col=0)
sccube_meta['point_x'] = 100 * (sccube_meta['point_x'] - sccube_meta['point_x'].min()) / (sccube_meta['point_x'].max() - sccube_meta['point_x'].min())
sccube_meta['point_y'] = 100 * (sccube_meta['point_y'] - sccube_meta['point_y'].min()) / (sccube_meta['point_y'].max() - sccube_meta['point_y'].min())
sccube_meta['point_x'] = sccube_meta['point_x'] + 100

concat_count = pd.concat([Xenium_count, simspace_count, sccube_count], axis=1)
Xenium_meta = Xenium_meta[['x_centroid', 'y_centroid', 'Cluster']]
simspace_meta = simspace_meta[['row', 'col', 'fitted_celltype']]
sccube_meta = sccube_meta[['point_x', 'point_y', 'Cell_type']]
Xenium_meta.columns = ['x_centroid', 'y_centroid', 'Cluster']
simspace_meta.columns = ['x_centroid', 'y_centroid', 'Cluster']
sccube_meta.columns = ['x_centroid', 'y_centroid', 'Cluster']
Xenium_meta['dataset'] = 'Xenium'
simspace_meta['dataset'] = 'SimSpace'
sccube_meta['dataset'] = 'scCube'
concat_meta = pd.concat([Xenium_meta, simspace_meta, sccube_meta], axis=0)
concat_meta.index = concat_count.columns
adata = sc.AnnData(X=concat_count.T, obs=concat_meta)

adata.layers['count'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.scale(adata)
adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values
graph_dict = SEDR.graph_construction(adata, 12)
from sklearn.decomposition import PCA  # sklearn PCA is used because PCA in scanpy is not stable.
adata_X = PCA(n_components=200, random_state=42).fit_transform(adata.X)
adata.obsm['X_pca'] = adata_X
sedr_net = SEDR.Sedr(adata.obsm['X_pca'], graph_dict, mode='clustering', device=device)
using_dec = True
if using_dec:
    sedr_net.train_with_dec(N=1)
else:
    sedr_net.train_without_dec(N=1)
sedr_feat, _, _, _ = sedr_net.process()
adata.obsm['SEDR'] = sedr_feat
import umap

reducer = umap.UMAP(random_state=random_seed)
umap_embedding = reducer.fit_transform(adata.obsm['SEDR'])
adata.obsm['SEDR_umap'] = umap_embedding
umap_df = pd.DataFrame(umap_embedding, index=concat_meta.index, columns=['UMAP1', 'UMAP2'])
umap_df['Dataset'] = concat_meta['dataset']
umap_df['Cluster'] = concat_meta['Cluster']
umap_df.to_csv(f'{script_dir}/umap_embedding.csv')