import scanpy as sc
import pandas as pd
import numpy as np
import os 
script_dir = os.path.dirname(os.path.abspath(__file__))
import pickle

import cell2location
import argparse

import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='cell2location deconvolution',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_meta_file', type=str, default=script_dir + '/spot_meta_kernel5_niche2_state8_seed0.csv',
                    help='input meta file')
parser.add_argument('--input_omics_file', type=str, default=script_dir + '/spot_omics_kernel5_niche2_state8_seed0.csv',
                    help='input omics file')
parser.add_argument('--ref_meta_file', type=str, default=script_dir + '/tile_23_meta_update.csv',
                    help='reference meta file')
parser.add_argument('--ref_omics_file', type=str, default=script_dir + '/tile_23_count.csv',
                    help='reference omics file')
parser.add_argument('--output_dir', type=str, default=script_dir,
                    help='output directory')
arg = parser.parse_args()

def main():
    count_data = pd.read_csv(arg.input_omics_file)
    meta_data = pd.read_csv(arg.input_meta_file)

    reference_data = pd.read_csv(arg.ref_omics_file, index_col=0)
    reference_data = reference_data.T
    reference_meta = pd.read_csv(arg.ref_meta_file, index_col=0)
    reference_meta['Batch'] = np.array(["A"] * reference_meta.shape[0])
    reference_adata = sc.AnnData(reference_data)
    reference_adata.obs = reference_meta

    adata = sc.AnnData(count_data)
    adata.obs = meta_data

    cell2location.models.RegressionModel.setup_anndata(
        adata=reference_adata,
        batch_key='Batch',
        labels_key='Cluster',
        )
    from cell2location.models import RegressionModel
    mod = RegressionModel(reference_adata)
    mod.train(max_epochs=3000)
    reference_adata = mod.export_posterior(
        reference_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 256}
    )
    if 'means_per_cluster_mu_fg' in reference_adata.varm.keys():
        inf_aver = reference_adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in reference_adata.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = reference_adata.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in reference_adata.uns['mod']['factor_names']]].copy()
    inf_aver.columns = reference_adata.uns['mod']['factor_names']

    inf_aver_path = os.path.join(script_dir, "inf_aver.pkl")
    mod_path = os.path.join(script_dir, "regression_model")
    
    # inf_aver_path = os.path.join("/Users/zhaotianxiao/Library/CloudStorage/Dropbox/FenyoLab/Project/Spatialsim/output/deconvolution/results/raw/inf_aver.pkl")
    # mod_path = os.path.join("/Users/zhaotianxiao/Library/CloudStorage/Dropbox/FenyoLab/Project/Spatialsim/output/deconvolution/results/raw/regression_model")
    
    # Try to load inf_aver and mod if they exist
    if os.path.exists(inf_aver_path) and os.path.exists(mod_path):
        with open(inf_aver_path, "rb") as f:
            inf_aver = pickle.load(f)
        mod = RegressionModel(reference_adata)
        mod.load(mod_path)
        # reference_adata = mod.export_posterior(
        #     reference_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 256}
        # )
    else:
        mod = RegressionModel(reference_adata)
        mod.train(max_epochs=3000)
        reference_adata = mod.export_posterior(
            reference_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 256}
        )
        if 'means_per_cluster_mu_fg' in reference_adata.varm.keys():
            inf_aver = reference_adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                            for i in reference_adata.uns['mod']['factor_names']]].copy()
        else:
            inf_aver = reference_adata.var[[f'means_per_cluster_mu_fg_{i}'
                                            for i in reference_adata.uns['mod']['factor_names']]].copy()
        inf_aver.columns = reference_adata.uns['mod']['factor_names']
        # Save inf_aver and mod
        with open(inf_aver_path, "wb") as f:
            pickle.dump(inf_aver, f)
        mod.save(mod_path, save_anndata=True, overwrite=True)



    cell2location.models.Cell2location.setup_anndata(adata=adata)
    mod = cell2location.models.Cell2location(
        adata, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=30,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20
    )
    mod.train(max_epochs=5000,
            # train using full data (batch_size=None)
            batch_size=None,
            # use all data points in training because
            # we need to estimate cell abundance at all locations
            train_size=1,
            )
    adata = mod.export_posterior(
        adata, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, }
    )
    # adata.obsm['means_cell_abundance_w_sf'].to_csv(f'{arg.output_dir}/c2l_res.csv')
    # Remove 'means_cell_abundance_w_sf_' prefix from column names before saving
    df = adata.obsm['q05_cell_abundance_w_sf'].copy()
    df.columns = [col.replace('q05cell_abundance_w_sf_', '') for col in df.columns]

    df = df[sorted(df.columns)]
    df['X'] = adata.obs['col']
    df['Y'] = adata.obs['row']

    df.to_csv(f'{arg.output_dir}/c2l_res.csv')


if __name__ == "__main__":
    main()




