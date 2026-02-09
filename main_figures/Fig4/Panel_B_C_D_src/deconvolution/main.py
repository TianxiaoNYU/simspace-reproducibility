import numpy as np
import pandas as pd
import scipy

from run_deconvolve import *

import os
script_dir = os.path.dirname(os.path.abspath(__file__))

def rowwise_corr(a, b):
    corrs = []
    for i in range(a.shape[0]):
        if np.std(a[i, :]) == 0 or np.std(b[i, :]) == 0:
            corrs.append(0)
        else:
            corrs.append(np.corrcoef(a[i, :], b[i, :])[0, 1])
    return corrs

def rowwise_smcorr(a, b):
    corrs = []
    for i in range(a.shape[0]):
        if np.std(a[i, :]) == 0 or np.std(b[i, :]) == 0:
            corrs.append(0)
        else:
            corrs.append(scipy.stats.spearmanr(a[i, :], b[i, :]).statistic)
    return corrs

def rowwise_rmse(a, b):
    rmse = []
    for i in range(a.shape[0]):
        if np.std(a[i, :]) == 0 or np.std(b[i, :]) == 0:
            rmse.append(max(np.max(np.sum(a[i, :])), np.sum(np.abs(b[i, :]))))
        else:
            a_row = a[i, :] / np.sum(a[i, :]) if np.sum(a[i, :]) != 0 else a[i, :]
            b_row = b[i, :] / np.sum(b[i, :]) if np.sum(b[i, :]) != 0 else b[i, :]
            tmp = np.sqrt(np.mean(np.square(a_row - b_row)))
            rmse.append(tmp)
    return rmse

def jaccard_index(a, b, threshold=1e-3):
    a = np.array(a)
    b = np.array(b)
    a_bin = np.abs(a) > threshold
    b_bin = np.abs(b) > threshold
    intersect = a_bin & b_bin
    union = a_bin | b_bin
    if np.sum(union) == 0:
        return 0.0  # or np.nan depending on your preference
    return np.sum(intersect) / np.sum(union)

def rowwise_jaccard(a, b):
    jaccards = []
    for i in range(a.shape[0]):
        tmp = jaccard_index(a[i, :], b[i, :])
        jaccards.append(tmp)
    return jaccards

def deconvolve(
    kernel_size = (5, 7, 10, 15),
    n_niche = (2, 3),
    n_state = (8, 9),
    seeds = (0, 1) 
):
    ref_meta_file = f"{script_dir}/../Panel_B_C_D_data/Xenium_reference_metadata.csv"
    ref_omics_file = f"{script_dir}/../Panel_B_C_D_data/Xenium_reference_count.csv"

    for kernel in kernel_size:
        for niches in n_niche:
            for states in n_state:
                for seed in seeds:
                    input_meta_file = f"{script_dir}/tmp_convolve_data/meta/spot_meta_kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv"
                    input_omics_file = f"{script_dir}/tmp_convolve_data/omics/spot_omics_kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv"
                    output_dir = f"{script_dir}/kernel{kernel}_niche{niches}_state{states}_seed{seed}/"
                    
                    run_RCTD(input_meta_file, input_omics_file, ref_meta_file, ref_omics_file, 8, output_dir)
                    run_CARD(input_meta_file, input_omics_file, ref_meta_file, ref_omics_file, output_dir)
                    run_Seurat(input_meta_file, input_omics_file, ref_meta_file, ref_omics_file, output_dir)
                    run_STdeconvolve(input_meta_file, input_omics_file, output_dir)
                    run_spatialDWLS(input_meta_file, input_omics_file, ref_meta_file, ref_omics_file, output_dir)
                    run_cell2location(input_meta_file, input_omics_file, ref_meta_file, ref_omics_file, output_dir)

def compare_results(
    kernel_size = (5, 7, 10, 15),
    n_niche = (2, 3),
    n_state = (8, 9),
    seeds = (0, 1),
    ):
    for kernel in kernel_size:
        for niches in n_niche:
            for states in n_state:
                for seed in seeds:
                    ground_truth_file = f"{script_dir}/tmp_convolve_data/meta/spot_meta_kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv"
                    ground_truth = pd.read_csv(ground_truth_file)
                    ground_truth = ground_truth.iloc[:, 2:]

                    output_dir = f"{script_dir}/kernel{kernel}_niche{niches}_state{states}_seed{seed}/"
                    save_dir = f"{script_dir}/results/summary/"

                    c2l_res = pd.read_csv(output_dir + "c2l_res.csv", index_col=0)
                    RCTD_full_res = pd.read_csv(output_dir + "RCTD_full_res.csv")
                    RCTD_multi_res = pd.read_csv(output_dir + "RCTD_multi_res.csv")
                    CARD_res = pd.read_csv(output_dir + "CARD_res.csv")
                    Seurat_res = pd.read_csv(output_dir + "Seurat_res.csv")
                    # STdeconvolve_res = pd.read_csv(output_dir + "STdeconvolve_res.csv")
                    spatialDWLS_res = pd.read_csv(output_dir + "spatialDWLS_res.csv")
                    ## Filter only the first 'states' columns of each result, ignoring the last two columns
                    
                    if states == 8:                 
                        for df in [c2l_res, RCTD_full_res, RCTD_multi_res, CARD_res, Seurat_res, spatialDWLS_res]:
                            cols_to_drop = [col for col in df.columns if "CD4" in col]
                            df.drop(columns=cols_to_drop, inplace=True)
                    if states == 9 and niches == 3 and seed == 1:
                        for df in [c2l_res, RCTD_full_res, RCTD_multi_res, CARD_res, Seurat_res, spatialDWLS_res]:
                            cols_to_drop = [col for col in df.columns if "CD4" in col]
                            df.drop(columns=cols_to_drop, inplace=True)

                    c2l_res = c2l_res.iloc[:, :-2]
                    RCTD_full_res = RCTD_full_res.iloc[:, :-2]
                    RCTD_multi_res = RCTD_multi_res.iloc[:, :-2]
                    CARD_res = CARD_res.iloc[:, :-2]
                    Seurat_res = Seurat_res.iloc[:, :-2]
                    # STdeconvolve_res = STdeconvolve_res.iloc[:, :-2]
                    spatialDWLS_res = spatialDWLS_res.iloc[:, :-2]

                    # For example, you can calculate the correlation between the results and the ground truth
                    # Calculate row-wise correlation between ground_truth and each result

                    c2l_corr = rowwise_corr(ground_truth.values, c2l_res.values)
                    RCTD_full_corr = rowwise_corr(ground_truth.values, RCTD_full_res.values)
                    RCTD_multi_corr = rowwise_corr(ground_truth.values, RCTD_multi_res.values)
                    CARD_corr = rowwise_corr(ground_truth.values, CARD_res.values)
                    Seurat_corr = rowwise_corr(ground_truth.values, Seurat_res.values)
                    spatialDWLS_corr = rowwise_corr(ground_truth.values, spatialDWLS_res.values)

                    # Save the results as a DataFrame
                    result_df = pd.DataFrame({
                        "c2l": c2l_corr,
                        "RCTD_full": RCTD_full_corr,
                        "RCTD_multi": RCTD_multi_corr,
                        "CARD": CARD_corr,
                        "Seurat": Seurat_corr,
                        "spatialDWLS": spatialDWLS_corr
                    })
                    result_df.to_csv(save_dir + f"pcc/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv", index=False)

                    # Calculate MSE
                    c2l_mse = rowwise_rmse(ground_truth.values, c2l_res.values)
                    RCTD_full_mse = rowwise_rmse(ground_truth.values, RCTD_full_res.values)
                    RCTD_multi_mse = rowwise_rmse(ground_truth.values, RCTD_multi_res.values)
                    CARD_mse = rowwise_rmse(ground_truth.values, CARD_res.values)
                    Seurat_mse = rowwise_rmse(ground_truth.values, Seurat_res.values)
                    spatialDWLS_mse = rowwise_rmse(ground_truth.values, spatialDWLS_res.values)
                    # Save the results as a DataFrame
                    result_df = pd.DataFrame({
                        "c2l": c2l_mse,
                        "RCTD_full": RCTD_full_mse,
                        "RCTD_multi": RCTD_multi_mse,
                        "CARD": CARD_mse,
                        "Seurat": Seurat_mse,
                        "spatialDWLS": spatialDWLS_mse
                    })
                    result_df.to_csv(save_dir + f"rmse/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv", index=False)

                    # Calculate Jaccard index
                    c2l_jaccard = rowwise_jaccard(ground_truth.values, c2l_res.values)
                    RCTD_full_jaccard = rowwise_jaccard(ground_truth.values, RCTD_full_res.values)
                    RCTD_multi_jaccard = rowwise_jaccard(ground_truth.values, RCTD_multi_res.values)
                    CARD_jaccard = rowwise_jaccard(ground_truth.values, CARD_res.values)
                    Seurat_jaccard = rowwise_jaccard(ground_truth.values, Seurat_res.values)
                    spatialDWLS_jaccard = rowwise_jaccard(ground_truth.values, spatialDWLS_res.values)
                    # Save the results as a DataFrame
                    result_df = pd.DataFrame({
                        "c2l": c2l_jaccard,
                        "RCTD_full": RCTD_full_jaccard,
                        "RCTD_multi": RCTD_multi_jaccard,
                        "CARD": CARD_jaccard,
                        "Seurat": Seurat_jaccard,
                        "spatialDWLS": spatialDWLS_jaccard
                    })
                    result_df.to_csv(save_dir + f"jaccard/kernel{kernel}_niche{niches}_state{states}_seed{seed}.csv", index=False)

if __name__ == "__main__":
    deconvolve()
    deconvolve(
        kernel_size = (5, 7, 10, 15),
        n_niche = (3, ),
        n_state = (11, ),
        seeds = (0, 1, 3, 4),
    )
    deconvolve(
        kernel_size = (5, 7, 10, 15),
        n_niche = (3, ),
        n_state = (14, ),
        seeds = (0, 1, 2, 3),
    )
    
    compare_results()
    compare_results(
        kernel_size = (5, 7, 10, 15),
        n_niche = (3, ),
        n_state = (11, ),
        seeds = (0, 1, 3, 4),
    )
    compare_results(
        kernel_size = (5, 7, 10, 15),
        n_niche = (3, ),
        n_state = (14, ),
        seeds = (0, 1, 2, 3),
    )