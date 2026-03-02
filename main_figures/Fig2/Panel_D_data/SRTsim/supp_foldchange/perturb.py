import pandas as pd
import numpy as np
from pathlib import Path

import os

script_dir_global = os.path.dirname(os.path.abspath(__file__))
script_dir = os.path.join(script_dir_global, '../../../')

base_dir = os.path.join(script_dir, "Panel_D_data/SRTsim")
output_dir_1 = Path(f"{base_dir}/supp_foldchange/perturbed_1")
output_dir_1.mkdir(parents=True, exist_ok=True)
output_dir_2 = Path(f"{base_dir}/supp_foldchange/perturbed_2")
output_dir_2.mkdir(parents=True, exist_ok=True)
# Set global seed to 42
np.random.seed(42)

# Process ref_free_0.csv to ref_free_8.csv
for i in range(9):
    input_file = f"{base_dir}/ref_free_{i}.csv"
    df = pd.read_csv(input_file)
    
    # Generate random fold changes from U(0.5, 2.5) for each unique group
    unique_groups = df['group'].unique()
    group_foldchange_map = {group: np.round(np.random.uniform(0.5, 2.5), 3) for group in unique_groups}
    
    # Apply the fold changes to each row based on its group
    df['foldchange'] = df['group'].map(group_foldchange_map)
    
    # Save the modified CSV
    output_file = output_dir_1 / f"ref_free_{i}.csv"
    df.to_csv(output_file, index=False)

    # Generate another set of random fold changes for the second perturbed dataset
    group_foldchange_map_2 = {group: np.round(np.random.uniform(0.5, 2.5), 3) for group in unique_groups}
    df['foldchange'] = df['group'].map(group_foldchange_map_2)
    output_file_2 = output_dir_2 / f"ref_free_{i}.csv"
    df.to_csv(output_file_2, index=False)

    print(f"Perturbed CSV saved to {output_file}")
    print(f"Perturbed CSV saved to {output_file_2}")
    # print(f"Group fold changes: {group_foldchange_map}")
    # print(f"Group fold changes 2: {group_foldchange_map_2}")

