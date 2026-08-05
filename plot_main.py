import warnings
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc
import os

import simspace as ss
import subprocess
import glob
# Find the exact entry script for each main-figure directory.
figure_dirs = glob.glob('main_figures/Fig*')
fig_scripts = [
    os.path.join(figure_dir, f'{os.path.basename(figure_dir)}.py')
    for figure_dir in figure_dirs
    if os.path.isfile(os.path.join(figure_dir, f'{os.path.basename(figure_dir)}.py'))
]

# Sort by manuscript figure number.
fig_scripts.sort(key=lambda path: int(os.path.basename(os.path.dirname(path))[3:]))

print(f"Found {len(fig_scripts)} figure scripts to run:")
for script in fig_scripts:
    print(f"  - {script}")

# Run each script
for script in fig_scripts:
    print(f"\nRunning {script}...")
    try:
        result = subprocess.run(['python', script], 
                              capture_output=True, 
                              text=True, 
                              check=True)
        print(f"✓ {script} completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"✗ Error running {script}:")
        print(f"  Return code: {e.returncode}")
        print(f"  Error output: {e.stderr}")
    except Exception as e:
        print(f"✗ Unexpected error running {script}: {e}")

print("\nAll figure scripts have been processed.")
