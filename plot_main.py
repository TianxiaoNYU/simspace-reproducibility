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
# Find all FigX.py scripts in the main_figures directory
fig_scripts = glob.glob('main_figures/Fig*/Fig*.py')

# Sort the scripts to ensure consistent execution order
fig_scripts.sort()

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
