import numpy as np
import pandas as pd
import typing
import subprocess
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

def run_CARD(
    input_meta_file: str,
    input_omics_file: str,
    ref_meta_file: str = None,
    ref_omics_file: str = None,
    # threads: int = 6,
    output_dir: str = f"{script_dir}",
    ):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    result = subprocess.run(
        [
            "Rscript", 
            f'{script_dir}/CARD.R', 
            input_meta_file,
            input_omics_file,
            ref_meta_file,
            ref_omics_file,
            output_dir
        ],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None
    else:
        print("CARD deconvolution complete.")
        return result.stdout
    
def run_RCTD(
    input_meta_file: str,
    input_omics_file: str,
    ref_meta_file: str = None,
    ref_omics_file: str = None,
    threads: int = 6,
    output_dir: str = f"{script_dir}",
    ):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    result = subprocess.run(
        [
            "Rscript", 
            f'{script_dir}/RCTD.R', 
            input_meta_file,
            input_omics_file,
            ref_meta_file,
            ref_omics_file,
            str(threads),
            output_dir
        ],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None
    else:
        print("RCTD deconvolution complete.")
        return result.stdout
    
def run_Seurat(
    input_meta_file: str,
    input_omics_file: str,
    ref_meta_file: str = None,
    ref_omics_file: str = None,
    output_dir: str = f"{script_dir}",
    ):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    result = subprocess.run(
        [
            "Rscript", 
            f'{script_dir}/Seurat.R', 
            input_meta_file,
            input_omics_file,
            ref_meta_file,
            ref_omics_file,
            output_dir
        ],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None
    else:
        print("Seurat deconvolution complete.")
        return result.stdout
    
def run_STdeconvolve(
    input_meta_file: str,
    input_omics_file: str,
    output_dir: str = f"{script_dir}",
    ):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    result = subprocess.run(
        [
            "Rscript", 
            f'{script_dir}/STdeconvolve.R', 
            input_meta_file,
            input_omics_file,
            output_dir
        ],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None
    else:
        print("STdeconvolve deconvolution complete.")
        return result.stdout
    
def run_spatialDWLS(
    input_meta_file: str,
    input_omics_file: str,
    ref_meta_file: str = None,
    ref_omics_file: str = None,
    output_dir: str = f"{script_dir}",
    ):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    result = subprocess.run(
        [
            "Rscript", 
            f'{script_dir}/spatialDWLS.R', 
            input_meta_file,
            input_omics_file,
            ref_meta_file,
            ref_omics_file,
            output_dir
        ],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None
    else:
        print("spatialDWLS deconvolution complete.")
        return result.stdout

def run_cell2location(
    input_meta_file: str,
    input_omics_file: str,
    ref_meta_file: str = None,
    ref_omics_file: str = None,
    output_dir: str = f"{script_dir}",
    ):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    result = subprocess.run(
        [
            "python", 
            f'{script_dir}/cell2loc.py', 
            "--input_meta_file", input_meta_file,
            "--input_omics_file", input_omics_file,
            "--ref_meta_file", ref_meta_file,
            "--ref_omics_file", ref_omics_file,
            "--output_dir", output_dir
        ],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("Python script failed:")
        print(result.stderr)
        return None
    else:
        print("cell2location deconvolution complete.")
        return result.stdout


def test():
    pass

if __name__ == "__main__":
    test()