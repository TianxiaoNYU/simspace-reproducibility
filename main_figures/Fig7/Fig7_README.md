# Figure 7

This directory reproduces the code-generated panels for Figure 7. The final
multi-panel manuscript figure is assembled manually from these outputs.

## Environment

Figure 7 is reproduced with SimSpace 0.4.2, which provides the
continuous-intensity support used by the scDesign3 adapter. Create or activate
the repository environment and verify the installed distribution:

```bash
conda env create -f environment.yml
conda activate simspace-repro
python -c "from importlib.metadata import version; print(version('simspace'))"
```

The figure's version guard uses the installed distribution metadata. The local
R environment must contain the scDesign3 dependencies described in the main
reproducibility README.

## Run

From the repository root:

```bash
python main_figures/Fig7/Fig7.py
```

The script uses spatial seed 1 and molecular seed 0. It explicitly selects the
Gaussian pathway for the continuous CODEX marker intensities; the execution
log should report `family 'gaussian'`, assay `logcounts`, and `DT=FALSE`.

## Inputs

- `Panel_A_B_C_data/simspace_fitted.json` contains the fitted spatial
  parameters. The figure script generates a new layout on a 50-by-50 lattice
  using four cell-type sweeps and six niche sweeps.
- `Panel_A_B_C_data/CODEX_expr_mat.csv` contains 31 nonnegative continuous
  protein-marker intensities for 206 cells.
- `Panel_A_B_C_data/CODEX_meta.csv` contains the coordinates and nine phenotype
  annotations used consistently for both model fitting and reference maps.
- `Panel_A_B_C_data/CODEX_ref.csv` is retained as a legacy visualization table
  but is not used by the current script because its labels differ from the
  fitting annotations for 41 cells.
- `Panel_A_B_C_data/perturb/ss_perturb_1.json` through
  `ss_perturb_7.json` vary only the fitted cytotoxic-T-cell--helper-T-cell
  affinity over 0.4, 0.7, 0.8, 0.9, 1.0, 1.2, and 1.4.

The processed CODEX/PhenoCycler-Fusion data and related project scripts are
available from <https://github.com/FenyoLab/EC_codeximaging>.

## Outputs

The script writes the following files into this directory:

- `Fig7_Panel_A1.png`: reference phenotype map.
- `Fig7_Panel_A2.png`: reference PD1-intensity map.
- `Fig7_Panel_A3.png`: simulated PD1-intensity map.
- `Fig7_Panel_A4.png`: simulated phenotype map.
- `Fig7_Panel_B1.png`: fitted affinity heatmap.
- `Fig7_Panel_B2.png` and `Fig7_Panel_B3.png`: simulated and reference T-cell
  maps.
- `Fig7_Panel_B4.png` and `Fig7_Panel_B5.png`: simulated and reference
  endothelial/stromal maps.
- `Fig7_Panel_C.png`: seven affinity-perturbation layouts.

The script verifies that simulated marker intensities are finite, nonnegative,
and continuous. Curated copies of the panel outputs used for assembly are
stored under `example_output/Fig7/`.
