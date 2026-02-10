# SimSpace Reproducibility Repository

This repository contains scripts and data to reproduce **all main and supplementary figures** for the SimSpace manuscript.

- Main figures: `main_figures/`
- Supplementary figures: `supp_figures/`
- Example outputs: `example_output/`
- Figure-to-script mapping: **`FIGURE_INDEX.md`**
- Default environment: `environment.yml`

---

## Quick start

### 1) Create the default environment
Using mamba/conda:
```bash
mamba env create -f environment.yml
mamba activate simspace-repro
```

### 2) Reproduce all main figures
```bash
python plot_main.py
```

### 3) Reproduce all supplementary figures
```bash
python plot_supp.py
```

---

## Reproducing a single figure

Each figure has its own folder and a single entry script.

Example (Figure 2):
```bash
python main_figures/Fig2/Fig2.py
```

See `FIGURE_INDEX.md` for the complete mapping.

---

## Repository layout

```
simspace-reproducibility/
  example_output/
    Fig1/
      Fig1_panel_A.png
      Fig1_panel_B1.png
      ...
    Fig2/
    ...
    SFig1/
      SFig1_panel_B1.png
      SFig1_panel_B2.png
    SFig2/
    ...
  main_figures/
    Fig1/
    Fig2/
      Fig2.py
      Panel_D_data/
      Panel_D_src/
      Panel_D.md
    Fig3/
      Fig3.py
      Panel_*_data/
      Panel_*_src/
      README.md
    ...
  supp_figures/
    SFig1/
    SFig2/
    ...
  plot_main.py
  plot_supp.py
  FIGURE_INDEX.md
  environment.yml
  README.md
  LICENSE
```

---

## Panel data vs. external code (`Panel_*_data` and `Panel_*_src`)

Some figures are organized by panel to keep dependencies and data clear.

- `Panel_*_data/`  
  Inputs and/or cached intermediates required to reproduce that panel.  
  These are intended to work with the **default** environment unless explicitly stated otherwise.

- `Panel_*_src/`  
  External code that is difficult to integrate into the default environment (e.g., torch/CUDA, specialized R packages).  
  Each `Panel_*_src/` directory should be self-contained and include:
  - a short `Fig*_README.md` with exact setup, run commands, and conda env/R session info
  - a single entrypoint script for each method/package where possible

Panels that require a separate environment are explicitly labeled in `FIGURE_INDEX.md` and in the figure-level `Fig*_README.md`.

---

## Outputs

By convention, scripts write outputs under the its own directory, and print the output paths upon completion. 

We also provided the expected results, storing in `example_output/`. One can compare the outputs with the example_output to ensure the reproducibility.

---

## Notes on determinism / reproducibility

We aim for reproducible figure regeneration under the provided environments.  
Some third-party tools may still exhibit minor run-to-run variation even with fixed seeds (e.g., parallelism, non-deterministic low-level libraries, or upstream package behavior).

To minimize variability:
- use the provided environment specs (`environment.yml` or per-panel env files)
- for methods/packages beyond SimSpace, follow the README in each figure folder
- keep seeds fixed where applicable
- prefer CPU execution if GPU nondeterminism is observed

Where relevant, we provide cached intermediate results / caches to avoid refitting and stabilize outputs.

---

## Troubleshooting

### Common issues
- **Package/version mismatch:** recreate the environment from the relevant env file.
- **Missing optional dependencies:** check the error message and install into the correct env.
- **Large files:** some intermediates may be excluded from git and generated/downloaded on first run.

### Getting help
If you encounter an issue reproducing a figure, please open a GitHub issue and include:
- figure name (e.g., Fig3, SFig2)
- the command you ran
- full traceback/error message
- OS + Python version
