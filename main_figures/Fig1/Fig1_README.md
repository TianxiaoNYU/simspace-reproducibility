# Figure 1

This folder contains everything needed to reproduce the code-generated panels in **Figure 1**

---

## Contents

- `Fig1.py`  
  Entry script to reproduce Figure 2 (generates the figure and/or panel outputs).

- `fitted_params.json`  
  SimSpace parameter file. Used for generating example plots in the workflow

---

## Environment

### Default (SimSpace)
```bash
# Run the line below if not creating the repro env before
# conda env create -f ../../environment.yml
conda activate simspace-repro
```

---

## Run Figure 1

From the repository root:
```bash
python main_figures/Fig1/Fig1.py
```

Or from this directory:
```bash
python Fig1.py
```

---
