# Figure Index (SimSpace reproducibility)

This index maps each manuscript figure to the exact script and data used to reproduce it.

**Conventions**
- Figure scripts are under `main_figures/` and `supp_figures/`.
- Panel-specific inputs/caches are under `Panel_*_data/`.
- Panel-specific external code (separate env) is under `Panel_*_src/`.
- Example outputs are written under `example_output/` (as printed by the scripts).

**Recommended entry points**
- Main figures: `python plot_main.py`
- Supplementary figures: `python plot_supp.py`

---

## Main Figures

### Figure 1
- **Path:** `main_figures/Fig1/`
- **Entry script:** `main_figures/Fig1/Fig1.py`
- **Figure README:** `main_figures/Fig1/Fig1_README.md`
- **Environment:** default (`environment.yml`)
- **Run:**
  ```bash
  python main_figures/Fig1/Fig1.py
  ```

### Figure 2
- **Path:** `main_figures/Fig2/`
- **Entry script:** `main_figures/Fig2/Fig2.py`
- **Figure README:** `main_figures/Fig2/Fig2_README.md`
- **Panels / resources:**
  - **Panel D**
    - Data: `main_figures/Fig2/Panel_D_data/`
    - External code (if any): `main_figures/Fig2/Panel_D_src/`
    - Environment: default (`environment.yml`) and sccube (`main_figures/Fig2/sccube_environment.yaml`)
- **Run:**
  ```bash
  python main_figures/Fig2/Fig2.py
  ```

### Figure 3
- **Path:** `main_figures/Fig3/`
- **Entry script:** `main_figures/Fig3/Fig3.py`
- **Figure README:** `main_figures/Fig3/Fig3_README.md`
- **Panels / resources:**
  - Panel-specific data under `main_figures/Fig3/Panel_*_data/` (or shared data folders)
  - External code under `main_figures/Fig3/Panel_*_src/` (separate env where noted)
- **Run:**
  ```bash
  python main_figures/Fig3/Fig3.py
  ```

### Figure 4
- **Path:** `main_figures/Fig4/`
- **Entry script:** `main_figures/Fig4/Fig4.py`
- **Figure README:** `main_figures/Fig4/Fig4_README.md`
- **Panels / resources:**
  - Panel-specific data under `main_figures/Fig4/Panel_*_data/` (or shared data folders)
  - External code under `main_figures/Fig4/Panel_*_src/` (separate env where noted)
- **Run:**
  ```bash
  python main_figures/Fig4/Fig4.py
  ```

### Figure 5
- **Path:** `main_figures/Fig5/`
- **Entry script:** `main_figures/Fig5/Fig5.py`
- **Figure README:** `main_figures/Fig5/Fig5_README.md`
- **Panels / resources:**
  - Panel-specific data under `main_figures/Fig5/Panel_*_data/` (or shared data folders)
  - External code under `main_figures/Fig5/Panel_*_src/` (separate env where noted)
- **Run:**
  ```bash
  python main_figures/Fig5/Fig5.py
  ```

### Figure 6
- **Path:** `main_figures/Fig6/`
- **Entry script:** `main_figures/Fig6/Fig6.py`
- **Figure README:** `main_figures/Fig6/Fig6_README.md`
- **Panels / resources:**
  - Panel-specific data under `main_figures/Fig6/Panel_*_data/` (or shared data folders)
  - External code under `main_figures/Fig6/Panel_*_src/` (separate env where noted)
- **Run:**
  ```bash
  python main_figures/Fig6/Fig6.py
  ```

### The new Figure 2 (provisional)
- **Manuscript placement:** after current Figure 6 pending final figure renumbering; existing figure numbers are unchanged
- **Path:** `main_figures/Fig2_new/`
- **Entry script:** `main_figures/Fig2_new/Fig2_new.py`
- **Figure README:** `main_figures/Fig2_new/Fig2_new_README.md`
- **Outputs:** `Panel_A.png`, `Panel_B.png`, and `Panel_C.png`; final assembly is maintained in the manuscript figure directory
- **Analysis:** four-niche cortex-inspired anchor, 5 x 5 thickness-by-layer-specificity design, and reference-free marker expression
- **Launcher integration:** included automatically by `plot_main.py` through its `Fig*` discovery; existing numeric figure labels remain unchanged
- **Run:**
  ```bash
  python main_figures/Fig2_new/Fig2_new.py
  ```

### The new Figure 5 (provisional)
- **Manuscript placement:** intended to replace current Figure 5; final numbering is pending
- **Path:** `main_figures/Fig5_new/`
- **Entry script:** `main_figures/Fig5_new/Fig5_new.py`
- **Figure README:** `main_figures/Fig5_new/Fig5_new_README.md`
- **Outputs:** complete `Fig5_new.png/.svg` plus standalone `Panel_A`, `Panel_B`, and `Panel_C` PNG/SVG files
- **Analysis:** two-pattern, three-signal-level, three-seed spatial-domain benchmark of six named methods plus an expression-only negative control
- **Panel data:** compact map cache under `Panel_A_B_data/` and frozen aggregate metrics under `Panel_C_data/`
- **External environments/source:** `Panel_A_B_C_src/`; full method reruns use two dedicated Python environments and the host R library
- **Launcher integration:** included automatically by `plot_main.py` through its `Fig*` discovery
- **Run:**
  ```bash
  python main_figures/Fig5_new/Fig5_new.py
  ```

---

## Supplementary Figures

### Supplementary Figure 1
- **Path:** `supp_figures/SFig1/`
- **Entry script:** `supp_figures/SFig1/SFig1.py`
- **Panels:** B1/B2, 220-gene runtime/peak RAM; C1/C2, 2,200-gene runtime/peak RAM
- **Example output:** `example_output/SFig1/`
- **Environment:** default unless otherwise noted
- **Run:**
  ```bash
  python supp_figures/SFig1/SFig1.py
  ```

### Supplementary Figure 2
- **Path:** `supp_figures/SFig2/`
- **Entry script:** `supp_figures/SFig2/SFig2.py`
- **Panels / resources:**
  - Panel-specific data under `supp_figures/SFig2/Panel_*_data/` (or shared data folders)
  - External code under `supp_figures/SFig2/Panel_*_src/` (separate env where noted)
- **Run:**
  ```bash
  python supp_figures/SFig2/SFig2.py
  ```

### Supplementary Figure 3–9
- **Path:** `supp_figures/SFig3/` … `supp_figures/SFig9/`
- **Entry scripts:** `SFig3.py` … `SFig9.py` (as implemented)
- **Environment:** default unless otherwise noted

### Supplementary Figure 8
- **Path:** `supp_figures/SFig8/`
- **Entry script:** `supp_figures/SFig8/SFig8.py`
- **Figure README:** `supp_figures/SFig8/SFig8_README.md`
- **Inputs:** frozen Figure 3 Xenium tile, fitted SimSpace parameters,
  seed-0 molecular realization, ten-seed scDesign3 molecular summaries,
  and BANKSY outputs
- **Panels:** A–G, local spatial and niche validation; H–I, seed-0
  gene-mean and variance fidelity; J–K, separate ten-seed PCC and RMSE
  molecular-robustness summaries
- **Run:**
  ```bash
  python supp_figures/SFig8/SFig8.py
  ```

### Supplementary Figure 9
- **Path:** `supp_figures/SFig9/`
- **Entry script:** `supp_figures/SFig9/SFig9.py`
- **Figure README:** `supp_figures/SFig9/SFig9_README.md`
- **Analysis:** phenotype Gibbs-sweep sensitivity for R1-8
- **Run:**
  ```bash
  python supp_figures/SFig9/SFig9.py
  ```

### Supplementary Figure 10
- **Path:** `supp_figures/SFig10/`
- **Entry script:** `supp_figures/SFig10/SFig10.py`
- **Full data-generation script:**
  `supp_figures/SFig10/SFig10_generate_data.py`
- **Figure README:** `supp_figures/SFig10/SFig10_README.md`
- **Analysis:** quantitative reference-free 3D spatial structure and
  truth-based niche recovery for R1-5
- **Example output:** `example_output/SFig10/SFig10.png`
- **Run:**
  ```bash
  # Fast: render from archived ten-seed outputs
  python supp_figures/SFig10/SFig10.py

  # Full: regenerate simulations, metrics, exports, and figure
  python supp_figures/SFig10/SFig10_generate_data.py
  ```

### Supplementary Figure 11
- **Path:** `supp_figures/SFig11/`
- **Entry script:** `supp_figures/SFig11/SFig11.py`
- **Download helper:** `supp_figures/SFig11/download_starmap.py`
- **Figure README:** `supp_figures/SFig11/SFig11_README.md`
- **Analysis:** post-generation comparison of the fixed reference-free
  Figure 2 cortex simulation with a compact mouse visual-cortex STARmap field
  for R2-M4
- **Panel data:** archived source, mapped cells, profiles, metrics, complete
  configuration, checksums, and software provenance under `Panel_A_E_data/`
- **Example output:** `example_output/SFig11/SFig11.png`
- **Launcher integration:** discovered automatically by `plot_supp.py`
- **Run:**
  ```bash
  python supp_figures/SFig11/SFig11.py
  ```

### Supplementary Figure 12
- **Path:** `supp_figures/SFig12/`
- **Entry script:** `supp_figures/SFig12/SFig12.py`
- **Figure README:** `supp_figures/SFig12/SFig12_README.md`
- **Analysis:** twenty-seed, 50 × 50-lattice and 500-gene validation of four
  optional spatial-effect bases, capture thinning, ambient background, and
  excess dropout for R2-M1/R2-m5; conditional gene truth is reusable for R2-M3
- **Outputs:** PNG, raw seed-level metrics, complete gene truth,
  backward-compatibility checks, and software/configuration provenance
- **Launcher integration:** discovered automatically by `plot_supp.py`
- **Run:**
  ```bash
  # Full: regenerate simulations, metrics, checks, and figure
  python supp_figures/SFig12/SFig12.py

  # Fast: render from archived outputs
  python supp_figures/SFig12/SFig12.py --render-only
  ```

---

## Entry-point scripts

### `plot_main.py`
Runs all main-figure scripts in order.
```bash
python plot_main.py
```

### `plot_supp.py`
Runs all supplementary-figure scripts in order.
```bash
python plot_supp.py
```
