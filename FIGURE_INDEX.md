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
- **Run:** `python main_figures/Fig1/Fig1.py`

### Figure 2
- **Path:** `main_figures/Fig2/`
- **Entry script:** `main_figures/Fig2/Fig2.py`
- **Figure README:** `main_figures/Fig2/Fig2_README.md`
- **Panel D resources:** `main_figures/Fig2/Panel_D_data/` and `main_figures/Fig2/Panel_D_src/`
- **Run:** `python main_figures/Fig2/Fig2.py`

### Figure 3
- **Path:** `main_figures/Fig3/`
- **Entry script:** `main_figures/Fig3/Fig3.py`
- **Figure README:** `main_figures/Fig3/Fig3_README.md`
- **Analysis:** four-niche cortex-inspired anchor, 5 × 5 thickness-by-layer-specificity design, and reference-free marker expression
- **Panel data:** `Panel_A_C_data/`, `Panel_B_data/`, and `Panel_C_data/`
- **Run:** `python main_figures/Fig3/Fig3.py`

### Figure 4
- **Path:** `main_figures/Fig4/`
- **Entry script:** `main_figures/Fig4/Fig4.py`
- **Figure README:** `main_figures/Fig4/Fig4_README.md`
- **Panels/resources:** `main_figures/Fig4/Panel_*_data/` and `main_figures/Fig4/Panel_*_src/`
- **Run:** `python main_figures/Fig4/Fig4.py`

### Figure 5
- **Path:** `main_figures/Fig5/`
- **Entry script:** `main_figures/Fig5/Fig5.py`
- **Figure README:** `main_figures/Fig5/Fig5_README.md`
- **Panels/resources:** `main_figures/Fig5/Panel_*_data/` and `main_figures/Fig5/Panel_*_src/`
- **Run:** `python main_figures/Fig5/Fig5.py`

### Figure 6
- **Path:** `main_figures/Fig6/`
- **Entry script:** `main_figures/Fig6/Fig6.py`
- **Figure README:** `main_figures/Fig6/Fig6_README.md`
- **Analysis:** two-pattern, three-signal-level, three-seed spatial-domain benchmark of six named methods plus an expression-only negative control
- **Panel data:** `Panel_A_B_data/` and `Panel_C_data/`
- **External environments/source:** `Panel_A_B_C_src/`
- **Run:** `python main_figures/Fig6/Fig6.py`

### Figure 7
- **Path:** `main_figures/Fig7/`
- **Entry script:** `main_figures/Fig7/Fig7.py`
- **Figure README:** `main_figures/Fig7/Fig7_README.md`
- **Panels/resources:** `main_figures/Fig7/Panel_*_data/` and `main_figures/Fig7/Panel_*_src/`
- **Run:** `python main_figures/Fig7/Fig7.py`

---

## Supplementary Figures

### Supplementary Figure 1
- **Path:** `supp_figures/SFig1/`
- **Entry script:** `supp_figures/SFig1/SFig1.py`
- **Figure README:** `supp_figures/SFig1/SFig1_README.md`
- **Panels:** B1/B2, 220-gene runtime/peak RAM; C1/C2, 2,200-gene runtime/peak RAM
- **Run:** `python supp_figures/SFig1/SFig1.py`

### Supplementary Figure 2
- **Path:** `supp_figures/SFig2/`
- **Entry script:** `supp_figures/SFig2/SFig2.py`
- **Figure README:** `supp_figures/SFig2/SFig2_README.md`
- **Run:** `python supp_figures/SFig2/SFig2.py`

### Supplementary Figure 3
- **Path:** `supp_figures/SFig3/`
- **Entry script:** `supp_figures/SFig3/SFig3.py`
- **Figure README:** `supp_figures/SFig3/SFig3_README.md`
- **Analysis:** validation of optional direct spatial-expression and observation layers across twenty seeds
- **Run:** `python supp_figures/SFig3/SFig3.py --render-only`

### Supplementary Figure 4
- **Path:** `supp_figures/SFig4/`
- **Entry script:** `supp_figures/SFig4/SFig4.py`
- **Figure README:** `supp_figures/SFig4/SFig4_README.md`
- **Analysis:** sensitivity to spatial smoothness, neighborhood size, and coordinate perturbation
- **Run:** `python supp_figures/SFig4/SFig4.py`

### Supplementary Figure 5
- **Path:** `supp_figures/SFig5/`
- **Entry script:** `supp_figures/SFig5/SFig5.py`
- **Figure README:** `supp_figures/SFig5/SFig5_README.md`
- **Analysis:** phenotype-level Gibbs-sweep sensitivity
- **Run:** `python supp_figures/SFig5/SFig5.py`

### Supplementary Figure 6
- **Path:** `supp_figures/SFig6/`
- **Entry script:** `supp_figures/SFig6/SFig6.py`
- **Download helper:** `supp_figures/SFig6/download_starmap.py`
- **Figure README:** `supp_figures/SFig6/SFig6_README.md`
- **Analysis:** post-generation comparison of the fixed Figure 3 cortex simulation with a mouse visual-cortex STARmap field
- **Run:** `python supp_figures/SFig6/SFig6.py`

### Supplementary Figure 7
- **Path:** `supp_figures/SFig7/`
- **Entry script:** `supp_figures/SFig7/SFig7.py`
- **Full data-generation script:** `supp_figures/SFig7/SFig7_generate_data.py`
- **Figure README:** `supp_figures/SFig7/SFig7_README.md`
- **Analysis:** quantitative reference-free 3-D spatial structure and truth-based niche recovery
- **Run:** `python supp_figures/SFig7/SFig7.py`

### Supplementary Figure 8
- **Path:** `supp_figures/SFig8/`
- **Entry script:** `supp_figures/SFig8/SFig8.py`
- **Figure README:** `supp_figures/SFig8/SFig8_README.md`
- **Analysis:** additional evaluation of reference-based SimSpace simulations
- **Run:** `python supp_figures/SFig8/SFig8.py`

### Supplementary Figure 9
- **Path:** `supp_figures/SFig9/`
- **Entry script:** `supp_figures/SFig9/SFig9.py`
- **Figure README:** `supp_figures/SFig9/SFig9_README.md`
- **Inputs:** Figure 4 Xenium tile; ten independent spatial fits (optimizer seeds 0–9); molecular realizations; and BANKSY outputs
- **Analysis:** local independent-calibration spatial robustness and molecular, co-localization, and niche fidelity
- **Run:** follow the staged Python–R–Python commands in `supp_figures/SFig9/SFig9_README.md`

### Supplementary Figure 10
- **Path:** `supp_figures/SFig10/`
- **Entry script:** `supp_figures/SFig10/SFig10.py`
- **Figure README:** `supp_figures/SFig10/SFig10_README.md`
- **Analysis:** example simulated datasets and their convolved spot-level representations
- **Run:** `python supp_figures/SFig10/SFig10.py`

### Supplementary Figure 11
- **Path:** `supp_figures/SFig11/`
- **Entry script:** `supp_figures/SFig11/SFig11.py`
- **Figure README:** `supp_figures/SFig11/SFig11_README.md`
- **Analysis:** cell-type deconvolution performance across kernel sizes
- **Run:** `python supp_figures/SFig11/SFig11.py`

### Supplementary Figure 12
- **Path:** `supp_figures/SFig12/`
- **Entry script:** `supp_figures/SFig12/SFig12.py`
- **Figure README:** `supp_figures/SFig12/SFig12_README.md`
- **Analysis:** SVG-detection benchmarking under reference-based and reference-free simulations
- **Run:** `python supp_figures/SFig12/SFig12.py`

### Supplementary Figure 13
- **Path:** `supp_figures/SFig13/`
- **Entry script:** `supp_figures/SFig13/SFig13.py`
- **Figure README:** `supp_figures/SFig13/SFig13_README.md`
- **Analysis:** SVG-rank concordance between SimSpace simulations and the Xenium reference
- **Run:** `python supp_figures/SFig13/SFig13.py`

---

## Entry-point scripts

### `plot_main.py`
Runs all main-figure scripts in numeric order.

```bash
python plot_main.py
```

### `plot_supp.py`
Runs all supplementary-figure scripts in numeric order.

```bash
python plot_supp.py
```
