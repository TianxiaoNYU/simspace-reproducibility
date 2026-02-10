# Figure Index (SimSpace reproducibility)

This index maps each manuscript figure to the exact script and data used to reproduce it.

**Conventions**
- Figure scripts are under `main_figures/` and `supp_figures/`.
- Panel-specific inputs/caches are under `Panel_*_data/`.
- Panel-specific external code (separate env) is under `Panel_*_src/`.
- Outputs are written under `example_output/` (or as printed by the scripts).

**Recommended entry points**
- Main figures: `python plot_main.py`
- Supplementary figures: `python plot_supp.py`

---

## Main Figures

### Figure 1
- **Path:** `main_figures/Fig1/`
- **Entry script:** `main_figures/Fig1/Fig1.py`
- **Fitted parameters:** `main_figures/Fig1/fitted_params.json`
- **Figure README:** `main_figures/Fig1/README.md`
- **Environment:** default (`environment.yml`)
- **Run:**
  ```bash
  python main_figures/Fig1/Fig1.py
  ```

### Figure 2
- **Path:** `main_figures/Fig2/`
- **Entry script:** `main_figures/Fig2/Fig2.py`
- **Panels / resources:**
  - **Panel D**
    - Data: `main_figures/Fig2/Panel_D_data/`
    - External code (if any): `main_figures/Fig2/Panel_D_src/`
    - Panel runbook: `main_figures/Fig2/Panel_D.md`
    - Environment: default **unless** Panel_D_src specifies otherwise
- **Run:**
  ```bash
  python main_figures/Fig2/Fig2.py
  ```
- **Notes:**
  - If Panel D requires a separate environment, follow `Panel_D_src/README.md` first.

### Figure 3
- **Path:** `main_figures/Fig3/`
- **Entry script:** `main_figures/Fig3/Fig3.py`
- **Figure README:** `main_figures/Fig3/README.md`
- **Panels / resources:**
  - Panel-specific data under `main_figures/Fig3/Panel_*_data/` (or shared data folders)
  - External code under `main_figures/Fig3/Panel_*_src/` (separate env where noted)
- **Run:**
  ```bash
  python main_figures/Fig3/Fig3.py
  ```

### Figure 4
- **Path:** `main_figures/Fig4/`
- **Entry script:** `main_figures/Fig4/Fig4.py` (or as implemented)
- **Environment:** default unless otherwise noted
- **Run:**
  ```bash
  python main_figures/Fig4/Fig4.py
  ```

### Figure 5
- **Path:** `main_figures/Fig5/`
- **Entry script:** `main_figures/Fig5/Fig5.py` (or as implemented)
- **Environment:** default unless otherwise noted
- **Run:**
  ```bash
  python main_figures/Fig5/Fig5.py
  ```

### Figure 6
- **Path:** `main_figures/Fig6/`
- **Entry script:** `main_figures/Fig6/Fig6.py` (or as implemented)
- **Environment:** default unless otherwise noted
- **Run:**
  ```bash
  python main_figures/Fig6/Fig6.py
  ```

---

## Supplementary Figures

### Supplementary Figure 1
- **Path:** `supp_figures/SFig1/`
- **Entry script:** `supp_figures/SFig1/SFig1.py`
- **Environment:** default unless otherwise noted
- **Run:**
  ```bash
  python supp_figures/SFig1/SFig1.py
  ```

### Supplementary Figure 2
- **Path:** `supp_figures/SFig2/`
- **Entry script:** `supp_figures/SFig2/SFig2.py`
- **Inputs:**  
  - `supp_figures/SFig2/Panel_D_data/`
  - `supp_figures/SFig2/Panel_E_F_data/`
- **Environment:** default unless otherwise noted
- **Run:**
  ```bash
  python supp_figures/SFig2/SFig2.py
  ```

### Supplementary Figure 3–7
- **Path:** `supp_figures/SFig3/` … `supp_figures/SFig7/`
- **Entry scripts:** `SFig3.py` … `SFig7.py` (as implemented)
- **Environment:** default unless otherwise noted

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
