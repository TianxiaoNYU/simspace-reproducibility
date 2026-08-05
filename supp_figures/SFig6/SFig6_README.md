# Supplementary Figure 6 — reference-free mouse-cortex comparison

This analysis addresses reviewer comment R2-M4 by comparing the fixed
reference-free cortex simulation from Figure 3 with an observed mouse
visual-cortex STARmap field. The observed data are used only after simulation;
they are not used to fit, tune, or select the simulated anchor.

## Dataset

The reviewer suggested STARmap, and a compact processed visual-cortex sample
is available under CC BY 4.0 at Zenodo DOI
[`10.5281/zenodo.10698912`](https://doi.org/10.5281/zenodo.10698912).
The 471,002-byte archive contains a 1,207-cell by 1,020-gene raw-count
matrix, spatial coordinates, cortical-domain labels, and cell annotations.
The original STARmap study is Wang et al., *Science* 2018, DOI
[`10.1126/science.aat5691`](https://doi.org/10.1126/science.aat5691).

The downloaded archive is retained in `Panel_A_E_data/` so the normal figure
entry point runs without network access. `download_starmap.py` independently
reproduces the download and verifies its byte count, MD5, and SHA-256 digest.

## Design

The analysis retains cells with an explicit excitatory-layer annotation
(`eL2/3`, `eL4`, `eL5`, `eL6-1`, or `eL6-2`) and an explicit cortical-domain
label (`L2/3`, `L4`, `L5`, or `L6`), leaving 601 STARmap cells. Observed and
simulated labels are collapsed to four declared classes: superficial, middle,
deep, and inner.

Spatial structure is compared through row-normalized 4 × 4
cell-class-by-domain matrices. Expression is compared for all eight named
Figure 3 markers present in STARmap: *Cux2*, *Rorb*, *Fezf2*, *Deptor*,
*Tshz2*, *Foxp2*, *Sulf2*, and *Ctgf*. For each marker, domain-wise mean counts
are normalized within each dataset to sum to one; absolute assay count scales
are not compared.

Panels A and B compare the spatial cell-class organization. Panel C shows
spatial expression maps for *Cux2* and *Rorb* as representative, visually
legible examples. Panels D and E retain the complete prespecified eight-gene
comparison, including discordant profiles.

Cells and spatial domains within the single observed field are not treated as
independent biological replicates. The correlations summarize agreement of
the declared coarse profiles rather than sampling uncertainty or biological
replication.

## Run

From the reproducibility-repository root:

```bash
conda activate simspace-repro
python main_figures/Fig3/Fig3.py
python supp_figures/SFig6/SFig6.py
```

If the archived STARmap input needs to be reacquired:

```bash
python supp_figures/SFig6/download_starmap.py
```

The entry script writes:

- `supp_figures/SFig6/SFig6.png`;
- `example_output/SFig6/SFig6.png`;
- mapped observed and simulated cells under `Panel_A_E_data/`;
- cell-class/domain and marker-profile tables;
- aggregate and per-gene comparison metrics; and
- the complete analysis configuration, input manifest, checksums, and
  software versions.

## Result and interpretation limits

The cell-class-by-domain matrices agree strongly (Pearson `r = 0.967`, RMSE
`= 0.079`). Across all eight shared markers, the profile-matrix Pearson
correlation is `r = 0.524`, the median per-marker Pearson correlation is
`r = 0.721`, and observed and simulated peak domains match for five of eight
markers.

The profiles of *Cux2*, *Rorb*, *Fezf2*, and *Deptor* are highly concordant,
and *Ctgf* has the same peak domain with weaker profile agreement. *Tshz2*,
*Foxp2*, and *Sulf2* are discordant in this observed field. The figure therefore
supports strong coarse laminar agreement and heterogeneous molecular agreement;
it does not claim reconstruction of the particular STARmap specimen or broad
agreement for every cortical marker.
