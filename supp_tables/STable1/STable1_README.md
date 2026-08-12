# Supplementary Table 1 capability audit

This directory contains the evidence audit supporting the finalized Supplementary Table 1 in the SimSpace revision manuscript.

## Files

- `STable1.tsv`: detailed evidence-audit source table. Each non-label cell ends with one or more evidence identifiers.
- `manuscript_matrix.tsv`: exact eight-column check/cross matrix used in Supplementary Table 1. Each cell maps to the evidence identifiers supporting the binary classification.
- `methods.tsv`: audited method/version metadata.
- `evidence.tsv`: traceability table linking evidence identifiers to versioned primary sources or inspected source code.
- `build_table.mjs`: reproducible workbook builder using the bundled `@oai/artifact-tool` runtime.
- `build_publication_pdf.py`: retired compatibility stub; the publication-layout table is maintained only in the manuscript TeX.

The detailed audit workbook is exported to `outputs/019f9047-2c3a-7852-b086-30ff9814c7d9/Supplementary_Table_1.xlsx`; it is raw supporting material rather than the publication-layout table.

## Interpretation rules

- `YES`: direct, documented native support in the audited version.
- `PARTIAL`: limited, indirect, optional-backend, or proxy support.
- `NO`: the audited source makes the capability incompatible with, or explicitly outside, the method's operating model.
- `NOT DESCRIBED`: the capability was not identified in the audited sources. This is evidence-bounded and is not a claim that implementation is impossible.
- `UNKNOWN`: sources were incomplete or conflicting.
- `N/A`: the capability does not conceptually apply to the method's documented simulation layer.
- `NOT RUNNABLE`: a capability could not be audited because the software could not be executed. This value is reserved but not used in the current table.

## Manuscript matrix derivation

The manuscript uses a deliberately compact binary matrix. `CHECK` requires direct support for the exact column definition in the cited article or audited software version. `CROSS` includes `PARTIAL`, `NO`, `NOT DESCRIBED`, `UNKNOWN`, and `N/A`, as well as capabilities that are available only indirectly or outside the exact column definition. This conservative collapse prevents a partial or adjacent feature from being displayed as equivalent direct support.

The exact column definitions are:

- **Reference-free mode:** the method can generate a spatial omics dataset without fitting or conditioning on an empirical reference dataset.
- **Reference-based mode:** empirical data, including observed coordinates or spatial covariates when required by the method, are used to fit or condition at least one spatial or molecular simulation component. This broad capability does not by itself require the method to generate new coordinates.
- **Molecular data simulation:** molecular profiles or attributes are generated within the framework or through an integrated simulation engine.
- **3-D spatial generation:** the documented simulation workflow natively generates three-coordinate spatial layouts rather than merely accepting or plotting existing 3-D coordinates.
- **Pairwise cell-type control:** users can directly specify or modify pair-specific spatial association, transition, affinity, or co-localization behavior between cell types.
- **Niche-to-cell model:** an explicit generated tissue-domain process is followed by a cell-type process conditioned on that domain layer.
- **Fitted-model perturbation:** fitted niche or niche-conditioned cell-type parameters can be modified while other components of that fitted model are held fixed; permutation of observed data is insufficient.
- **Validated reference-based layouts:** a reference-calibrated two-stage niche-to-cell-type model generates new coordinates, and the resulting layouts are evaluated quantitatively for reference spatial similarity and biologically meaningful domain and cell-type organization. A fitted tissue outline, transferred expression ranks, cell-type point processes without a generated niche layer, or a global transition matrix without a fitted niche layer is insufficient.

The exact displayed column order is: reference-free mode, reference-based mode, molecular data simulation, 3-D spatial generation, pairwise cell-type control,
niche-to-cell model, fitted-model perturbation, and validated reference-based layouts.

The matrix is a capability audit, not a global performance ranking. Quantitative performance questions are addressed separately under their matched reviewer-comment analyses.

## Scope

This is a capability and interface audit, not an overall performance ranking. The compact headers in the manuscript are shorthand for the exact rules above. A cross therefore means that the audited method does not meet that full column definition; it does not imply that the method lacks every related or adjacent function.

SimSpace is audited at version 0.4.0. This release adds the documented headless command-line interface and stable CLI outputs while leaving the production simulation functions used by the manuscript workflows unchanged.
