## SimSpace conda env
python bench.py \
  --methods simspace \
  --grid-sizes 10 20 30 50 75 100 150 200 \
  --repeats 5 \
  --threads 1 \
  --n-genes 2200 \
  --outfile results_simspace_2200.csv

## sccube conda env
python generate_sccube_runtime_inputs.py
# For an existing generated_data/ directory:
python generate_sccube_runtime_inputs.py --verify-existing

python bench.py \
  --methods sccube \
  --grid-sizes 10 20 30 50 75 100 150 200 \
  --repeats 5 \
  --threads 1 \
  --n-genes 2200 \
  --sccube-data-dir generated_data \
  --outfile results_sccube_2200.csv

python bench.py \
  --methods scmultisim \
  --grid-sizes 10 12 15 \
  --repeats 5 \
  --threads 1 \
  --n-genes 2200 \
  --outfile results_scmultisim_2200.csv

The planned 18 x 18 scMultiSim benchmark was stopped manually because the
15 x 15 benchmark already required approximately 500 seconds per run. No
18 x 18 result is included in the analysis.

## Make plots
python SFig1.py \
  --csv-220 results.csv \
  --csv-2200 results_2200_gene.csv \
  --outdir ./

## Large generated scCube inputs

`data/` contains the tracked 220-gene source matrices.
`generate_sccube_runtime_inputs.py` deterministically creates the
2,200-feature runtime-only matrices under the Git-ignored
`generated_data/` directory and writes `generated_data_manifest.tsv`.
The expanded profiles are synthetic computational stress-test inputs,
not a biologically realistic expression panel.
