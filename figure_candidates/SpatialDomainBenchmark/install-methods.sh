#!/usr/bin/env bash
set -euo pipefail

expected_env="spatial-domain-benchmark"
if [[ "${CONDA_DEFAULT_ENV:-}" != "$expected_env" ]]; then
  echo "Activate $expected_env before running this script." >&2
  exit 2
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Host R is required on PATH; do not install r-base in this Conda environment." >&2
  exit 2
fi

# mclust belongs to the user's normal R library, not to Conda.
if ! Rscript -e 'quit(status = if (requireNamespace("mclust", quietly = TRUE)) 0 else 1)'; then
  Rscript -e 'install.packages("mclust", repos = "https://cloud.r-project.org")'
fi

# ABI mode avoids compiling/linking a private R runtime and binds rpy2 to
# the host R discovered on PATH. It also works on Apple Silicon.
# Use the platform's generic compiler names for rpy2's build-time probe;
# this also avoids stale user-level CC/CXX paths leaking into the build.
env RPY2_CFFI_MODE=ABI CC=cc CXX=c++ python -m pip install \
  'rpy2==3.5.17' \
  'scanyuan==0.0.8' \
  'mygene==3.2.2'

source_root="${CONDA_PREFIX}/src/spatial-domain-methods"
mkdir -p "$source_root"

clone_at_commit() {
  local name="$1"
  local url="$2"
  local commit="$3"
  local destination="${source_root}/${name}"

  if [[ ! -d "${destination}/.git" ]]; then
    git clone "$url" "$destination"
  fi
  git -C "$destination" fetch --depth 1 origin "$commit"
  git -C "$destination" checkout --detach "$commit"
}

clone_at_commit spCLUE \
  https://github.com/EnchantedJoy/spCLUE.git \
  bbd2c342e7a67c1617275f721cec2e3f4c23a799

# Install directly from immutable commits without allowing upstream's broad
# dependency declarations to replace the tested Conda pins. Pip uses temporary
# clones, so repeated runs cannot dirty the persistent source checkouts.
python -m pip install --no-deps --force-reinstall \
  'GraphST @ git+https://github.com/JinmiaoChenLab/GraphST.git@d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74' \
  'STAGATE_pyG @ git+https://github.com/QIFEIDKN/STAGATE_pyG.git@ae1158ca8cf1eb6bb8ee198298552d44c9ac21db' \
  'SpaGCN @ git+https://github.com/jianhuupenn/SpaGCN.git@dc7a1c26ea0fdf4dfe7064adc7699be141b4871f#subdirectory=SpaGCN_package'

# spCLUE is source-only upstream. Expose its package directory without
# altering upstream source or inventing an unpublished wheel.
site_packages="$(python -c 'import site; print(site.getsitepackages()[0])')"
if [[ -e "${site_packages}/spCLUE" && ! -L "${site_packages}/spCLUE" ]]; then
  echo "Refusing to overwrite existing ${site_packages}/spCLUE" >&2
  exit 2
fi
ln -sfn "${source_root}/spCLUE/spCLUE" "${site_packages}/spCLUE"

# `pip check` is not authoritative for a mixed Conda/PyPI environment:
# Conda's python-igraph distribution intentionally imports as `igraph`, and
# several Conda builds carry placeholder wheel metadata. The direct tests
# below exercise the imports and compiled extensions that the methods use.
python "$(dirname "$0")/smoke-test.py"
