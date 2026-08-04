#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "$(dirname "$0")/.."
if ! command -v micromamba >/dev/null 2>&1; then
  export PATH="$PWD/.tools/bin:$PATH"
fi

RUNNER=(micromamba run -p .mamba_vina)

"${RUNNER[@]}" python scripts/01_check_inputs.py
"${RUNNER[@]}" python scripts/02_fetch_receptors.py
"${RUNNER[@]}" python scripts/03_prepare_ligands.py --limit 30000
"${RUNNER[@]}" python scripts/04_run_vina.py
"${RUNNER[@]}" python scripts/05_analyze_enopt.py

echo "Full β2-AR screening completed. Check results/tables/enopt_style_top_hits.csv"

