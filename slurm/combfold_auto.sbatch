#!/bin/bash
#SBATCH --job-name=combfold
#SBATCH --account=rp24
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.log
# CPU-only job. Usage:
#   sbatch slurm/combfold_auto.sbatch <out_dir> <structure1.cif> <structure2.cif> [...]
# Extra combfold_auto.py flags can be passed via the COMBFOLD_ARGS env var, e.g.
#   COMBFOLD_ARGS="--names names.fasta --top-n 10" sbatch slurm/combfold_auto.sbatch out/ jobs/*/*_model.cif
set -euo pipefail

out_dir=$1
shift

eval "$(micromamba shell hook --shell bash)"
micromamba activate combfold

repo_dir=${COMBFOLD_DIR:-$(dirname "$(dirname "$(readlink -f "$0")")")}
python "$repo_dir/scripts/combfold_auto.py" "$@" -o "$out_dir" --threads "$SLURM_CPUS_PER_TASK" ${COMBFOLD_ARGS:-}
