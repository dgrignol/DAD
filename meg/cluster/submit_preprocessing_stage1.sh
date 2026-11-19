#!/bin/bash
#
# submit_preprocessing_stage1.sh
#
# Purpose
# Run the first cohort review gate: source preflight, session-note and FIF
# header inventory, STI101 trigger extraction/alignment, and the full deep
# trigger report for every participant. Visual timing, including P24 timing,
# is deferred to the cohort-wide visual-timing gate. The job stops before
# Maxwell filtering, artifact decisions, ICA, epoching, or visual-time
# correction application.
#
# Usage from the deployed DAD root on the DIPSEN login node
#   mkdir -p cluster_logs
#   sbatch --array=4-39 meg/cluster/submit_preprocessing_stage1.sh
# Optional P24-only rerun (not needed solely to remove the old timing proposal):
#   sbatch --array=24 meg/cluster/submit_preprocessing_stage1.sh
# The default SLURM array is 4-39 with at most four tasks running at once.
# Inputs: read-only Attention/data files and repository configurations.
# Outputs: Attention/derivatives/reports and preprocessing/events only.
# Data flow: SLURM task ID -> subject config -> preflight/Stage 02 -> deep report.
# Rerun behavior: existing selected-subject Stage-1 outputs are rewritten by
# the runner, including inventory/provenance and deep plots. --resume does not
# protect trigger-stage outputs. Deploying this scope change requires no rerun.

#SBATCH --job-name=DAD_meg_stage1
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --export=ALL
#SBATCH --array=4-39

set -euo pipefail

# Resolve project paths outside the read-only raw-data tree.
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?Submit this script as an SLURM array job.}"
if (( TASK_ID < 4 || TASK_ID > 39 )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID must be 4-39; got $TASK_ID" >&2
    exit 2
fi
SUBJECT_ID="$(printf '%02d' "$TASK_ID")"
if [[ ! -x "$MAMBA_BIN" ]]; then
    echo "ERROR: micromamba not found: $MAMBA_BIN" >&2
    exit 2
fi

# Activate the cluster environment with isolated per-task caches.
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "$MAMBA_ENV"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
JOB_TMP="${TMPDIR:-/tmp}/dad_stage1_${SLURM_JOB_ID:-manual}_${TASK_ID}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib" XDG_CACHE_HOME="$JOB_TMP/cache"

# Stage 02 writes canonical alignment plus the full deep report, anomaly
# plots, and transition plots for independent review before the next stage.
cd "$PROJECT_ROOT"
mkdir -p cluster_logs
CONFIG="$PROJECT_ROOT/meg/params/pipeline_default.yaml"
SITE="$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
SUBJECT="$PROJECT_ROOT/meg/params/sub-${SUBJECT_ID}.yaml"
for required in "$CONFIG" "$SITE" "$SUBJECT"; do
    [[ -f "$required" ]] || { echo "ERROR: missing $required" >&2; exit 2; }
done
python "$PROJECT_ROOT/meg/code/preprocessing_basic/run_pipeline.py" \
    --config "$CONFIG" --site-config "$SITE" --subject-config "$SUBJECT" \
    --stages trigger_extraction --write-inventory

