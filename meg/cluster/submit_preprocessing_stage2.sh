#!/bin/bash
#
# Purpose: Stage 2 trial accounting for P04-P39, before epoch creation.
# Usage from the deployed DAD root:
#   mkdir -p cluster_logs
#   stage2_job=$(sbatch --parsable --array=4-39 meg/cluster/submit_preprocessing_stage2.sh)
#   sbatch --dependency="afterok:${stage2_job}" meg/cluster/submit_preprocessing_stage2_summary.sh
# Inputs: accepted Stage 1 derivatives, read-only raw FIF headers, configuration,
# and PREPROCESSING_DECISIONS.md. Outputs: dedicated trial_preparation bundles.
# Data flow: array task -> participant selector -> trial_preparation.py ledger/QC.
# No runner/preflight regeneration, timing correction, signal cleaning, or epochs.
# Existing outputs cause failure. DAD_STAGE2_OVERWRITE=1 explicitly replaces only
# named Stage 2 outputs for the selected tasks; default submissions never do so.
# DAD_STAGE2_VALIDATE_ONLY=1 builds accounting without writing any derivatives.
# The summary wrapper sets DAD_STAGE2_MODE=summary and delegates here.

#SBATCH --job-name=DAD_meg_stage2
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --export=ALL
#SBATCH --array=4-39

set -euo pipefail

# Resolve the operation before environment setup, with a strict cohort guard.
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
MODE="${DAD_STAGE2_MODE:-subject}"
ARGS=()
case "$MODE" in
    subject)
        TASK_ID="${SLURM_ARRAY_TASK_ID:?Submit Stage 2 as a SLURM array.}"
        if (( TASK_ID < 4 || TASK_ID > 39 )); then
            echo "ERROR: Stage 2 participant must be 4-39; got $TASK_ID" >&2
            exit 2
        fi
        SUBJECT_ID="$(printf '%02d' "$TASK_ID")"
        ARGS+=(--subject-config "$PROJECT_ROOT/meg/params/sub-${SUBJECT_ID}.yaml")
        ;;
    summary) TASK_ID="summary"; ARGS+=(--summarize-cohort) ;;
    *) echo "ERROR: unknown DAD_STAGE2_MODE: $MODE" >&2; exit 2 ;;
esac
for control in DAD_STAGE2_OVERWRITE DAD_STAGE2_VALIDATE_ONLY; do
    value="${!control:-0}"
    [[ "$value" == 0 || "$value" == 1 ]] || { echo "ERROR: $control must be 0 or 1" >&2; exit 2; }
done
[[ "${DAD_STAGE2_OVERWRITE:-0}" == 1 ]] && ARGS+=(--overwrite)
[[ "${DAD_STAGE2_VALIDATE_ONLY:-0}" == 1 ]] && ARGS+=(--validate-only)

# Activate the existing environment; isolate caches from shared input folders.
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
[[ -x "$MAMBA_BIN" ]] || { echo "ERROR: micromamba missing: $MAMBA_BIN" >&2; exit 2; }
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "${DAD_MAMBA_ENV:-dad-meg}"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
JOB_TMP="${TMPDIR:-/tmp}/dad_stage2_${SLURM_JOB_ID:-manual}_${TASK_ID}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib" XDG_CACHE_HOME="$JOB_TMP/cache"

# The standalone entry point leaves all Stage 1 provenance and tables intact.
cd "$PROJECT_ROOT"
python "$PROJECT_ROOT/meg/code/preprocessing_basic/trial_preparation.py" \
    --config "$PROJECT_ROOT/meg/params/pipeline_default.yaml" \
    --site-config "$PROJECT_ROOT/meg/params/pipeline_cluster.yaml" "${ARGS[@]}"
