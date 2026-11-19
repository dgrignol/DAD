#!/bin/bash
# Purpose: user-submitted Stage 3 visual timing, P04-P39; no cleaning/epoching.
# Usage from the deployed DAD root:
#   DAD_STAGE3_VERSION=s3-pilot-v001 DAD_STAGE3_PILOT=1 \
#     sbatch --array=8,9,11,24 meg/cluster/submit_preprocessing_stage3.sh
#   DAD_STAGE3_VERSION=s3-v001 sbatch --array=4-39 \
#     meg/cluster/submit_preprocessing_stage3.sh
# Inputs: accepted Stage 1/2, raw FIFs (read-only), saved Conf, candidate method.
# Outputs: derivatives/visual_timing/VERSION/sub-XX/ses-01 immutable proposals.
# Data flow: array participant -> standalone timing CLI -> explicit later review.
# DAD_STAGE3_VALIDATE_ONLY=1 checks inputs/headers without signal reads/writes.
# DAD_STAGE3_RESUME=1 verifies matching complete bundles. No overwrite exists.
# Pilot tasks read at most 12 x 250 ms windows each; P08/P09/P11/P24 only.
# Resource requests are conservative starting bounds; inspect pilot accounting
# before full cohort submission. This wrapper never submits another job.

#SBATCH --job-name=DAD_meg_stage3
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --export=ALL
#SBATCH --array=4-39

set -euo pipefail

# Resolve version/mode/cohort before any environment or output work.
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
VERSION="${DAD_STAGE3_VERSION:?Set an explicit immutable Stage 3 version.}"
MODE="${DAD_STAGE3_MODE:-subject}"
[[ "$VERSION" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$ && "$VERSION" != *..* ]] || {
    echo "ERROR: invalid Stage 3 version" >&2; exit 2;
}
ARGS=(--version "$VERSION" --method-config "$PROJECT_ROOT/meg/params/stage3_visual_timing.yaml")
for control in DAD_STAGE3_VALIDATE_ONLY DAD_STAGE3_PILOT DAD_STAGE3_RESUME; do
    value="${!control:-0}"
    [[ "$value" == 0 || "$value" == 1 ]] || { echo "ERROR: $control must be 0 or 1" >&2; exit 2; }
done
case "$MODE" in
    subject)
        TASK_ID="${SLURM_ARRAY_TASK_ID:?Submit Stage 3 as a participant array.}"
        [[ "$TASK_ID" =~ ^[0-9]+$ ]] && (( TASK_ID >= 4 && TASK_ID <= 39 )) || {
            echo "ERROR: Stage 3 cohort is P04-P39" >&2; exit 2;
        }
        SUBJECT_ID="$(printf '%02d' "$TASK_ID")"
        ARGS+=(--subject-config "$PROJECT_ROOT/meg/params/sub-${SUBJECT_ID}.yaml")
        if [[ "${DAD_STAGE3_PILOT:-0}" == 1 ]]; then
            case "$TASK_ID" in 8|9|11|24) ARGS+=(--pilot) ;;
                *) echo "ERROR: bounded pilot is P08/P09/P11/P24 only" >&2; exit 2;; esac
        fi
        [[ "${DAD_STAGE3_RESUME:-0}" == 1 ]] && ARGS+=(--resume)
        ;;
    summary)
        [[ "${DAD_STAGE3_PILOT:-0}" == 0 && "${DAD_STAGE3_RESUME:-0}" == 0 ]] || {
            echo "ERROR: full cohort summary cannot use pilot/resume" >&2; exit 2;
        }
        TASK_ID=summary; ARGS+=(--summarize-cohort)
        ;;
    *) echo "ERROR: unknown DAD_STAGE3_MODE: $MODE" >&2; exit 2 ;;
esac
[[ "${DAD_STAGE3_VALIDATE_ONLY:-0}" == 1 ]] && ARGS+=(--validate-only)

# Use the existing environment with caches isolated from accepted evidence.
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
[[ -x "$MAMBA_BIN" ]] || { echo "ERROR: micromamba missing: $MAMBA_BIN" >&2; exit 2; }
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "${DAD_MAMBA_ENV:-dad-meg}"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}" OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
JOB_TMP="${TMPDIR:-/tmp}/dad_stage3_${SLURM_JOB_ID:-manual}_${TASK_ID}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib" XDG_CACHE_HOME="$JOB_TMP/cache"

# Standalone timing never refreshes the runner, inventories or source manifest.
cd "$PROJECT_ROOT"
python "$PROJECT_ROOT/meg/code/preprocessing_basic/visual_timing.py" \
    --config "$PROJECT_ROOT/meg/params/pipeline_default.yaml" \
    --site-config "$PROJECT_ROOT/meg/params/pipeline_cluster.yaml" "${ARGS[@]}"
