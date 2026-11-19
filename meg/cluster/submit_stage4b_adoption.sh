#!/bin/bash
# Adopt/reconcile current BAD masks under the explicit 2026-11-07 policy.
#
# Usage from the cluster DAD root, after a validated immutable plan exists:
#   DAD_ADOPTION_LABEL=current-masks-20261109-v1 \
#     sbatch --array=4-14,16-39 meg/cluster/submit_stage4b_adoption.sh
#   DAD_ADOPTION_LABEL=current-masks-20261109-v1 DAD_ADOPTION_MODE=reconcile \
#     sbatch --cpus-per-task=35 --mem=64G meg/cluster/submit_stage4b_adoption.sh
# For a fresh plan: DAD_ADOPTION_MODE=plan with the same cohort allocation.
#
# Data flow: plan (whole-cohort validation) -> disjoint participant apply array ->
# independent cohort reconciliation -> policy processing acceptance, joined before
# publication. Reconcile must be ONE job, without --array. Default apply requires
# --array explicitly so a missing selector cannot process the wrong participant.
# DAD_ADOPTION_ACCEPT=1 (default) records the user's authorized processing policy
# only after every participant audit succeeds; 0 performs a read-only audit.
# This never submits jobs, invents manual receipts, fits ICA or creates epochs.
# P15 remains in storage and is excluded from every operation in this wrapper.
#
# Resources: one CPU/8 GB per apply task; cohort plan/reconcile uses the requested
# CPUs with one numerical thread per participant process. Override --mem and
# --cpus-per-task at submission; DAD_PARTICIPANT_WORKERS may reduce workers.
# Other controls: DAD_PROJECT_ROOT, DAD_DERIVATIVES_ROOT, DAD_MAMBA_ENV,
# MICROMAMBA_BIN, DAD_ADOPTION_SNAPSHOT. Per-run progress is printed unbuffered.
# Old working annotations/reports are archived by rename/copy before replacement;
# all earlier triage snapshots and manual files/receipts remain immutable.
#SBATCH --job-name=DAD_stage4b_adopt
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL

set -euo pipefail

# Resolve one authorized phase and the shared read-only scientific inputs.
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
LABEL="${DAD_ADOPTION_LABEL:?Set DAD_ADOPTION_LABEL to the immutable adoption plan label.}"
MODE="${DAD_ADOPTION_MODE:-apply}"
ACCEPT="${DAD_ADOPTION_ACCEPT:-1}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
if [[ ! "$LABEL" =~ ^[A-Za-z0-9_.-]+$ || "$LABEL" == . || "$LABEL" == .. ]]; then
    echo 'ERROR: invalid adoption label.' >&2; exit 2
fi
if [[ "$ACCEPT" != 0 && "$ACCEPT" != 1 ]]; then
    echo 'ERROR: DAD_ADOPTION_ACCEPT must be 0 or 1.' >&2; exit 2
fi
case "$MODE" in
    apply)
        TASK="${SLURM_ARRAY_TASK_ID:?Submit apply as --array=4-14,16-39.}"
        if [[ ! "$TASK" =~ ^[0-9]+$ ]] || (( TASK < 4 || TASK > 39 || TASK == 15 )); then
            echo 'ERROR: apply requires a working-cohort participant (P15 excluded).' >&2; exit 2
        fi
        ;;
    plan|reconcile)
        if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
            echo 'ERROR: submit ONE cohort plan/reconcile job without --array.' >&2; exit 2
        fi
        ;;
    *) echo 'ERROR: mode must be plan, apply or reconcile.' >&2; exit 2 ;;
esac

# Activate existing runtime and bound nested numerical work explicitly.
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "${DAD_MAMBA_ENV:-dad-meg}"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export MPLBACKEND=Agg MNE_BROWSER_BACKEND=matplotlib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ADOPTION_TMP="${TMPDIR:-/tmp}/dad_adoption_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-cohort}"
mkdir -p "$ADOPTION_TMP/matplotlib" "$ADOPTION_TMP/cache"
export MPLCONFIGDIR="$ADOPTION_TMP/matplotlib" XDG_CACHE_HOME="$ADOPTION_TMP/cache"
cd "$PROJECT_ROOT"
mkdir -p cluster_logs
COMMAND=(python -B "$PROJECT_ROOT/meg/code/preprocessing_basic/stage4b_adoption.py" "$MODE"
    --config "$PROJECT_ROOT/meg/params/pipeline_default.yaml"
    --site-config "$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
    --subject-config-dir "$PROJECT_ROOT/meg/params" --label "$LABEL")
if [[ "$MODE" == apply ]]; then
    COMMAND+=(--subject "$TASK")
else
    COMMAND+=(--jobs "${DAD_PARTICIPANT_WORKERS:-${SLURM_CPUS_PER_TASK:-1}}")
    if [[ "$MODE" == plan ]]; then
        COMMAND+=(--snapshot "${DAD_ADOPTION_SNAPSHOT:-retention-20261107-v1}")
    elif [[ "$ACCEPT" == 1 ]]; then
        COMMAND+=(--accept-policy)
    fi
fi
echo ">>> Stage 4B policy phase=$MODE label=$LABEL CPUs=${SLURM_CPUS_PER_TASK:-1} participant=${TASK:-cohort}"
"${COMMAND[@]}"
echo ">>> Finished $MODE; final component decisions and cleaned-output acceptance remain open."
