#!/bin/bash
# Rebuild FINAL ICA review figures/HTML from exact existing models; NEVER refit.
#
# Usage from cluster DAD root after deployment:
#   DAD_ICA_REPORT_LABEL=continuous-task-20261110-v1 \
#     sbatch --array=4-14,16-39 meg/cluster/submit_final_ica_reports.sh
# Inputs: exact prefit models/metadata, accepted 4A/4B and unchanged annotated
# inputs. Each task validates one participant and archives its old reports and
# figures before adding default native continuous AND task-locked panels for
# every IC, BAD labels and per-plot commands. Default native mode reuses complete
# validated full-fit properties. DAD_ICA_REPORT_MODE=full also rebuilds fit views.
# All components share one chunked native projection per run, never repeated
# per-IC raw reads. Model files/tuning are read-only.
# No ICA fit, exclusions, production apply, cleaned FIF or Stage 5 epochs occur.
# Existing final-model decisions cause a refusal to replace their bound reports;
# use independent read-only plotting commands for those participants instead.
#
# Resources: independent participant array tasks, default 4 CPUs/24 GB. Library
# threads use this allocation; --cpus-per-task/--mem and an optional array cap
# expose resource control. Within each task run one participant only. Outputs
# and immutable refresh labels are participant-specific; no cohort acceptance or
# scientific summary is published by this wrapper. Join the array before review.
# Other overrides: DAD_PROJECT_ROOT, MICROMAMBA_BIN, DAD_MAMBA_ENV.
#SBATCH --job-name=DAD_final_ica_reports
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL
#SBATCH --array=4-14,16-39

set -euo pipefail

# Resolve one participant and the immutable report-only refresh label.
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
LABEL="${DAD_ICA_REPORT_LABEL:?Set a new DAD_ICA_REPORT_LABEL for this report refresh.}"
MODE="${DAD_ICA_REPORT_MODE:-native}"
case "$MODE" in
    native) REFRESH_FLAGS=(--native-only) ;;
    full) REFRESH_FLAGS=(--regenerate) ;;
    *) echo 'ERROR: DAD_ICA_REPORT_MODE must be native or full.' >&2; exit 2 ;;
esac
TASK="${SLURM_ARRAY_TASK_ID:?Submit as --array=4-14,16-39.}"
if [[ ! "$TASK" =~ ^[0-9]+$ ]] || (( TASK < 4 || TASK > 39 || TASK == 15 )); then
    echo 'ERROR: final report refresh requires a working-cohort participant.' >&2; exit 2
fi
if [[ ! "$LABEL" =~ ^[A-Za-z0-9_-]+$ ]]; then
    echo 'ERROR: report label must contain only letters, numbers, underscores or hyphens.' >&2; exit 2
fi
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "${DAD_MAMBA_ENV:-dad-meg}"

# Isolate plotting/cache state and bound numerical-library threads per task.
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export MPLBACKEND=Agg MNE_BROWSER_BACKEND=matplotlib
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS" MKL_NUM_THREADS="$OMP_NUM_THREADS" NUMEXPR_NUM_THREADS="$OMP_NUM_THREADS"
REPORT_TMP="${TMPDIR:-/tmp}/dad_final_reports_${SLURM_JOB_ID:-manual}_${TASK}"
mkdir -p "$REPORT_TMP/matplotlib" "$REPORT_TMP/cache"
export MPLCONFIGDIR="$REPORT_TMP/matplotlib" XDG_CACHE_HOME="$REPORT_TMP/cache"
cd "$PROJECT_ROOT"
echo ">>> Report-only sub-${TASK}; original final model must remain unchanged."
python -B "$PROJECT_ROOT/meg/code/preprocessing_basic/review_final_ica.py" \
    --config "$PROJECT_ROOT/meg/params/pipeline_default.yaml" \
    --site-config "$PROJECT_ROOT/meg/params/pipeline_cluster.yaml" \
    --subject-config-dir "$PROJECT_ROOT/meg/params" \
    --refresh --subject "$TASK" "${REFRESH_FLAGS[@]}" --label "${LABEL}_sub${TASK}" --jobs 1
echo '>>> Final ICA report refresh complete; component decisions remain explicit.'
