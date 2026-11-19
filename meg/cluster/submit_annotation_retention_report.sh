#!/bin/bash
#
# Render embedded HTML retention books, optionally recounting revised annotations.
#
# First report from the running count job (submit ONE dependent job, no array):
#   DAD_TRIAGE_LABEL=retention-20261103-v1 \
#     sbatch --dependency=afterok:10590 meg/cluster/submit_annotation_retention_report.sh
#
# After saving targeted annotation reviews, recompute the full cohort and compare:
#   DAD_TRIAGE_LABEL=retention-20261103-v2 DAD_TRIAGE_RECOUNT=1 \
#   DAD_TRIAGE_COMPARE_TO=retention-20261103-v1 \
#     sbatch meg/cluster/submit_annotation_retention_report.sh
#
# Data flow: optional existing annotation_triage.py recomputation -> validated
# snapshot tables -> render_annotation_triage_report.py -> one cohort HTML and
# one HTML per participant, with embedded scientific figures. Render-only mode
# reads no MEG/ICA samples. Recount mode uses current manual work, existing choices
# and compatible masks; it does not repeat the completed notes import. Old count
# and HTML snapshots remain immutable. Counts precede final ICA/amplitude rejection.
#
# Required: DAD_TRIAGE_LABEL (input snapshot, or NEW snapshot when recounting).
# Optional: DAD_TRIAGE_RECOUNT=0/1 (default 0), DAD_TRIAGE_COMPARE_TO (old snapshot
# label), DAD_TRIAGE_VISUAL_LABEL (default same as DAD_TRIAGE_LABEL; select a new
# one to rerender existing counts), DAD_PROJECT_ROOT, DAD_DERIVATIVES_ROOT,
# MICROMAMBA_BIN, DAD_MAMBA_ENV, DAD_PARTICIPANT_WORKERS (default allocated CPUs).
# Independent participant recounts and books run in parallel (up to all 36 by
# default); cohort summaries/publication wait for every participant. Numerical
# libraries use one thread per child. Path overrides must match pipeline_cluster.yaml.
# No notes import, production cleaning/acceptance, final ICA, epoching, edits to
# manual annotations, dependency installation or nested job submission occurs.
#
#SBATCH --job-name=DAD_retention_report
#SBATCH --output=cluster_logs/%x.%j.log
#SBATCH --error=cluster_logs/%x.%j.err
#SBATCH --cpus-per-task=36
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL

set -euo pipefail

# Resolve independent immutable count and visual labels; avoid duplicate arrays.
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
LABEL="${DAD_TRIAGE_LABEL:?Set DAD_TRIAGE_LABEL to the count snapshot name.}"
VISUAL_LABEL="${DAD_TRIAGE_VISUAL_LABEL:-$LABEL}"
COMPARE_TO="${DAD_TRIAGE_COMPARE_TO:-}"
RECOUNT="${DAD_TRIAGE_RECOUNT:-0}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
PARTICIPANT_WORKERS="${DAD_PARTICIPANT_WORKERS:-${SLURM_CPUS_PER_TASK:-36}}"
if [[ ! "$PARTICIPANT_WORKERS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: DAD_PARTICIPANT_WORKERS must be a positive integer." >&2
    exit 2
fi
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" || ! "$RECOUNT" =~ ^[01]$ ]]; then
    echo "ERROR: use one job without --array; DAD_TRIAGE_RECOUNT must be 0 or 1." >&2
    exit 2
fi
for VALUE in "$LABEL" "$VISUAL_LABEL" ${COMPARE_TO:+"$COMPARE_TO"}; do
    if [[ ! "$VALUE" =~ ^[A-Za-z0-9_.-]+$ || "$VALUE" == "." || "$VALUE" == ".." ]]; then
        echo "ERROR: snapshot labels must be simple folder names." >&2
        exit 2
    fi
done
TRIAGE_ROOT="$DERIVATIVES_ROOT/annotation_triage"
REPORT_DIR="$TRIAGE_ROOT/reports/$LABEL"
if [[ -e "$TRIAGE_ROOT/visual_reports/$VISUAL_LABEL" ]]; then
    echo "ERROR: HTML snapshot already exists; choose a new DAD_TRIAGE_VISUAL_LABEL." >&2
    exit 2
fi
if [[ "$RECOUNT" == "1" && -e "$REPORT_DIR" ]]; then
    echo "ERROR: count snapshot already exists; choose a new DAD_TRIAGE_LABEL." >&2
    exit 2
fi
if [[ "$RECOUNT" == "0" && ! -f "$REPORT_DIR/report.json" ]]; then
    echo "ERROR: count snapshot not ready; wait for counting to finish or use an afterok dependency." >&2
    exit 2
fi
if [[ -n "$COMPARE_TO" && ! -f "$TRIAGE_ROOT/reports/$COMPARE_TO/report.json" ]]; then
    echo "ERROR: comparison snapshot not found." >&2
    exit 2
fi
if [[ ! -x "$MAMBA_BIN" ]]; then
    echo "ERROR: micromamba executable not found: $MAMBA_BIN" >&2
    exit 2
fi

# Use the existing environment and independent noninteractive plot cache.
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "$MAMBA_ENV"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
REPORT_TMP="${TMPDIR:-/tmp}/dad_retention_report_${SLURM_JOB_ID:-manual}"
mkdir -p "$REPORT_TMP/matplotlib"
export MPLCONFIGDIR="$REPORT_TMP/matplotlib"
cd "$PROJECT_ROOT"

# Recount only when requested; a count failure stops HTML publication.
if [[ "$RECOUNT" == "1" ]]; then
    python -B "$PROJECT_ROOT/meg/code/preprocessing_basic/annotation_triage.py" \
        --config "$PROJECT_ROOT/meg/params/pipeline_default.yaml" \
        --site-config "$PROJECT_ROOT/meg/params/pipeline_cluster.yaml" \
        --subject-config-dir "$PROJECT_ROOT/meg/params" --subjects {4..39} --label "$LABEL" --jobs "$PARTICIPANT_WORKERS"
fi

# Rendering consumes only saved tables and embeds all PNG bytes in each HTML.
COMMAND=(python -B "$PROJECT_ROOT/meg/code/preprocessing_basic/render_annotation_triage_report.py"
    --report-dir "$REPORT_DIR" --label "$VISUAL_LABEL" --jobs "$PARTICIPANT_WORKERS")
if [[ -n "$COMPARE_TO" ]]; then
    COMMAND+=(--compare-to "$TRIAGE_ROOT/reports/$COMPARE_TO")
fi
"${COMMAND[@]}"
