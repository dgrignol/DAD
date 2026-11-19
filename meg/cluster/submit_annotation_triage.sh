#!/bin/bash
#
# Import completed component choices and count preliminary annotation retention.
#
# Usage from the cluster DAD root, after deployment (ONE job, not an array):
#   DAD_COMPONENT_REPORT_LABEL=ocular-reports-20261030-v1 \
#   DAD_TRIAGE_DECISIONS=/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives/annotation_triage/review_notes/ocular_component_choices_completed_20261103.csv \
#   DAD_TRIAGE_LABEL=retention-20261103-v1 \
#     sbatch meg/cluster/submit_annotation_triage.sh
#
# Inputs: completed CSV choices, exact offline report manifests/HTML, compatible
# preliminary ICA caches, current source/muscle/manual evidence and accepted
# Stage 2/3 timing/eligibility. All CSV rows are validated before the importer
# writes any choices. Existing compatible choices are retained byte-for-byte;
# conflicting choices or stale inputs stop the job before retention counting.
#
# Data flow: validate/import the complete decision table -> run annotation_triage
# for P04-P39 -> immutable masks plus a cohort retention/ranking snapshot. Earlier
# component selections also contribute; current manual work always takes priority.
# Outputs stay under derivatives/annotation_triage. No production annotations,
# manual reviews, tuning, acceptance records, final ICA or epochs are created.
# Counts are conditional annotation-only estimates, not final clean epoch counts.
# EDF/SQUID candidates are not generated. This wrapper submits no other job.
#
# Required: DAD_COMPONENT_REPORT_LABEL, DAD_TRIAGE_DECISIONS, DAD_TRIAGE_LABEL.
# Independent participants run in parallel (up to 36 by default). Whole-table
# validation/recheck phases finish before parallel per-participant choice saves;
# cohort publication follows completion of all counting workers. Each child uses
# one numerical-library thread. --cpus-per-task can change the resource budget.
# Optional: DAD_PARTICIPANT_WORKERS (default allocated CPUs), DAD_PROJECT_ROOT,
# DAD_DERIVATIVES_ROOT, MICROMAMBA_BIN, DAD_MAMBA_ENV.
# After a failed job, compatible choices/masks can be reused; select a NEW output
# label if the previous attempt already published a report snapshot.
#
#SBATCH --job-name=DAD_annotation_triage
#SBATCH --output=cluster_logs/%x.%j.log
#SBATCH --error=cluster_logs/%x.%j.err
#SBATCH --cpus-per-task=36
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL

set -euo pipefail

# Resolve this exact decision/report batch; reject accidental duplicate arrays.
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
REPORT_LABEL="${DAD_COMPONENT_REPORT_LABEL:?Set DAD_COMPONENT_REPORT_LABEL to the reviewed report batch.}"
DECISIONS="${DAD_TRIAGE_DECISIONS:?Set DAD_TRIAGE_DECISIONS to the completed CSV path.}"
TRIAGE_LABEL="${DAD_TRIAGE_LABEL:?Set DAD_TRIAGE_LABEL to a new retention snapshot name.}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
PARTICIPANT_WORKERS="${DAD_PARTICIPANT_WORKERS:-${SLURM_CPUS_PER_TASK:-36}}"
if [[ ! "$PARTICIPANT_WORKERS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: DAD_PARTICIPANT_WORKERS must be a positive integer." >&2
    exit 2
fi
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: submit ONE cohort job without --array." >&2
    exit 2
fi
for LABEL in "$REPORT_LABEL" "$TRIAGE_LABEL"; do
    if [[ ! "$LABEL" =~ ^[A-Za-z0-9_.-]+$ || "$LABEL" == "." || "$LABEL" == ".." ]]; then
        echo "ERROR: report labels must be simple folder names." >&2
        exit 2
    fi
done
REPORTS_DIR="$DERIVATIVES_ROOT/annotation_triage/component_reports/$REPORT_LABEL"
OUTPUT_DIR="$DERIVATIVES_ROOT/annotation_triage/reports/$TRIAGE_LABEL"
if [[ ! -f "$DECISIONS" || ! -d "$REPORTS_DIR" ]]; then
    echo "ERROR: missing completed CSV or reviewed component-report folder." >&2
    exit 2
fi
if [[ -e "$OUTPUT_DIR" ]]; then
    echo "ERROR: retention snapshot already exists; choose a new DAD_TRIAGE_LABEL: $OUTPUT_DIR" >&2
    exit 2
fi
if [[ ! -x "$MAMBA_BIN" ]]; then
    echo "ERROR: micromamba executable not found: $MAMBA_BIN" >&2
    exit 2
fi

# Activate the established compute environment; keep all plotting noninteractive.
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "$MAMBA_ENV"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg
# Parallelism is across participants, not nested BLAS/OpenMP thread teams.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
TRIAGE_TMP="${TMPDIR:-/tmp}/dad_annotation_triage_${SLURM_JOB_ID:-manual}"
mkdir -p "$TRIAGE_TMP/matplotlib"
export MPLCONFIGDIR="$TRIAGE_TMP/matplotlib"
cd "$PROJECT_ROOT"
CONFIG_ARGS=(
    --config "$PROJECT_ROOT/meg/params/pipeline_default.yaml"
    --site-config "$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
    --subject-config-dir "$PROJECT_ROOT/meg/params"
)

# Import uses whole-table validation internally. A failure stops this shell, so
# retention cannot silently continue with a missing or partially imported table.
echo ">>> Validating and importing completed component choices: $DECISIONS"
python -B "$PROJECT_ROOT/meg/code/preprocessing_basic/import_ocular_component_choices.py" \
    "${CONFIG_ARGS[@]}" --reports-dir "$REPORTS_DIR" --decisions "$DECISIONS" --reviewer Damiano --apply --jobs "$PARTICIPANT_WORKERS"

# Count all actual acquisitions, preserving manual reviews and earlier selections.
# No production pipeline runner is invoked, and no acceptance gate is modified.
echo ">>> Generating provisional BAD masks and cohort retention: $TRIAGE_LABEL"
python -B "$PROJECT_ROOT/meg/code/preprocessing_basic/annotation_triage.py" \
    "${CONFIG_ARGS[@]}" --subjects {4..39} --label "$TRIAGE_LABEL" --jobs "$PARTICIPANT_WORKERS"
echo ">>> Preliminary retention report ready: $OUTPUT_DIR/report.md"
