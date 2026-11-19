#!/bin/bash
#
# Build ONE offline ocular-component HTML report per participant.
#
# Usage from the cluster DAD root, after deployment:
#   DAD_COMPONENT_REPORT_LABEL=ocular-reports-20261030 sbatch --array=15-39 meg/cluster/submit_ocular_component_reports.sh
# For a new full-cohort batch choose a new label and --array=4-39.
# Reuse the same label for guarded resume of compatible report generation.
#
# Inputs: accepted 4A/current provisional muscle evidence, existing cached
# preliminary ICA and saved choices/manual reviews. Only pending runs get
# display traces. No new ICA fitting, artifact detection, annotations,
# component-choice persistence, scientific acceptance or epoching is performed.
# Outputs: component_reports/LABEL beneath the derivative annotation_triage
# namespace. Each participant has one offline HTML with internal run navigation,
# a source/model-bound manifest and a blank decision template. Internal per-run
# assets support interrupted generation; the cohort index updates under a lock.
# Existing scientific files remain read-only. This wrapper submits no other job.
#
# Optional overrides: DAD_PROJECT_ROOT, MICROMAMBA_BIN, DAD_MAMBA_ENV.
# The label is mandatory so every array task writes to the same batch folder.
#
#SBATCH --job-name=DAD_ocular_reports
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL
#SBATCH --array=4-39

set -euo pipefail

# Resolve exact cohort members and report-only runtime controls.
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
LABEL="${DAD_COMPONENT_REPORT_LABEL:?Set DAD_COMPONENT_REPORT_LABEL to the shared report batch name.}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?Submit this wrapper as a SLURM array.}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
if [[ ! "$TASK_ID" =~ ^[0-9]+$ ]] || (( TASK_ID < 4 || TASK_ID > 39 )); then
    echo "ERROR: participant must be P04-P39." >&2
    exit 2
fi
if [[ ! "$LABEL" =~ ^[A-Za-z0-9_.-]+$ || "$LABEL" == "." || "$LABEL" == ".." ]]; then
    echo "ERROR: report label must be a simple folder name." >&2
    exit 2
fi
if [[ ! -x "$MAMBA_BIN" ]]; then
    echo "ERROR: micromamba executable not found: $MAMBA_BIN" >&2
    exit 2
fi

# Activate the existing environment; do not install dependencies or open GUI.
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "$MAMBA_ENV"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
REPORT_TMP="${TMPDIR:-/tmp}/dad_ocular_reports_${SLURM_JOB_ID:-manual}_$TASK_ID"
mkdir -p "$REPORT_TMP/matplotlib"
export MPLCONFIGDIR="$REPORT_TMP/matplotlib"

# The builder validates source/model lineage and writes only dedicated reports.
cd "$PROJECT_ROOT"
COMMAND=(
    python -B "$PROJECT_ROOT/meg/code/preprocessing_basic/build_ocular_component_reports.py"
    --config "$PROJECT_ROOT/meg/params/pipeline_default.yaml"
    --site-config "$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
    --subject-config-dir "$PROJECT_ROOT/meg/params"
    --subjects "$TASK_ID" --label "$LABEL" --resume
)
"${COMMAND[@]}"
