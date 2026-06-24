#!/bin/bash
#
# submit_dot_side_group_report.sh
#
# Purpose
# Build the combined group-level dot-side decoding sanity-check report on the
# CIMeC cluster after the per-subject dot-side array jobs have completed. This
# is a one-task aggregation job; it does not rerun participant decoding. The
# report places all group analyses before participant-specific figures.
#
# Submission examples (run from the cluster head node)
#   cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
#   mkdir -p cluster_logs
#   qsub meg/cluster/submit_dot_side_group_report.sh
#
# Use SEM instead of the default 95% CI display band:
#   qsub -v DAD_DOT_SIDE_INTERVAL=sem \
#     meg/cluster/submit_dot_side_group_report.sh
#
# Validate inputs without writing outputs:
#   qsub -v DAD_VALIDATE_ONLY=1 \
#     meg/cluster/submit_dot_side_group_report.sh
#
# Inputs
# - Completed subject-level dot-side summary JSON, scores CSV, and PNG outputs
#   for every requested subject.
#
# Outputs
# - One group PNG per analysis under derivatives/sanity_checks/group/ses-01/
#   dot_side_decoding/<condition>/<model-source>.
# - One three-condition trial-count PNG per participant in the combined report
#   directory, using source counts by identity before pseudo-trial averaging.
# - Combined JSON, linked HTML, and embedded HTML under derivatives/reports/
#   group/ses-01/dot_side_decoding/combined.
#
# Environment variables
# - DAD_PROJECT_ROOT: cluster DAD root.
# - DAD_MAMBA_ENV, MICROMAMBA_BIN: Python environment selection.
# - DAD_DOT_SIDE_ANALYSES: comma-separated analysis keys, default all. Valid
#   keys are always_visible_observed, occluded_nondeviant_observed,
#   occluded_deviant_observed, and occluded_deviant_predicted.
# - DAD_DOT_SIDE_GROUP_SUBJECTS: comma-separated subjects, default 01..07.
# - DAD_DOT_SIDE_METRIC: scores CSV metric column, default balanced accuracy.
# - DAD_DOT_SIDE_INTERVAL: ci95 or sem; default ci95.
# - DAD_CLUSTER_FORMING_P: pointwise p for cluster formation; default 0.05.
# - DAD_CLUSTER_ALPHA: cluster-level significance alpha; default 0.05.
# - DAD_CLUSTER_TAIL: greater or two-sided; default greater.
# - DAD_N_PERMUTATIONS: 0 means exact sign flips when feasible; default 0.
# - DAD_RANDOM_SEED: randomized permutation seed; default 42.
# - DAD_VALIDATE_ONLY: 1 validates inputs and exits before writing outputs.
# - DAD_ALLOW_MISSING: 1 builds from available subjects and lists missing ones.
# - DAD_REBUILD: 1 overwrites existing group outputs; 0 resumes existing outputs.

#$ -N DAD_dot_side_group
#$ -S /bin/bash
#$ -cwd
#$ -o cluster_logs/$JOB_NAME.$JOB_ID.log
#$ -j y
#$ -pe smp 2
#$ -l h_vmem=8G

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve paths and analysis settings
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
REBUILD="${DAD_REBUILD:-0}"
VALIDATE_ONLY="${DAD_VALIDATE_ONLY:-0}"
ALLOW_MISSING="${DAD_ALLOW_MISSING:-0}"
ANALYSES_CSV="${DAD_DOT_SIDE_ANALYSES:-all}"
SUBJECTS_CSV="${DAD_DOT_SIDE_GROUP_SUBJECTS:-01,02,03,04,05,06,07}"
METRIC="${DAD_DOT_SIDE_METRIC:-accuracy_balanced_undersampled_sets}"
INTERVAL="${DAD_DOT_SIDE_INTERVAL:-ci95}"
CLUSTER_FORMING_P="${DAD_CLUSTER_FORMING_P:-0.05}"
CLUSTER_ALPHA="${DAD_CLUSTER_ALPHA:-0.05}"
CLUSTER_TAIL="${DAD_CLUSTER_TAIL:-greater}"
N_PERMUTATIONS="${DAD_N_PERMUTATIONS:-0}"
RANDOM_SEED="${DAD_RANDOM_SEED:-42}"

if [[ "$REBUILD" != "0" && "$REBUILD" != "1" ]]; then
    echo "ERROR: DAD_REBUILD must be 0 or 1; got ${REBUILD}." >&2
    exit 2
fi
if [[ "$VALIDATE_ONLY" != "0" && "$VALIDATE_ONLY" != "1" ]]; then
    echo "ERROR: DAD_VALIDATE_ONLY must be 0 or 1; got ${VALIDATE_ONLY}." >&2
    exit 2
fi
if [[ "$ALLOW_MISSING" != "0" && "$ALLOW_MISSING" != "1" ]]; then
    echo "ERROR: DAD_ALLOW_MISSING must be 0 or 1; got ${ALLOW_MISSING}." >&2
    exit 2
fi

DEFAULT_CONFIG="$PROJECT_ROOT/meg/params/pipeline_default.yaml"
SITE_CONFIG="$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
SUBJECT_CONFIG="$PROJECT_ROOT/meg/params/sub-01.yaml"
SCRIPT="$PROJECT_ROOT/meg/code/sanity_checks/02_dot_side_group_report.py"
IFS=',' read -r -a SUBJECTS <<< "$SUBJECTS_CSV"
IFS=',' read -r -a ANALYSES <<< "$ANALYSES_CSV"


# -----------------------------------------------------------------------------
# Activate the shared headless Python environment
# -----------------------------------------------------------------------------

if [[ ! -x "$MAMBA_BIN" ]]; then
    echo "ERROR: micromamba executable not found at $MAMBA_BIN." >&2
    exit 2
fi

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "$MAMBA_ENV"

export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS="${NSLOTS:-2}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-2}"
export MKL_NUM_THREADS="${NSLOTS:-2}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-2}"

JOB_TMP="${TMPDIR:-/tmp}/dad_dot_side_group_${JOB_ID:-manual}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib"
export XDG_CACHE_HOME="$JOB_TMP/cache"


# -----------------------------------------------------------------------------
# Validate deployment and build the group report
# -----------------------------------------------------------------------------

cd "$PROJECT_ROOT"
mkdir -p cluster_logs

for required_path in \
    "$DEFAULT_CONFIG" \
    "$SITE_CONFIG" \
    "$SUBJECT_CONFIG" \
    "$SCRIPT" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/drsa_result_utils.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/pipeline_config.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/preprocessing_stage_utils.py"; do
    if [[ ! -f "$required_path" ]]; then
        echo "ERROR: required file is missing: $required_path" >&2
        exit 2
    fi
done

COMMAND=(
    python
    "$SCRIPT"
    --config "$DEFAULT_CONFIG"
    --site-config "$SITE_CONFIG"
    --subject-config "$SUBJECT_CONFIG"
    --analyses "${ANALYSES[@]}"
    --group-subjects "${SUBJECTS[@]}"
    --metric "$METRIC"
    --interval "$INTERVAL"
    --cluster-forming-p "$CLUSTER_FORMING_P"
    --cluster-alpha "$CLUSTER_ALPHA"
    --tail "$CLUSTER_TAIL"
    --n-permutations "$N_PERMUTATIONS"
    --random-seed "$RANDOM_SEED"
    --print-summary
)

if [[ "$ALLOW_MISSING" == "1" ]]; then
    COMMAND+=(--allow-missing)
fi

if [[ "$VALIDATE_ONLY" == "1" ]]; then
    COMMAND+=(--validate-only)
    echo ">>> Validate-only mode: no group report outputs will be written."
elif [[ "$REBUILD" == "1" ]]; then
    COMMAND+=(--overwrite)
    echo ">>> Rebuild mode: group dot-side report outputs may be replaced."
else
    COMMAND+=(--resume)
    echo ">>> Resume mode: completed group dot-side report outputs are reused."
fi

echo ">>> Analyses: ${ANALYSES[*]}"
echo ">>> Subjects: ${SUBJECTS[*]}"
echo ">>> Metric: $METRIC"
echo ">>> Interval: $INTERVAL"
echo ">>> Cluster-forming p: $CLUSTER_FORMING_P"
echo ">>> Cluster alpha: $CLUSTER_ALPHA"
echo ">>> Tail: $CLUSTER_TAIL"
echo ">>> Host: $(hostname)"
echo ">>> Slots: ${NSLOTS:-unknown}"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '>>> Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ">>> Embedded report: $DERIVATIVES_ROOT/reports/group/ses-01/dot_side_decoding/combined/dot_side_all_analyses_group_report_embedded.html"
