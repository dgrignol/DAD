#!/bin/bash
#
# submit_preprocessing_review.sh
#
# Purpose
# Run the first, noninteractive DAD preprocessing review pass as one Sun Grid
# Engine array task per participant. Each task runs preflight, Maxwell
# filtering, physiological artifact annotation, and manual-ICA diagnostics.
# It creates the reports and plots needed to choose bad-channel policy,
# per-run muscle thresholds, and manual ICA exclusions. It never runs the
# deprecated automatic ICA branch and never applies ICA rejection.
#
# Submission examples (run from the DAD repository on the cluster head node)
#   mkdir -p cluster_logs
#   qsub -t 2-7 meg/cluster/submit_preprocessing_review.sh
#
# Include sub01 as well:
#   qsub -t 1-7 meg/cluster/submit_preprocessing_review.sh
#
# Select a different micromamba environment:
#   qsub -v DAD_MAMBA_ENV=dad-meg -t 2-7 \
#     meg/cluster/submit_preprocessing_review.sh
#
# Resume one failed participant:
#   qsub -t 4 meg/cluster/submit_preprocessing_review.sh
#
# Resume only selected stages after an infrastructure failure:
#   qsub -v DAD_STAGES=artifact_annotation:ica_manual -t 1 \
#     meg/cluster/submit_preprocessing_review.sh
#   qsub -v DAD_STAGES=ica_manual -t 2-7 \
#     meg/cluster/submit_preprocessing_review.sh
#
# Intentionally rebuild stage-03 derivatives after a scientific/configuration
# change. This switches from guarded resume to explicit overwrite:
#   qsub -v DAD_REBUILD=1 -t 1-7 \
#     meg/cluster/submit_preprocessing_review.sh
#
# ``submit_preprocessing_apply.sh`` sets the apply-specific variables and
# should be used instead of setting DAD_PROCESSING_PASS manually.
#
# Inputs
# - Read-only raw data below
#   /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/data
# - Scientific defaults in meg/params/pipeline_default.yaml
# - Cluster paths in meg/params/pipeline_cluster.yaml
# - Subject selector meg/params/sub-XX.yaml
# - Existing 03_tuning_choices.yaml only when that subject was already reviewed
#
# Outputs
# - Guarded derivatives below
#   /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives
# - Per-subject JSON/HTML reports and PNG decision plots
# - A compatible prefit manual ICA model, but no cleaned *_ica.fif files
# - SGE console logs below cluster_logs/
#
# Data flow
# SGE_TASK_ID -> subject selector -> config/preflight -> 03_A review ->
# 03_B review with FLUX threshold 10 unless tuned -> 03_C manual ICA review ->
# reports and diagnostic plots for human inspection.
#
# Safety and restart behavior
# - Raw data are never written.
# - Normal mode uses --resume and validates existing derivatives.
# - DAD_REBUILD=1 uses --overwrite and intentionally recreates selected stages.
# - Existing tuning is reused only when 03_tuning_choices.yaml already exists.
# - ICA manual_exclude=null remains review_required and is never auto-filled.

#$ -N DAD_meg_review
#$ -S /bin/bash
#$ -cwd
#$ -o cluster_logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -j y
#$ -pe smp 4
#$ -l h_vmem=24G
#$ -t 1-7

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve the subject and shared cluster paths
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
TASK_ID="${SGE_TASK_ID:?SGE_TASK_ID is required; submit this as an SGE array job.}"
STAGE_SPEC="${DAD_STAGES:-maxwell_filter:artifact_annotation:ica_manual}"
REBUILD="${DAD_REBUILD:-0}"
PROCESSING_PASS="${DAD_PROCESSING_PASS:-review}"

if (( TASK_ID < 1 || TASK_ID > 7 )); then
    echo "ERROR: SGE_TASK_ID must be between 1 and 7; got ${TASK_ID}." >&2
    exit 2
fi
if [[ "$REBUILD" != "0" && "$REBUILD" != "1" ]]; then
    echo "ERROR: DAD_REBUILD must be 0 or 1; got ${REBUILD}." >&2
    exit 2
fi
if [[ "$PROCESSING_PASS" != "review" && "$PROCESSING_PASS" != "apply" ]]; then
    echo "ERROR: DAD_PROCESSING_PASS must be review or apply." >&2
    exit 2
fi
if [[ "$PROCESSING_PASS" == "apply" && "$REBUILD" != "1" ]]; then
    echo "ERROR: apply requires DAD_REBUILD=1 and explicit overwrite." >&2
    exit 2
fi

# A colon-separated environment variable avoids both shell quoting ambiguity
# and SGE's comma delimiter for separate ``qsub -v`` variable assignments.
IFS=':' read -r -a SELECTED_STAGES <<< "$STAGE_SPEC"
if (( ${#SELECTED_STAGES[@]} == 0 )); then
    echo "ERROR: DAD_STAGES must select at least one stage." >&2
    exit 2
fi
for stage in "${SELECTED_STAGES[@]}"; do
    case "$stage" in
        maxwell_filter|artifact_annotation|ica_manual)
            ;;
        *)
            echo "ERROR: unsupported DAD_STAGES value: $stage" >&2
            exit 2
            ;;
    esac
done

SUBJECT_ID="$(printf "%02d" "$TASK_ID")"
SUBJECT_LABEL="sub-${SUBJECT_ID}"
DEFAULT_CONFIG="$PROJECT_ROOT/meg/params/pipeline_default.yaml"
SITE_CONFIG="$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
SUBJECT_CONFIG="$PROJECT_ROOT/meg/params/${SUBJECT_LABEL}.yaml"
REPORT_DIR="$DERIVATIVES_ROOT/reports/${SUBJECT_LABEL}/ses-01"
TUNING_CONFIG="$REPORT_DIR/03_tuning_choices.yaml"


# -----------------------------------------------------------------------------
# Activate the reproducible Python environment
# -----------------------------------------------------------------------------

if [[ ! -x "$MAMBA_BIN" ]]; then
    echo "ERROR: micromamba executable not found at $MAMBA_BIN." >&2
    echo "Set MICROMAMBA_BIN or install it as documented in CIMeC_cluster_docs.md." >&2
    exit 2
fi

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "$MAMBA_ENV"

export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
# MNE raw/source browsers have a backend independent of MPLBACKEND. Force the
# matplotlib implementation for this noninteractive headless review job.
export MNE_BROWSER_BACKEND=matplotlib
export OMP_NUM_THREADS="${NSLOTS:-4}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-4}"
export MKL_NUM_THREADS="${NSLOTS:-4}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-4}"

# Matplotlib and other libraries must not race on a shared home-directory
# cache when all seven array tasks start together.
JOB_TMP="${TMPDIR:-/tmp}/dad_meg_${JOB_ID:-manual}_${TASK_ID}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib"
export XDG_CACHE_HOME="$JOB_TMP/cache"


# -----------------------------------------------------------------------------
# Validate required paths and runtime packages before expensive processing
# -----------------------------------------------------------------------------

cd "$PROJECT_ROOT"
mkdir -p cluster_logs

for required_path in \
    "$DEFAULT_CONFIG" \
    "$SITE_CONFIG" \
    "$SUBJECT_CONFIG" \
    "$PROJECT_ROOT/meg/data/fwmegmaxfilter/ct_sparse_triux2.fif" \
    "$PROJECT_ROOT/meg/data/fwmegmaxfilter/sss_cal_3182_250410.dat"; do
    if [[ ! -f "$required_path" ]]; then
        echo "ERROR: required configuration file is missing: $required_path" >&2
        exit 2
    fi
done

python - <<'PY'
import matplotlib
import mne
import numpy
import pandas
import sklearn
import yaml
from packaging.version import Version

if Version(mne.__version__) < Version("1.10.0"):
    raise RuntimeError(
        "DAD task FIF head-shape points require MNE >= 1.10.0; "
        f"found {mne.__version__}. Update the dad-meg environment."
    )

print(
    "Runtime imports OK:",
    f"mne={mne.__version__}",
    f"numpy={numpy.__version__}",
    f"sklearn={sklearn.__version__}",
)
PY


# -----------------------------------------------------------------------------
# Build and execute the review command
# -----------------------------------------------------------------------------

COMMAND=(
    python
    "$PROJECT_ROOT/meg/code/preprocessing_basic/run_pipeline.py"
    --config "$DEFAULT_CONFIG"
    --site-config "$SITE_CONFIG"
    --subject-config "$SUBJECT_CONFIG"
    --stages "${SELECTED_STAGES[@]}"
    --processing-pass "$PROCESSING_PASS"
)
if [[ "$REBUILD" == "1" ]]; then
    COMMAND+=(--overwrite)
    echo ">>> Rebuild mode: existing selected-stage derivatives will be replaced."
else
    COMMAND+=(--resume)
    echo ">>> Resume mode: existing derivatives must pass stage validation."
fi

# A previous interactive review may already have accepted decisions. Reusing
# that subject-specific file makes restart behavior deterministic. A new
# participant has no tuning file and therefore uses review defaults.
if [[ -f "$TUNING_CONFIG" ]]; then
    COMMAND+=(--tuning-config "$TUNING_CONFIG")
    echo ">>> Reusing existing tuning: $TUNING_CONFIG"
elif [[ "$PROCESSING_PASS" == "apply" ]]; then
    echo "ERROR: apply tuning file not found: $TUNING_CONFIG" >&2
    exit 2
else
    echo ">>> No tuning file found; using first-pass review defaults."
fi

echo ">>> Subject: $SUBJECT_LABEL"
echo ">>> Processing pass: $PROCESSING_PASS"
echo ">>> Stages: ${SELECTED_STAGES[*]}"
echo ">>> Host: $(hostname)"
echo ">>> Slots: ${NSLOTS:-unknown}"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '>>> Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ">>> Processing reports: $REPORT_DIR"
echo ">>> Decision plots: $DERIVATIVES_ROOT/preprocessing/$SUBJECT_LABEL/ses-01/figures"
