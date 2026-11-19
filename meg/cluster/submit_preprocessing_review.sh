#!/bin/bash
#
# submit_preprocessing_review.sh
#
# Purpose
# Run one explicitly selected Stage 4 review step per participant P04-P39.
# Default: 4A only (maxwell_filter). 4B requires accepted 4A; 4C requires
# accepted 4A and 4B. Jobs never infer scientific acceptance.
#
# Usage from remote DAD root
#   sbatch --array=4-14,16-39 meg/cluster/submit_preprocessing_review.sh
#   DAD_STAGES=artifact_annotation sbatch --array=4-14,16-39 meg/cluster/submit_preprocessing_review.sh
#   DAD_STAGES=ica_manual sbatch --array=4-14,16-39 meg/cluster/submit_preprocessing_review.sh
# Prepare only supplementary 4B blink models/sources, without rerunning muscle:
#   DAD_BLINK_AID_ONLY=1 sbatch --array=4-14,16-39 meg/cluster/submit_preprocessing_review.sh
# Pilot one exact source; omit DAD_SOURCE for all participant acquisitions:
#   DAD_BLINK_AID_ONLY=1 DAD_SOURCE=sub04_run01_220526.fif sbatch --array=4 meg/cluster/submit_preprocessing_review.sh
# DAD_REBUILD=1 replaces only selected aid caches in this mode. No final ICA,
# annotations or acceptance receipts are written by aid preparation.
# Plot-only refresh of saved 4A scores and one 20-second raw window per run:
#   DAD_MAXWELL_QC_ONLY=1 sbatch --array=4-14,16-39 meg/cluster/submit_preprocessing_review.sh
# Prepare SQUID candidates from accepted channels, without SSS or 4B writes:
#   DAD_SQUID_PREPARE_ONLY=1 sbatch --array=4 meg/cluster/submit_preprocessing_review.sh
# Optional DAD_SOURCE selects one exact acquisition. Review locally afterward.
# PSDs additionally read matching first-up-to-120-second raw/SSS windows.
# This writes channel PNGs, plot metadata and HTML only; it never invokes the
# pipeline runner, detection, SSS, annotations, ICA or acceptance machinery.
# DAD_QC_WINDOW_SEC defaults to 20; DAD_QC_SEED defaults to 97. Exact native
# sample windows and the saved per-run channel decisions are recorded in plots.
# Rebuild only selected affected review outputs with DAD_REBUILD=1. Ordinary
# runs use guarded resume. Independent review and review_cleaning.py decisions
# precede downstream steps. Raw data are read-only. No epoching/export runs.
# Tuning is read from reports/sub-XX/ses-01/03_tuning_choices.yaml; the new
# calibration location is meg/meg_calibration_files/fwmegmaxfilter.
# The cluster environment/backend/thread/cache setup below is shared with apply.
#
#SBATCH --job-name=DAD_meg_review
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --export=ALL
#SBATCH --array=4-14,16-39

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve the subject and shared cluster paths
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required; submit this as an SLURM array job.}"
STAGE_SPEC="${DAD_STAGES:-maxwell_filter}"
REBUILD="${DAD_REBUILD:-0}"
PROCESSING_PASS="${DAD_PROCESSING_PASS:-review}"
MAXWELL_QC_ONLY="${DAD_MAXWELL_QC_ONLY:-0}"
BLINK_AID_ONLY="${DAD_BLINK_AID_ONLY:-0}"
SQUID_PREPARE_ONLY="${DAD_SQUID_PREPARE_ONLY:-0}"
RAW_ROOT="${DAD_RAW_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/data}"
if [[ "$BLINK_AID_ONLY" == "1" && -z "${DAD_STAGES+x}" ]]; then
    STAGE_SPEC=artifact_annotation
fi

if (( TASK_ID < 4 || TASK_ID > 39 )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID must be between 4 and 39; got ${TASK_ID}." >&2
    exit 2
fi
# The working cohort excludes P15 by explicit 2026-11-07 participant decision.
if (( TASK_ID == 15 )); then
    echo "ERROR: P15 is excluded from current processing; preserve its data/provenance." >&2
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
if [[ "$MAXWELL_QC_ONLY" != "0" && "$MAXWELL_QC_ONLY" != "1" ]]; then
    echo "ERROR: DAD_MAXWELL_QC_ONLY must be 0 or 1." >&2
    exit 2
fi
if [[ "$BLINK_AID_ONLY" != "0" && "$BLINK_AID_ONLY" != "1" ]]; then
    echo "ERROR: DAD_BLINK_AID_ONLY must be 0 or 1." >&2
    exit 2
fi
if [[ "$SQUID_PREPARE_ONLY" != "0" && "$SQUID_PREPARE_ONLY" != "1" ]]; then
    echo "ERROR: DAD_SQUID_PREPARE_ONLY must be 0 or 1." >&2
    exit 2
fi
if [[ "$SQUID_PREPARE_ONLY" == "1" && ( "$BLINK_AID_ONLY" == "1" || "$MAXWELL_QC_ONLY" == "1" || "$STAGE_SPEC" != "maxwell_filter" || "$PROCESSING_PASS" != "review" ) ]]; then
    echo "ERROR: SQUID preparation requires maxwell_filter review only and cannot combine preparation/QC modes." >&2
    exit 2
fi
if [[ "$BLINK_AID_ONLY" == "1" && ( "$MAXWELL_QC_ONLY" == "1" || "$STAGE_SPEC" != "artifact_annotation" || "$PROCESSING_PASS" != "review" ) ]]; then
    echo "ERROR: blink preparation requires artifact_annotation review only and cannot be combined with Maxwell QC." >&2
    exit 2
fi
if [[ "$MAXWELL_QC_ONLY" == "1" && "$STAGE_SPEC" != "maxwell_filter" ]]; then
    echo "ERROR: plot-only mode requires DAD_STAGES=maxwell_filter." >&2
    exit 2
fi
if [[ "$PROCESSING_PASS" == "apply" && "$REBUILD" != "1" && "$MAXWELL_QC_ONLY" != "1" ]]; then
    echo "ERROR: apply requires DAD_REBUILD=1 and explicit overwrite." >&2
    exit 2
fi

# Stage names are colon-separated so the wrapper receives one unambiguous
# environment variable even when several stages are selected.
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
    echo "Set MICROMAMBA_BIN or follow the environment notes in meg/cluster/README.md." >&2
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
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

# Matplotlib and other libraries must not race on a shared home-directory
# cache when many array tasks start together.
JOB_TMP="${TMPDIR:-/tmp}/dad_meg_${SLURM_JOB_ID:-manual}_${TASK_ID}"
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
    "$PROJECT_ROOT/meg/meg_calibration_files/fwmegmaxfilter/ct_sparse_triux2.fif" \
    "$PROJECT_ROOT/meg/meg_calibration_files/fwmegmaxfilter/sss_cal_3182_250410.dat"; do
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

if [[ "$SQUID_PREPARE_ONLY" == "1" ]]; then
    if [[ ! -f "$TUNING_CONFIG" ]]; then
        echo "ERROR: SQUID preparation requires the accepted 4A tuning file: $TUNING_CONFIG" >&2
        exit 2
    fi
    COMMAND=(python "$PROJECT_ROOT/meg/code/preprocessing_basic/review_squid_jumps.py"
        --mode prepare --config "$DEFAULT_CONFIG" --site-config "$SITE_CONFIG"
        --subject-config "$SUBJECT_CONFIG" --tuning-config "$TUNING_CONFIG")
    if [[ -n "${DAD_SOURCE:-}" ]]; then
        COMMAND+=(--source "$DAD_SOURCE")
    fi
    echo ">>> SQUID candidates only: accepted SSS, active 4B reviews and manual annotation files remain read-only."
elif [[ "$BLINK_AID_ONLY" == "1" ]]; then
    if [[ ! -f "$TUNING_CONFIG" ]]; then
        echo "ERROR: blink preparation requires the accepted 4A tuning file: $TUNING_CONFIG" >&2
        exit 2
    fi
    COMMAND=(python "$PROJECT_ROOT/meg/code/preprocessing_basic/blink_annotation.py"
        --config "$DEFAULT_CONFIG" --site-config "$SITE_CONFIG"
        --subject-config "$SUBJECT_CONFIG" --tuning-config "$TUNING_CONFIG")
    if [[ -n "${DAD_SOURCE:-}" ]]; then
        COMMAND+=(--source "$DAD_SOURCE")
    fi
    if [[ "$REBUILD" == "1" ]]; then
        COMMAND+=(--overwrite)
    else
        COMMAND+=(--resume)
    fi
    echo ">>> Supplementary blink preparation only; existing muscle/SSS/manual annotations are read-only."
elif [[ "$MAXWELL_QC_ONLY" == "1" ]]; then
    # Plot settings do not alter scientific YAML or rerun the stage. Existing
    # report source stamps and detected/applied lists are used verbatim.
    COMMAND=(python "$PROJECT_ROOT/meg/code/preprocessing_basic/maxwell_qc.py"
        --derivatives-root "$DERIVATIVES_ROOT" --raw-root "$RAW_ROOT"
        --subjects "$SUBJECT_ID" --processing-pass "$PROCESSING_PASS"
        --window-sec "${DAD_QC_WINDOW_SEC:-20}" --seed "${DAD_QC_SEED:-97}" --overwrite)
    echo ">>> Plot-only refresh: existing scores plus bounded raw windows; no preprocessing."
else
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
fi

echo ">>> Subject: $SUBJECT_LABEL"
echo ">>> Processing pass: $PROCESSING_PASS"
echo ">>> Stages: ${SELECTED_STAGES[*]}"
echo ">>> Host: $(hostname)"
echo ">>> Allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '>>> Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ">>> Processing reports: $REPORT_DIR"
echo ">>> Decision plots: $DERIVATIVES_ROOT/preprocessing/$SUBJECT_LABEL/ses-01/figures"
