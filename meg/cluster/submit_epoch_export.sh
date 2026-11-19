#!/bin/bash
#
# submit_epoch_export.sh
#
# Purpose
# Run production Stage 04 epoching followed by Stage 06 dRSA neural export on
# the DIPSEN cluster. One SLURM array task processes one participant.
# The two guarded stage CLIs run sequentially, so export starts only after
# epoching succeeds. This preserves earlier Stage 1 inventory/config snapshots.
# Accepted subject tuning and the hash-bound s3-v002 review are passed explicitly;
# DAD_TIMING_REVIEW_DIR can select another already accepted review directory.
#
# Submission examples (run from the DAD root on the DIPSEN login node)
#   mkdir -p cluster_logs
#   sbatch --array=4-14,16-39 meg/cluster/submit_epoch_export.sh
# Retry selected subjects with guarded resume:
#   sbatch --array=4,7,24 meg/cluster/submit_epoch_export.sh
# Intentionally rebuild existing Stage-04/06 outputs:
#   DAD_REBUILD=1 sbatch --array=4-14,16-39 \
#     meg/cluster/submit_epoch_export.sh
# Environment values can be set inline before sbatch, for example:
#   DAD_MAMBA_ENV=dad-meg sbatch --array=4-14,16-39 \
#     meg/cluster/submit_epoch_export.sh
# Storage-saving modes (accepted annotated inputs; no continuous output FIF):
#   DAD_EPOCH_PREPARATION=inline_no_ica DAD_PROCESSING_RUN=noica-20261112-v1 \
#     sbatch --cpus-per-task=1 --array=4-14,16-39 meg/cluster/submit_epoch_export.sh
# Future reviewed ICA: inline_ica with a NEW processing-run label; ICA runs
# before the same continuous spectrum-fit. Model-bound Stage 4C is required.
# DAD_VALIDATE_ONLY=1 stops after read-only input/timing validation.
# DAD_SUBJECT_CONFIG_DIR selects an explicit directory of sub-XX.yaml overrides
# for recovery arrays; every participant still has its own independent config.
# Labelled epochs are immutable: resume checks exact contracts; use a new
# label after changing source/model/choices/settings or an incomplete write.
# Inputs
# - Completed manual-ICA apply lineage in source_manifest.json.
# - One lineage-matched <exact-source-stem>_ica.fif per acquisition run.
# - Canonical 02_trigger_events.csv.
# - Subject-specific observed trajectory MAT containing all three conditions.
# - Shared/site/subject configuration and fixed canonical MEG channel axis.
#
# Outputs
# - Stage-04 epochs FIF, epoch metadata, complete trial decisions, JSON/HTML.
# - Stage-06 condition MAT/JSON exports and aggregate JSON/HTML.
# - Updated epoching/export_drsa derivative lineage in source_manifest.json.
#
# Safety and restart behavior
# - Raw data are read-only and all outputs pass the shared raw-tree guard.
# - Default mode uses --resume. If Stage 04 is complete and Stage 06 failed,
#   Stage 04 is skipped and export is retried.
# - DAD_REBUILD=1 uses --overwrite and intentionally replaces Stage-04/06
#   derivatives. It does not rerun Stage 03.

#SBATCH --job-name=DAD_epoch_export
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL
#SBATCH --array=4-14,16-39

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve subject, paths, and restart mode
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required; submit this as an SLURM array job.}"
REBUILD="${DAD_REBUILD:-0}"
VALIDATE_ONLY="${DAD_VALIDATE_ONLY:-0}"
if [[ "$VALIDATE_ONLY" != "0" && "$VALIDATE_ONLY" != "1" ]]; then
    echo "ERROR: DAD_VALIDATE_ONLY must be 0 or 1." >&2
    exit 2
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

SUBJECT_ID="$(printf "%02d" "$TASK_ID")"
SUBJECT_LABEL="sub-${SUBJECT_ID}"
DEFAULT_CONFIG="$PROJECT_ROOT/meg/params/pipeline_default.yaml"
SITE_CONFIG="$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
SUBJECT_CONFIG_DIR="${DAD_SUBJECT_CONFIG_DIR:-$PROJECT_ROOT/meg/params}"
SUBJECT_CONFIG="$SUBJECT_CONFIG_DIR/${SUBJECT_LABEL}.yaml"
REPORT_DIR="$DERIVATIVES_ROOT/reports/${SUBJECT_LABEL}/ses-01"
SOURCE_MANIFEST="$REPORT_DIR/source_manifest.json"
TUNING_CONFIG="$REPORT_DIR/03_tuning_choices.yaml"
TIMING_REVIEW_DIR="${DAD_TIMING_REVIEW_DIR:-$DERIVATIVES_ROOT/visual_timing/s3-v002/reviews}"
TIMING_REVIEW="$TIMING_REVIEW_DIR/${SUBJECT_LABEL}_review.json"
DOWNSTREAM_ROOT="$DERIVATIVES_ROOT"
if [[ -n "${DAD_PROCESSING_RUN:-}" ]]; then
    DOWNSTREAM_ROOT="$DERIVATIVES_ROOT/processing_runs/$DAD_PROCESSING_RUN"
fi
DOWNSTREAM_MANIFEST="$DOWNSTREAM_ROOT/reports/${SUBJECT_LABEL}/ses-01/source_manifest.json"


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
export MNE_BROWSER_BACKEND=matplotlib
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

JOB_TMP="${TMPDIR:-/tmp}/dad_epoch_export_${SLURM_JOB_ID:-manual}_${TASK_ID}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib"
export XDG_CACHE_HOME="$JOB_TMP/cache"


# -----------------------------------------------------------------------------
# Validate deployment and Stage-03 prerequisites before processing
# -----------------------------------------------------------------------------

cd "$PROJECT_ROOT"
mkdir -p cluster_logs

for required_path in \
    "$DEFAULT_CONFIG" \
    "$SITE_CONFIG" \
    "$SUBJECT_CONFIG" \
    "$SOURCE_MANIFEST" \
    "$TUNING_CONFIG" \
    "$TIMING_REVIEW" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/04_epoching.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/06_export_neural_for_drsa.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/downstream_stage_utils.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/drsa_canonical_meg_channels.txt"; do
    if [[ ! -f "$required_path" ]]; then
        echo "ERROR: required file is missing: $required_path" >&2
        exit 2
    fi
done

python - <<'PY'
import mne
import numpy
import pandas
import scipy
import yaml

print(
    "Runtime imports OK:",
    f"mne={mne.__version__}",
    f"numpy={numpy.__version__}",
    f"scipy={scipy.__version__}",
    f"pandas={pandas.__version__}",
)
PY

python "$PROJECT_ROOT/meg/code/preprocessing_basic/04_epoching.py" \
    --config "$DEFAULT_CONFIG" \
    --site-config "$SITE_CONFIG" \
    --subject-config "$SUBJECT_CONFIG" \
    --tuning-config "$TUNING_CONFIG" \
    --source-manifest "$SOURCE_MANIFEST" \
    --timing-review "$TIMING_REVIEW" \
    --validate-only
if [[ "$VALIDATE_ONLY" == "1" ]]; then
    echo ">>> ${SUBJECT_LABEL}: input/timing validation passed; no derivatives produced."
    exit 0
fi


# -----------------------------------------------------------------------------
# Run Stage 04 then Stage 06 with accepted tuning/timing and guarded resume
# -----------------------------------------------------------------------------

COMMAND=(
    python
    "$PROJECT_ROOT/meg/code/preprocessing_basic/04_epoching.py"
    --config "$DEFAULT_CONFIG"
    --site-config "$SITE_CONFIG"
    --subject-config "$SUBJECT_CONFIG"
    --tuning-config "$TUNING_CONFIG"
    --source-manifest "$SOURCE_MANIFEST"
    --timing-review "$TIMING_REVIEW"
    --print-summary
)
if [[ "$REBUILD" == "1" ]]; then
    COMMAND+=(--overwrite)
    echo ">>> Rebuild mode: Stage-04/06 outputs will be replaced."
else
    COMMAND+=(--resume)
    echo ">>> Resume mode: completed Stage-04/06 reports are reused."
fi

echo ">>> Subject: $SUBJECT_LABEL"
echo ">>> Host: $(hostname)"
echo ">>> Allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '>>> Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

# Export consumes only the completed, lineage-matched epochs from the first CLI.
EXPORT_COMMAND=(python "$PROJECT_ROOT/meg/code/preprocessing_basic/06_export_neural_for_drsa.py"
    --config "$DEFAULT_CONFIG" --site-config "$SITE_CONFIG"
    --subject-config "$SUBJECT_CONFIG" --tuning-config "$TUNING_CONFIG"
    --source-manifest "$DOWNSTREAM_MANIFEST" --print-summary)
if [[ "$REBUILD" == "1" ]]; then
    EXPORT_COMMAND+=(--overwrite)
else
    EXPORT_COMMAND+=(--resume)
fi
"${EXPORT_COMMAND[@]}"

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ">>> Reports: $DOWNSTREAM_ROOT/reports/$SUBJECT_LABEL/ses-01"
echo ">>> Epochs: $DOWNSTREAM_ROOT/analysis/$SUBJECT_LABEL/ses-01/meg"
echo ">>> dRSA exports: $DOWNSTREAM_ROOT/analysis/$SUBJECT_LABEL/ses-01/drsa"
