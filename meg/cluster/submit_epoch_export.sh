#!/bin/bash
#
# submit_epoch_export.sh
#
# Purpose
# Run production Stage 04 epoching followed by Stage 06 dRSA neural export on
# the CIMeC cluster. One Sun Grid Engine array task processes one participant.
# The two stages run sequentially through run_pipeline.py, so export starts only
# after epoching succeeds.
#
# Submission examples (run from the cluster head node)
#   cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
#   mkdir -p cluster_logs
#   qsub -t 1-7 meg/cluster/submit_epoch_export.sh
#
# Retry selected subjects with guarded resume:
#   qsub -t 2,4,7 meg/cluster/submit_epoch_export.sh
#
# Intentionally rebuild existing Stage-04/06 outputs:
#   qsub -v DAD_REBUILD=1 -t 1-7 meg/cluster/submit_epoch_export.sh
#
# Optional environment/project overrides:
#   qsub -v DAD_MAMBA_ENV=dad-meg -t 1-7 \
#     meg/cluster/submit_epoch_export.sh
#
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

#$ -N DAD_epoch_export
#$ -S /bin/bash
#$ -cwd
#$ -o cluster_logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -j y
#$ -pe smp 4
#$ -l h_vmem=24G
#$ -t 1-7

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve subject, paths, and restart mode
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
TASK_ID="${SGE_TASK_ID:?SGE_TASK_ID is required; submit this as an SGE array job.}"
REBUILD="${DAD_REBUILD:-0}"

if (( TASK_ID < 1 || TASK_ID > 7 )); then
    echo "ERROR: SGE_TASK_ID must be between 1 and 7; got ${TASK_ID}." >&2
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
SUBJECT_CONFIG="$PROJECT_ROOT/meg/params/${SUBJECT_LABEL}.yaml"
REPORT_DIR="$DERIVATIVES_ROOT/reports/${SUBJECT_LABEL}/ses-01"
SOURCE_MANIFEST="$REPORT_DIR/source_manifest.json"


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
export OMP_NUM_THREADS="${NSLOTS:-4}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-4}"
export MKL_NUM_THREADS="${NSLOTS:-4}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-4}"

JOB_TMP="${TMPDIR:-/tmp}/dad_epoch_export_${JOB_ID:-manual}_${TASK_ID}"
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
    --source-manifest "$SOURCE_MANIFEST" \
    --validate-only


# -----------------------------------------------------------------------------
# Run Stage 04 then Stage 06 through the shared orchestration entry point
# -----------------------------------------------------------------------------

COMMAND=(
    python
    "$PROJECT_ROOT/meg/code/preprocessing_basic/run_pipeline.py"
    --config "$DEFAULT_CONFIG"
    --site-config "$SITE_CONFIG"
    --subject-config "$SUBJECT_CONFIG"
    --stages epoching export_drsa
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
echo ">>> Slots: ${NSLOTS:-unknown}"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '>>> Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ">>> Reports: $REPORT_DIR"
echo ">>> Epochs: $DERIVATIVES_ROOT/analysis/$SUBJECT_LABEL/ses-01/meg"
echo ">>> dRSA exports: $DERIVATIVES_ROOT/analysis/$SUBJECT_LABEL/ses-01/drsa"
