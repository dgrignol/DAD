#!/bin/bash
#
# submit_drsa_analysis.sh
#
# Purpose
# Run the optional Stage-07 MATLAB dRSA analysis layer on the DIPSEN cluster.
# One SLURM array task processes one participant. This stage consumes
# completed Stage-06 identity-averaged neural exports and writes MATLAB dRSA
# PNG/MAT outputs plus JSON/HTML reports.
# The guarded Stage-07 CLI runs directly, preserving earlier Stage-1 inventory
# and resolved-configuration snapshots instead of re-entering orchestration.
#
# Submission examples (run from the DAD root on the DIPSEN login node)
#   mkdir -p cluster_logs
#   sbatch --array=4-14,16-39 meg/cluster/submit_drsa_analysis.sh
# DIPSEN provides MATLAB through modules; the wrapper prefers R2020a and can
# be given another executable if needed:
#   DAD_MATLAB_BIN=/absolute/path/to/matlab sbatch --array=4-39 \
#     meg/cluster/submit_drsa_analysis.sh
# Optional analysis selection can be exported before submission:
#   DAD_DRSA_ANALYSES=real \
#   DAD_DRSA_SENSOR_GROUPS=wholeBrain \
#   DAD_DRSA_MODEL_SOURCES=observed,predicted,prediction_error,pe_prior_constrained_50 \
#     sbatch --array=4-39 meg/cluster/submit_drsa_analysis.sh
# Keep DAD_EPOCH_PREPARATION and DAD_PROCESSING_RUN exported from epoching.
# They select matching exports, reports and analysis outputs.
# Reviewed participant tuning is loaded from the original Stage-03 reports,
# including P38's exact-source filtering exception. Override its path with
# DAD_DRSA_TUNING_CONFIG if the reviewed report tree is stored elsewhere.
# DAD_SUBJECT_CONFIG_DIR must match any explicit epoch-recovery selector folder.
# Inputs
# - Completed Stage-06 export_drsa lineage in source_manifest.json.
# - Stage-06 condition MAT/JSON neural exports.
# - Subject-specific observed trajectory MAT. Matching occluded_nondeviant
#   trials provide the expected continuation for occluded_deviant and the
#   baseline for prediction_error / pe_* alt model sources.
# - meg/code/analysis MATLAB scripts and simulations/functions dRSA helpers.
#
# Outputs
# - MATLAB dRSA PNG/MAT outputs under derivatives/analysis/sub-XX/ses-01/drsa.
# - 07_matlab_drsa_handoff.json, JSON report, and HTML report under reports.

#SBATCH --job-name=DAD_drsa
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL
#SBATCH --array=4-14,16-39
# MATLAB R2020a Update 4 reserves up to 1 GB of JVM metaspace during startup.
# Keep the validated 24 GB task limit; a 2 GB probe aborted before MATLAB ran.

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve subject, paths, and executable settings
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required; submit this as an SLURM array job.}"
REBUILD="${DAD_REBUILD:-0}"

resolve_matlab_bin() {
    local requested="${DAD_MATLAB_BIN:-}"
    local candidate=""
    local init_script=""

    if [[ -n "$requested" ]]; then
        if [[ "$requested" == /* && -x "$requested" ]]; then
            printf '%s\n' "$requested"
            return 0
        fi
        echo "ERROR: DAD_MATLAB_BIN must be an executable absolute path: $requested" >&2
        return 2
    fi

    # Load the DIPSEN compatibility module before consulting PATH so an
    # unrelated MATLAB installation inherited from the login shell is not
    # selected accidentally.
    for init_script in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
        if [[ -f "$init_script" ]]; then
            # shellcheck source=/dev/null
            source "$init_script"
            break
        fi
    done
    if command -v module >/dev/null 2>&1; then
        module load matlab/R2020a >/dev/null 2>&1 \
            || module load matlab/R2025b >/dev/null 2>&1 \
            || module load matlab >/dev/null 2>&1 \
            || true
    fi

    if command -v matlab >/dev/null 2>&1; then
        candidate="$(command -v matlab)"
        if [[ "$candidate" == /* && -x "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    fi

    for candidate in \
        /usr/local/MATLAB/R2020a/bin/matlab \
        /usr/local/MATLAB/R2021a/bin/matlab \
        /usr/local/MATLAB/R2022a/bin/matlab \
        /usr/local/MATLAB/R2022b/bin/matlab \
        /usr/local/MATLAB/R2023a/bin/matlab \
        /opt/MATLAB/R2020a/bin/matlab \
        /opt/matlab/R2020a/bin/matlab \
        /cm/shared/apps/matlab/R2020a/bin/matlab \
        /apps/matlab/R2020a/bin/matlab; do
        if [[ -x "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done

    echo "ERROR: unable to locate MATLAB." >&2
    echo "       Set DAD_MATLAB_BIN=/absolute/path/to/matlab before sbatch." >&2
    echo "       Useful checks on the login node: command -v matlab; module avail matlab." >&2
    return 2
}

if (( TASK_ID < 4 || TASK_ID > 39 )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID must be between 4 and 39; got ${TASK_ID}." >&2
    exit 2
fi
if (( TASK_ID == 15 )); then
    echo "ERROR: P15 is excluded from current analysis; preserve its data/provenance." >&2
    exit 2
fi
if [[ "$REBUILD" != "0" && "$REBUILD" != "1" ]]; then
    echo "ERROR: DAD_REBUILD must be 0 or 1; got ${REBUILD}." >&2
    exit 2
fi
MATLAB_BIN="$(resolve_matlab_bin)"

SUBJECT_ID="$(printf "%02d" "$TASK_ID")"
SUBJECT_LABEL="sub-${SUBJECT_ID}"
DEFAULT_CONFIG="$PROJECT_ROOT/meg/params/pipeline_default.yaml"
SITE_CONFIG="$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
SUBJECT_CONFIG_DIR="${DAD_SUBJECT_CONFIG_DIR:-$PROJECT_ROOT/meg/params}"
SUBJECT_CONFIG="$SUBJECT_CONFIG_DIR/${SUBJECT_LABEL}.yaml"
DOWNSTREAM_ROOT="$DERIVATIVES_ROOT"
if [[ -n "${DAD_PROCESSING_RUN:-}" ]]; then
    DOWNSTREAM_ROOT="$DERIVATIVES_ROOT/processing_runs/$DAD_PROCESSING_RUN"
fi
REPORT_DIR="$DOWNSTREAM_ROOT/reports/${SUBJECT_LABEL}/ses-01"
SOURCE_MANIFEST="$REPORT_DIR/source_manifest.json"
TUNING_CONFIG="${DAD_DRSA_TUNING_CONFIG:-$DERIVATIVES_ROOT/reports/$SUBJECT_LABEL/ses-01/03_tuning_choices.yaml}"


# -----------------------------------------------------------------------------
# Activate Python environment and configure headless MATLAB/plotting
# -----------------------------------------------------------------------------

if [[ ! -x "$MAMBA_BIN" ]]; then
    echo "ERROR: micromamba executable not found at $MAMBA_BIN." >&2
    exit 2
fi

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "$MAMBA_ENV"

export DAD_MATLAB_BIN="$MATLAB_BIN"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export DAD_MATLAB_MAX_THREADS="${SLURM_CPUS_PER_TASK:-4}"

JOB_TMP="${TMPDIR:-/tmp}/dad_drsa_${SLURM_JOB_ID:-manual}_${TASK_ID}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib"
export XDG_CACHE_HOME="$JOB_TMP/cache"


# -----------------------------------------------------------------------------
# Validate deployment before launching MATLAB
# -----------------------------------------------------------------------------

cd "$PROJECT_ROOT"
mkdir -p cluster_logs

for required_path in \
    "$DEFAULT_CONFIG" \
    "$SITE_CONFIG" \
    "$SUBJECT_CONFIG" \
    "$SOURCE_MANIFEST" \
    "$TUNING_CONFIG" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/07_run_matlab_drsa.py" \
    "$PROJECT_ROOT/meg/code/analysis/A_dRSA.m" \
    "$PROJECT_ROOT/meg/code/analysis/A_dRSA_alt.m" \
    "$PROJECT_ROOT/meg/code/analysis/dRSA_simulations.m" \
    "$PROJECT_ROOT/meg/code/analysis/drsa_stage_utils.m" \
    "$PROJECT_ROOT/simulations/functions/dRSA_coreFunction.m"; do
    if [[ ! -f "$required_path" ]]; then
        echo "ERROR: required file is missing: $required_path" >&2
        exit 2
    fi
done

python "$PROJECT_ROOT/meg/code/preprocessing_basic/07_run_matlab_drsa.py" \
    --config "$DEFAULT_CONFIG" \
    --site-config "$SITE_CONFIG" \
    --subject-config "$SUBJECT_CONFIG" \
    --tuning-config "$TUNING_CONFIG" \
    --source-manifest "$SOURCE_MANIFEST" \
    --validate-only


# -----------------------------------------------------------------------------
# Run Stage 07
# -----------------------------------------------------------------------------

COMMAND=(
    python
    "$PROJECT_ROOT/meg/code/preprocessing_basic/07_run_matlab_drsa.py"
    --config "$DEFAULT_CONFIG"
    --site-config "$SITE_CONFIG"
    --subject-config "$SUBJECT_CONFIG"
    --tuning-config "$TUNING_CONFIG"
    --source-manifest "$SOURCE_MANIFEST"
    --print-summary
)
if [[ "$REBUILD" == "1" ]]; then
    COMMAND+=(--overwrite)
    echo ">>> Rebuild mode: MATLAB dRSA outputs may be replaced."
else
    COMMAND+=(--resume)
    echo ">>> Resume mode: completed MATLAB dRSA outputs are reused."
fi

echo ">>> Subject: $SUBJECT_LABEL"
echo ">>> MATLAB: $MATLAB_BIN"
echo ">>> Analyses: ${DAD_DRSA_ANALYSES:-config-default}"
echo ">>> Conditions: ${DAD_DRSA_CONDITIONS:-config-default}"
echo ">>> Sensor groups: ${DAD_DRSA_SENSOR_GROUPS:-config-default}"
echo ">>> Model sources: ${DAD_DRSA_MODEL_SOURCES:-config-default}"
echo ">>> Smoothing plots: ${DAD_DRSA_SMOOTHING_PLOTS:-config-default}"
echo ">>> Host: $(hostname)"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '>>> Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ">>> Reports: $REPORT_DIR"
echo ">>> dRSA outputs: $DERIVATIVES_ROOT/analysis/$SUBJECT_LABEL/ses-01/drsa"
