#!/bin/bash
#
# submit_drsa_analysis.sh
#
# Purpose
# Run the optional Stage-07 MATLAB dRSA analysis layer on the CIMeC cluster.
# One Sun Grid Engine array task processes one participant. This stage consumes
# completed Stage-06 identity-averaged neural exports and writes MATLAB dRSA
# PNG/MAT outputs plus JSON/HTML reports.
#
# Submission examples (run from the cluster head node)
#   cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
#   mkdir -p cluster_logs
#   qsub -t 1-7 meg/cluster/submit_drsa_analysis.sh
#
# The wrapper first checks the confirmed CIMeC R2020a installation, then
# common MATLAB locations, PATH, and modules. An explicit override remains
# available for cluster changes or alternate installations:
#   qsub -v DAD_MATLAB_BIN=/state/partition1/MATLAB/R2020a/bin/matlab \
#     -t 1-7 meg/cluster/submit_drsa_analysis.sh
#
# Optional analysis selection can be exported before qsub when needed:
#   export DAD_DRSA_ANALYSES=real
#   export DAD_DRSA_SENSOR_GROUPS=wholeBrain
#   export DAD_DRSA_MODEL_SOURCES=observed,predicted
#   qsub -v DAD_MATLAB_BIN,DAD_DRSA_ANALYSES,DAD_DRSA_SENSOR_GROUPS,DAD_DRSA_MODEL_SOURCES \
#     -t 1-7 meg/cluster/submit_drsa_analysis.sh
#
# Inputs
# - Completed Stage-06 export_drsa lineage in source_manifest.json.
# - Stage-06 condition MAT/JSON neural exports.
# - Subject-specific observed trajectory MAT. Matching occluded_nondeviant
#   trials provide the expected continuation for occluded_deviant.
# - meg/code/analysis MATLAB scripts and simulations/functions dRSA helpers.
#
# Outputs
# - MATLAB dRSA PNG/MAT outputs under derivatives/analysis/sub-XX/ses-01/drsa.
# - 07_matlab_drsa_handoff.json, JSON report, and HTML report under reports.

#$ -N DAD_drsa
#$ -S /bin/bash
#$ -cwd
#$ -o cluster_logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -j y
#$ -pe smp 4
# MATLAB R2020a Update 4 reserves up to 1 GB of JVM metaspace during startup.
# Keep the validated 24 GB task limit; a 2 GB probe aborted before MATLAB ran.
#$ -l h_vmem=24G
#$ -t 1-7

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve subject, paths, and executable settings
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
TASK_ID="${SGE_TASK_ID:?SGE_TASK_ID is required; submit this as an SGE array job.}"
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

    # This is the MATLAB installation exposed on CIMeC compute nodes. Keep it
    # first so normal submissions do not depend on login-shell PATH or modules.
    for candidate in \
        /state/partition1/MATLAB/R2020a/bin/matlab \
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

    if command -v matlab >/dev/null 2>&1; then
        candidate="$(command -v matlab)"
        if [[ "$candidate" == /* && -x "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    fi

    for init_script in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
        if [[ -f "$init_script" ]]; then
            # shellcheck source=/dev/null
            source "$init_script"
            break
        fi
    done
    if command -v module >/dev/null 2>&1; then
        module load matlab/R2020a >/dev/null 2>&1 || module load matlab >/dev/null 2>&1 || true
        if command -v matlab >/dev/null 2>&1; then
            candidate="$(command -v matlab)"
            if [[ "$candidate" == /* && -x "$candidate" ]]; then
                printf '%s\n' "$candidate"
                return 0
            fi
        fi
    fi

    echo "ERROR: unable to locate MATLAB." >&2
    echo "       Resubmit with -v DAD_MATLAB_BIN=/absolute/path/to/matlab." >&2
    echo "       Useful checks on the head node: command -v matlab; module avail matlab." >&2
    return 2
}

if (( TASK_ID < 1 || TASK_ID > 7 )); then
    echo "ERROR: SGE_TASK_ID must be between 1 and 7; got ${TASK_ID}." >&2
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
SUBJECT_CONFIG="$PROJECT_ROOT/meg/params/${SUBJECT_LABEL}.yaml"
REPORT_DIR="$DERIVATIVES_ROOT/reports/${SUBJECT_LABEL}/ses-01"
SOURCE_MANIFEST="$REPORT_DIR/source_manifest.json"


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
export OMP_NUM_THREADS="${NSLOTS:-4}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-4}"
export MKL_NUM_THREADS="${NSLOTS:-4}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-4}"

JOB_TMP="${TMPDIR:-/tmp}/dad_drsa_${JOB_ID:-manual}_${TASK_ID}"
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
    --source-manifest "$SOURCE_MANIFEST" \
    --validate-only


# -----------------------------------------------------------------------------
# Run Stage 07
# -----------------------------------------------------------------------------

COMMAND=(
    python
    "$PROJECT_ROOT/meg/code/preprocessing_basic/run_pipeline.py"
    --config "$DEFAULT_CONFIG"
    --site-config "$SITE_CONFIG"
    --subject-config "$SUBJECT_CONFIG"
    --source-manifest "$SOURCE_MANIFEST"
    --stages matlab_drsa
    --skip-preflight
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
