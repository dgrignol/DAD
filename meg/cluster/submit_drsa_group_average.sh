#!/bin/bash
#
# submit_drsa_group_average.sh
#
# Purpose
# Build the optional Stage-07 numeric group-average dRSA report on the DIPSEN
# cluster from existing subject-level Stage-07 MATLAB result MAT files. This is
# a one-task aggregation job; it does not run MATLAB and does not rerun Stage 04,
# Stage 06, or the subject-level Stage 07 dRSA jobs.
#
# Submission examples (run from the DAD root on the DIPSEN login node)
#   mkdir -p cluster_logs
#   sbatch meg/cluster/submit_drsa_group_average.sh
# Submit after the Stage-07 subject array finishes:
#   sbatch --dependency=afterok:<stage07_job_id> \
#     meg/cluster/submit_drsa_group_average.sh
# Override included subjects and disable smoothed group figures:
#   DAD_DRSA_GROUP_SUBJECTS=01,02,03 \
#   DAD_DRSA_GROUP_SMOOTHING_PLOTS=0 \
#     sbatch meg/cluster/submit_drsa_group_average.sh
# Inputs
# - Existing Stage-07 result MAT files below
#   derivatives/analysis/sub-XX/ses-01/drsa/.
#
# Outputs
# - Group-average PNG outputs below
#   derivatives/analysis/group/ses-01/drsa/group_average/.
# - Group raw/smoothed JSON and HTML reports below
#   derivatives/reports/group/ses-01/.
#
# Environment variables
# - DAD_PROJECT_ROOT: cluster DAD root.
# - DAD_MAMBA_ENV, MICROMAMBA_BIN: Python environment selection.
# - DAD_DRSA_GROUP_CONFIG: pipeline default YAML path.
# - DAD_DRSA_GROUP_SITE_CONFIG: cluster/site YAML path.
# - DAD_DRSA_GROUP_SUBJECT_CONFIG: representative subject YAML path.
# - DAD_DRSA_GROUP_SUBJECTS: comma-separated subjects, default 04..39.
# - DAD_DRSA_GROUP_SMOOTHING_PLOTS: 1 writes smoothed group report, 0 disables it.
# - DAD_DRSA_GROUP_PLOT_SMOOTHING_WINDOW: optional smoothing window in samples.

#SBATCH --job-name=DAD_drsa_groupavg
#SBATCH --output=cluster_logs/%x.%j.log
#SBATCH --error=cluster_logs/%x.%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH --export=ALL

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve paths and group-average settings
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
DEFAULT_CONFIG="${DAD_DRSA_GROUP_CONFIG:-$PROJECT_ROOT/meg/params/pipeline_default.yaml}"
SITE_CONFIG="${DAD_DRSA_GROUP_SITE_CONFIG:-$PROJECT_ROOT/meg/params/pipeline_cluster.yaml}"
SUBJECT_CONFIG="${DAD_DRSA_GROUP_SUBJECT_CONFIG:-$PROJECT_ROOT/meg/params/sub-04.yaml}"
SUBJECTS_CSV="${DAD_DRSA_GROUP_SUBJECTS:-04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39}"
SMOOTHING_PLOTS="${DAD_DRSA_GROUP_SMOOTHING_PLOTS:-1}"
PLOT_SMOOTHING_WINDOW="${DAD_DRSA_GROUP_PLOT_SMOOTHING_WINDOW:-}"
SCRIPT="$PROJECT_ROOT/meg/code/preprocessing_basic/07_run_matlab_drsa.py"
IFS=',' read -r -a GROUP_SUBJECTS <<< "$SUBJECTS_CSV"

if [[ "$SMOOTHING_PLOTS" != "0" && "$SMOOTHING_PLOTS" != "1" ]]; then
    echo "ERROR: DAD_DRSA_GROUP_SMOOTHING_PLOTS must be 0 or 1; got ${SMOOTHING_PLOTS}." >&2
    exit 2
fi
if [[ "${#GROUP_SUBJECTS[@]}" -eq 0 ]]; then
    echo "ERROR: DAD_DRSA_GROUP_SUBJECTS resolved to an empty subject list." >&2
    exit 2
fi


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
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

JOB_TMP="${TMPDIR:-/tmp}/dad_drsa_groupavg_${SLURM_JOB_ID:-manual}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib"
export XDG_CACHE_HOME="$JOB_TMP/cache"


# -----------------------------------------------------------------------------
# Validate deployment and build the group-average report
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

python - <<'PY'
import h5py
import scipy
print(f">>> Python MAT readers: h5py={h5py.__version__}, scipy={scipy.__version__}")
PY

COMMAND=(
    python
    "$SCRIPT"
    --config "$DEFAULT_CONFIG"
    --site-config "$SITE_CONFIG"
    --subject-config "$SUBJECT_CONFIG"
    --group-average-only
    --group-subjects "${GROUP_SUBJECTS[@]}"
    --print-summary
)

if [[ "$SMOOTHING_PLOTS" == "1" ]]; then
    COMMAND+=(--smoothing-plots)
else
    COMMAND+=(--no-smoothing-plots)
fi

if [[ -n "$PLOT_SMOOTHING_WINDOW" ]]; then
    COMMAND+=(--plot-smoothing-window "$PLOT_SMOOTHING_WINDOW")
fi

echo ">>> Group subjects: ${GROUP_SUBJECTS[*]}"
echo ">>> Config: $DEFAULT_CONFIG"
echo ">>> Site config: $SITE_CONFIG"
echo ">>> Subject config: $SUBJECT_CONFIG"
echo ">>> Smoothing plots: $SMOOTHING_PLOTS"
echo ">>> Host: $(hostname)"
echo ">>> Allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '>>> Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
