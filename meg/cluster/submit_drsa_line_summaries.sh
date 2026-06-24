#!/bin/bash
#
# submit_drsa_line_summaries.sh
#
# Purpose
# Build the optional Stage-08 group dRSA line-summary report on the CIMeC
# cluster from existing Stage-07 result MAT files. This is a one-task group
# aggregation job; it does not run MATLAB and does not rerun Stage 04, Stage 06,
# or Stage 07 dRSA.
#
# Submission examples (run from the cluster head node)
#   cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
#   mkdir -p cluster_logs
#   qsub meg/cluster/submit_drsa_line_summaries.sh
#
# Rebuild existing Stage-08 group line-summary outputs:
#   qsub -v DAD_REBUILD=1 meg/cluster/submit_drsa_line_summaries.sh
#
# Rebuild with descriptive within-subject diagonal-SD shading:
#   qsub -v DAD_REBUILD=1,DAD_DRSA_LINE_SHADING=within-subject-diagonal-sd \
#     meg/cluster/submit_drsa_line_summaries.sh
#
# Rebuild with strict same-interval square support for interval lag curves:
#   qsub -v DAD_REBUILD=1,DAD_DRSA_LINE_INTERVAL_LAG_MODE=same-interval-square \
#     meg/cluster/submit_drsa_line_summaries.sh
#
# Restrict to occipital Stage-07 result files and use a +/-0.5 s lag window:
#   qsub -v DAD_DRSA_LINE_SENSOR_GROUPS=occipital,DAD_DRSA_LINE_LAG_WINDOW_SEC=0.5 \
#     meg/cluster/submit_drsa_line_summaries.sh
#
# Inputs
# - Existing Stage-07 result MAT files below
#   derivatives/analysis/sub-XX/ses-01/drsa/.
#
# Outputs
# - Group PNG outputs below
#   derivatives/analysis/group/ses-01/drsa/line_summaries/.
# - Group linked and embedded HTML reports below
#   derivatives/reports/group/ses-01/.
#
# Environment variables
# - DAD_PROJECT_ROOT: cluster DAD root.
# - DAD_MAMBA_ENV, MICROMAMBA_BIN: Python environment selection.
# - DAD_DRSA_LINE_SUBJECTS: comma-separated subjects, default 01..07.
# - DAD_DRSA_LINE_CONDITIONS: optional comma-separated condition filter.
# - DAD_DRSA_LINE_SENSOR_GROUPS: optional comma-separated sensor-group filter.
# - DAD_DRSA_LINE_LAG_WINDOW_SEC: lag search half-window, default 1.0.
# - DAD_DRSA_LINE_INCLUDE_SMOOTHED: 1 also writes display-smoothed summaries.
# - DAD_DRSA_LINE_PLOT_SMOOTHING_WINDOW: optional smoothing window in samples.
# - DAD_DRSA_LINE_SHADING: interval lag-curve band mode, default none.
# - DAD_DRSA_LINE_INTERVAL_LAG_MODE: model-anchored or same-interval-square.
# - DAD_REBUILD: 1 overwrites existing outputs; 0 resumes existing outputs.

#$ -N DAD_drsa_lines
#$ -S /bin/bash
#$ -cwd
#$ -o cluster_logs/$JOB_NAME.$JOB_ID.log
#$ -j y
#$ -pe smp 2
#$ -l h_vmem=12G

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve paths and Stage-08 settings
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
REBUILD="${DAD_REBUILD:-0}"
SUBJECTS_CSV="${DAD_DRSA_LINE_SUBJECTS:-01,02,03,04,05,06,07}"
CONDITIONS_CSV="${DAD_DRSA_LINE_CONDITIONS:-}"
SENSOR_GROUPS_CSV="${DAD_DRSA_LINE_SENSOR_GROUPS:-}"
LAG_WINDOW_SEC="${DAD_DRSA_LINE_LAG_WINDOW_SEC:-1.0}"
INCLUDE_SMOOTHED="${DAD_DRSA_LINE_INCLUDE_SMOOTHED:-1}"
PLOT_SMOOTHING_WINDOW="${DAD_DRSA_LINE_PLOT_SMOOTHING_WINDOW:-}"
SHADING="${DAD_DRSA_LINE_SHADING:-none}"
INTERVAL_LAG_MODE="${DAD_DRSA_LINE_INTERVAL_LAG_MODE:-model-anchored}"

if [[ "$REBUILD" != "0" && "$REBUILD" != "1" ]]; then
    echo "ERROR: DAD_REBUILD must be 0 or 1; got ${REBUILD}." >&2
    exit 2
fi
if [[ "$INCLUDE_SMOOTHED" != "0" && "$INCLUDE_SMOOTHED" != "1" ]]; then
    echo "ERROR: DAD_DRSA_LINE_INCLUDE_SMOOTHED must be 0 or 1; got ${INCLUDE_SMOOTHED}." >&2
    exit 2
fi
case "$SHADING" in
    none|within-subject-diagonal-sd|across-subject-sem)
        ;;
    *)
        echo "ERROR: DAD_DRSA_LINE_SHADING must be none, within-subject-diagonal-sd, or across-subject-sem; got ${SHADING}." >&2
        exit 2
        ;;
esac
case "$INTERVAL_LAG_MODE" in
    model-anchored|same-interval-square)
        ;;
    *)
        echo "ERROR: DAD_DRSA_LINE_INTERVAL_LAG_MODE must be model-anchored or same-interval-square; got ${INTERVAL_LAG_MODE}." >&2
        exit 2
        ;;
esac

DEFAULT_CONFIG="$PROJECT_ROOT/meg/params/pipeline_default.yaml"
SITE_CONFIG="$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
SCRIPT="$PROJECT_ROOT/meg/code/preprocessing_basic/08_plot_drsa_line_summaries.py"
IFS=',' read -r -a SUBJECTS <<< "$SUBJECTS_CSV"


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

JOB_TMP="${TMPDIR:-/tmp}/dad_drsa_lines_${JOB_ID:-manual}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib"
export XDG_CACHE_HOME="$JOB_TMP/cache"


# -----------------------------------------------------------------------------
# Validate deployment and build the group line-summary report
# -----------------------------------------------------------------------------

cd "$PROJECT_ROOT"
mkdir -p cluster_logs

for required_path in \
    "$DEFAULT_CONFIG" \
    "$SITE_CONFIG" \
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
    --subjects "${SUBJECTS[@]}"
    --group-only
    --lag-window-sec "$LAG_WINDOW_SEC"
    --shading "$SHADING"
    --interval-lag-mode "$INTERVAL_LAG_MODE"
    --print-summary
)

if [[ "$INCLUDE_SMOOTHED" == "1" ]]; then
    COMMAND+=(--include-smoothed)
fi

if [[ -n "$PLOT_SMOOTHING_WINDOW" ]]; then
    COMMAND+=(--plot-smoothing-window "$PLOT_SMOOTHING_WINDOW")
fi

if [[ -n "$CONDITIONS_CSV" ]]; then
    IFS=',' read -r -a CONDITIONS <<< "$CONDITIONS_CSV"
    COMMAND+=(--conditions "${CONDITIONS[@]}")
fi

if [[ -n "$SENSOR_GROUPS_CSV" ]]; then
    IFS=',' read -r -a SENSOR_GROUPS <<< "$SENSOR_GROUPS_CSV"
    COMMAND+=(--sensor-groups "${SENSOR_GROUPS[@]}")
fi

if [[ "$REBUILD" == "1" ]]; then
    COMMAND+=(--overwrite)
    echo ">>> Rebuild mode: Stage-08 line-summary outputs may be replaced."
else
    COMMAND+=(--resume)
    echo ">>> Resume mode: existing Stage-08 line-summary outputs are reused."
fi

echo ">>> Subjects: ${SUBJECTS[*]}"
echo ">>> Conditions: ${CONDITIONS_CSV:-all available}"
echo ">>> Sensor groups: ${SENSOR_GROUPS_CSV:-all available}"
echo ">>> Include smoothed-display summaries: $INCLUDE_SMOOTHED"
echo ">>> Interval lag-curve shading: $SHADING"
echo ">>> Interval lag-curve mode: $INTERVAL_LAG_MODE"
echo ">>> Lag window sec: $LAG_WINDOW_SEC"
echo ">>> Host: $(hostname)"
echo ">>> Slots: ${NSLOTS:-unknown}"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '>>> Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
