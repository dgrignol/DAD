#!/bin/bash
#
# submit_drsa_stats.sh
#
# Purpose
# Build the optional Stage-10 group dRSA statistical report on the CIMeC
# cluster from existing subject-level Stage-07 result MAT files. This is a
# one-task group aggregation job; it does not run MATLAB and does not rerun
# Stage 04, Stage 06, Stage 07, Stage 08, or Stage 09.
#
# Submission examples (run from the cluster head node)
#   cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
#   mkdir -p cluster_logs
#   qsub meg/cluster/submit_drsa_stats.sh
#
# Rebuild existing Stage-10 stats outputs:
#   qsub -v DAD_REBUILD=1 meg/cluster/submit_drsa_stats.sh
#
# Restrict to occipital A_dRSA_alt line-summary stats:
#   qsub -v DAD_DRSA_STATS_TARGETS=line-summaries,DAD_DRSA_STATS_ANALYSIS_FAMILIES=A_dRSA_alt,DAD_DRSA_STATS_SENSOR_GROUPS=occipital \
#     meg/cluster/submit_drsa_stats.sh
#
# Include display-smoothed inputs as explicitly labeled additional tests:
#   qsub -v DAD_REBUILD=1,DAD_DRSA_STATS_INCLUDE_SMOOTHED=1 \
#     meg/cluster/submit_drsa_stats.sh
#
# Inputs
# - Existing Stage-07 result MAT files below
#   derivatives/analysis/sub-XX/ses-01/drsa/.
#
# Outputs
# - Statistical arrays, figures, and machine-readable tables below
#   derivatives/analysis/group/ses-01/drsa/stats/.
# - Linked and embedded HTML reports below
#   derivatives/reports/group/ses-01/.
#
# Environment variables
# - DAD_PROJECT_ROOT: cluster DAD root.
# - DAD_MAMBA_ENV, MICROMAMBA_BIN: Python environment selection.
# - DAD_DRSA_STATS_SUBJECTS: comma-separated subjects, default 01..07.
# - DAD_DRSA_STATS_CONDITIONS: optional comma-separated condition filter.
# - DAD_DRSA_STATS_MODELS: optional comma-separated model/feature selectors.
# - DAD_DRSA_STATS_MODEL_SOURCES: optional comma-separated model-source filter.
# - DAD_DRSA_STATS_SENSOR_GROUPS: optional comma-separated sensor-group filter.
# - DAD_DRSA_STATS_ANALYSIS_FAMILIES: comma-separated family filter,
#   default A_dRSA,A_dRSA_alt.
# - DAD_DRSA_STATS_TARGETS: comma-separated matrices,line-summaries by default.
# - DAD_DRSA_STATS_LINE_INTERVALS: comma-separated interval filter.
# - DAD_DRSA_STATS_LINE_LAG_WINDOW_SEC: line-summary lag half-window, default 1.0.
# - DAD_DRSA_STATS_LINE_INTERVAL_LAG_MODE: model-anchored or same-interval-square.
# - DAD_DRSA_STATS_INCLUDE_SMOOTHED: 1 also tests smoothed-display inputs.
# - DAD_DRSA_STATS_PLOT_SMOOTHING_WINDOW: optional smoothing window in samples.
# - DAD_DRSA_STATS_N_PERMUTATIONS: sign-flip permutations, default 1024.
# - DAD_DRSA_STATS_CLUSTER_ALPHA: cluster-forming p-value, default 0.05.
# - DAD_DRSA_STATS_SIGNIFICANCE_ALPHA: cluster p threshold, default 0.05.
# - DAD_DRSA_STATS_TAIL: two-sided, positive, or negative.
# - DAD_DRSA_STATS_SEED: random seed, default 97.
# - DAD_DRSA_STATS_BACKEND: auto, mne, or numpy.
# - DAD_DRSA_STATS_PAIRED_CONTRASTS: none or actual_deviant_position_vs_nondeviant.
# - DAD_DRSA_STATS_REPORT_CLUSTER_TABLES: significant, all, or none.
# - DAD_REBUILD: 1 overwrites existing outputs; 0 resumes existing outputs.

#$ -N DAD_drsa_stats
#$ -S /bin/bash
#$ -cwd
#$ -o cluster_logs/$JOB_NAME.$JOB_ID.log
#$ -j y
#$ -pe smp 4
#$ -l h_vmem=16G

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve paths and Stage-10 settings
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
REBUILD="${DAD_REBUILD:-0}"
SUBJECTS_CSV="${DAD_DRSA_STATS_SUBJECTS:-01,02,03,04,05,06,07}"
CONDITIONS_CSV="${DAD_DRSA_STATS_CONDITIONS:-}"
MODELS_CSV="${DAD_DRSA_STATS_MODELS:-}"
MODEL_SOURCES_CSV="${DAD_DRSA_STATS_MODEL_SOURCES:-}"
SENSOR_GROUPS_CSV="${DAD_DRSA_STATS_SENSOR_GROUPS:-}"
ANALYSIS_FAMILIES_CSV="${DAD_DRSA_STATS_ANALYSIS_FAMILIES:-A_dRSA,A_dRSA_alt}"
TARGETS_CSV="${DAD_DRSA_STATS_TARGETS:-matrices,line-summaries}"
LINE_INTERVALS_CSV="${DAD_DRSA_STATS_LINE_INTERVALS:-pre_occlusion,during_occlusion,post_occlusion}"
LINE_LAG_WINDOW_SEC="${DAD_DRSA_STATS_LINE_LAG_WINDOW_SEC:-1.0}"
LINE_INTERVAL_LAG_MODE="${DAD_DRSA_STATS_LINE_INTERVAL_LAG_MODE:-model-anchored}"
INCLUDE_SMOOTHED="${DAD_DRSA_STATS_INCLUDE_SMOOTHED:-0}"
PLOT_SMOOTHING_WINDOW="${DAD_DRSA_STATS_PLOT_SMOOTHING_WINDOW:-}"
N_PERMUTATIONS="${DAD_DRSA_STATS_N_PERMUTATIONS:-1024}"
CLUSTER_ALPHA="${DAD_DRSA_STATS_CLUSTER_ALPHA:-0.05}"
SIGNIFICANCE_ALPHA="${DAD_DRSA_STATS_SIGNIFICANCE_ALPHA:-0.05}"
TAIL="${DAD_DRSA_STATS_TAIL:-two-sided}"
SEED="${DAD_DRSA_STATS_SEED:-97}"
STATS_BACKEND="${DAD_DRSA_STATS_BACKEND:-auto}"
PAIRED_CONTRASTS_CSV="${DAD_DRSA_STATS_PAIRED_CONTRASTS:-actual_deviant_position_vs_nondeviant}"
REPORT_CLUSTER_TABLES="${DAD_DRSA_STATS_REPORT_CLUSTER_TABLES:-significant}"

if [[ "$REBUILD" != "0" && "$REBUILD" != "1" ]]; then
    echo "ERROR: DAD_REBUILD must be 0 or 1; got ${REBUILD}." >&2
    exit 2
fi
if [[ "$INCLUDE_SMOOTHED" != "0" && "$INCLUDE_SMOOTHED" != "1" ]]; then
    echo "ERROR: DAD_DRSA_STATS_INCLUDE_SMOOTHED must be 0 or 1; got ${INCLUDE_SMOOTHED}." >&2
    exit 2
fi
case "$LINE_INTERVAL_LAG_MODE" in
    model-anchored|same-interval-square)
        ;;
    *)
        echo "ERROR: DAD_DRSA_STATS_LINE_INTERVAL_LAG_MODE must be model-anchored or same-interval-square; got ${LINE_INTERVAL_LAG_MODE}." >&2
        exit 2
        ;;
esac
case "$TAIL" in
    two-sided|positive|negative)
        ;;
    *)
        echo "ERROR: DAD_DRSA_STATS_TAIL must be two-sided, positive, or negative; got ${TAIL}." >&2
        exit 2
        ;;
esac
case "$STATS_BACKEND" in
    auto|mne|numpy)
        ;;
    *)
        echo "ERROR: DAD_DRSA_STATS_BACKEND must be auto, mne, or numpy; got ${STATS_BACKEND}." >&2
        exit 2
        ;;
esac
case "$REPORT_CLUSTER_TABLES" in
    significant|all|none)
        ;;
    *)
        echo "ERROR: DAD_DRSA_STATS_REPORT_CLUSTER_TABLES must be significant, all, or none; got ${REPORT_CLUSTER_TABLES}." >&2
        exit 2
        ;;
esac

DEFAULT_CONFIG="$PROJECT_ROOT/meg/params/pipeline_default.yaml"
SITE_CONFIG="$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
SCRIPT="$PROJECT_ROOT/meg/code/preprocessing_basic/10_drsa_stats.py"
IFS=',' read -r -a SUBJECTS <<< "$SUBJECTS_CSV"
IFS=',' read -r -a TARGETS <<< "$TARGETS_CSV"
IFS=',' read -r -a LINE_INTERVALS <<< "$LINE_INTERVALS_CSV"
IFS=',' read -r -a PAIRED_CONTRASTS <<< "$PAIRED_CONTRASTS_CSV"


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
export OMP_NUM_THREADS="${NSLOTS:-4}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-4}"
export MKL_NUM_THREADS="${NSLOTS:-4}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-4}"

JOB_TMP="${TMPDIR:-/tmp}/dad_drsa_stats_${JOB_ID:-manual}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib"
export XDG_CACHE_HOME="$JOB_TMP/cache"


# -----------------------------------------------------------------------------
# Validate deployment and build the group stats report
# -----------------------------------------------------------------------------

cd "$PROJECT_ROOT"
mkdir -p cluster_logs

for required_path in \
    "$DEFAULT_CONFIG" \
    "$SITE_CONFIG" \
    "$SCRIPT" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/08_plot_drsa_line_summaries.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/drsa_result_utils.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/pipeline_config.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/preprocessing_stage_utils.py"; do
    if [[ ! -f "$required_path" ]]; then
        echo "ERROR: required file is missing: $required_path" >&2
        exit 2
    fi
done

python - <<'PY'
import importlib.util
import numpy
import scipy
print(f">>> Python stats deps: numpy={numpy.__version__}, scipy={scipy.__version__}")
for name in ("mne", "h5py", "matplotlib"):
    spec = importlib.util.find_spec(name)
    print(f">>> Optional dep {name}: {'available' if spec else 'missing'}")
PY

COMMAND=(
    python
    "$SCRIPT"
    --config "$DEFAULT_CONFIG"
    --site-config "$SITE_CONFIG"
    --subjects "${SUBJECTS[@]}"
    --stats-targets "${TARGETS[@]}"
    --line-intervals "${LINE_INTERVALS[@]}"
    --line-lag-window-sec "$LINE_LAG_WINDOW_SEC"
    --line-interval-lag-mode "$LINE_INTERVAL_LAG_MODE"
    --n-permutations "$N_PERMUTATIONS"
    --cluster-alpha "$CLUSTER_ALPHA"
    --significance-alpha "$SIGNIFICANCE_ALPHA"
    --tail "$TAIL"
    --seed "$SEED"
    --stats-backend "$STATS_BACKEND"
    --paired-contrasts "${PAIRED_CONTRASTS[@]}"
    --report-cluster-tables "$REPORT_CLUSTER_TABLES"
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

if [[ -n "$MODELS_CSV" ]]; then
    IFS=',' read -r -a MODELS <<< "$MODELS_CSV"
    COMMAND+=(--models "${MODELS[@]}")
fi

if [[ -n "$MODEL_SOURCES_CSV" ]]; then
    IFS=',' read -r -a MODEL_SOURCES <<< "$MODEL_SOURCES_CSV"
    COMMAND+=(--model-sources "${MODEL_SOURCES[@]}")
fi

if [[ -n "$SENSOR_GROUPS_CSV" ]]; then
    IFS=',' read -r -a SENSOR_GROUPS <<< "$SENSOR_GROUPS_CSV"
    COMMAND+=(--sensor-groups "${SENSOR_GROUPS[@]}")
fi

if [[ -n "$ANALYSIS_FAMILIES_CSV" ]]; then
    IFS=',' read -r -a ANALYSIS_FAMILIES <<< "$ANALYSIS_FAMILIES_CSV"
    COMMAND+=(--analysis-families "${ANALYSIS_FAMILIES[@]}")
fi

if [[ "$REBUILD" == "1" ]]; then
    COMMAND+=(--overwrite)
    echo ">>> Rebuild mode: Stage-10 stats outputs may be replaced."
else
    COMMAND+=(--resume)
    echo ">>> Resume mode: existing Stage-10 stats outputs are reused."
fi

echo ">>> Subjects: ${SUBJECTS[*]}"
echo ">>> Stats targets: ${TARGETS[*]}"
echo ">>> Conditions: ${CONDITIONS_CSV:-all available}"
echo ">>> Models: ${MODELS_CSV:-position,prediction_error default}"
echo ">>> Model sources: ${MODEL_SOURCES_CSV:-all selected}"
echo ">>> Sensor groups: ${SENSOR_GROUPS_CSV:-all available}"
echo ">>> Analysis families: ${ANALYSIS_FAMILIES_CSV}"
echo ">>> Paired contrasts: ${PAIRED_CONTRASTS[*]}"
echo ">>> HTML cluster tables: $REPORT_CLUSTER_TABLES"
echo ">>> Line intervals: ${LINE_INTERVALS[*]}"
echo ">>> Line lag window sec: $LINE_LAG_WINDOW_SEC"
echo ">>> Line interval lag mode: $LINE_INTERVAL_LAG_MODE"
echo ">>> Include smoothed-display inputs: $INCLUDE_SMOOTHED"
echo ">>> N permutations: $N_PERMUTATIONS"
echo ">>> Cluster alpha: $CLUSTER_ALPHA"
echo ">>> Significance alpha: $SIGNIFICANCE_ALPHA"
echo ">>> Tail: $TAIL"
echo ">>> Backend: $STATS_BACKEND"
echo ">>> Host: $(hostname)"
echo ">>> Slots: ${NSLOTS:-unknown}"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '>>> Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
