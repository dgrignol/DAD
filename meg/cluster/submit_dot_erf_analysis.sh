#!/bin/bash
#
# submit_dot_erf_analysis.sh
#
# Purpose
# Run the optional group dot-condition ERF sanity-check analysis on the DIPSEN
# cluster. This wrapper executes ``meg/code/sanity_checks/05_dot_erf_analysis.py``
# as one SLURM array task per Stage-07 sensor group. Each task uses
# the same requested subject group, but writes only one sensor-group-specific
# set of subject/group ERF sanity-check outputs plus linked/embedded HTML
# reports.
#
# Submission examples (run from the DAD root on the DIPSEN login node)
#   mkdir -p cluster_logs
#   sbatch --array=1-4%4 meg/cluster/submit_dot_erf_analysis.sh
# Validate inputs without writing outputs:
#   DAD_VALIDATE_ONLY=1 sbatch meg/cluster/submit_dot_erf_analysis.sh
# Restrict to occipital magnetometers:
#   DAD_REBUILD=1 DAD_DOT_ERF_SENSOR_GROUPS=occipital \
#     sbatch --array=1 meg/cluster/submit_dot_erf_analysis.sh
# Change the analysis window or baseline (quote comma-separated values):
#   DAD_REBUILD=1 DAD_DOT_ERF_TMIN=-0.5 DAD_DOT_ERF_BASELINE='-0.300,-0.100' \
#     sbatch meg/cluster/submit_dot_erf_analysis.sh
# Inputs
# - Completed Stage-04 epoch FIF files under:
#   derivatives/analysis/sub-XX/ses-01/meg/
# - Stage-04 source manifests under:
#   derivatives/reports/sub-XX/ses-01/source_manifest.json
#
# Outputs
# - Per-subject ERF arrays, trial-count CSVs, summaries, and diagnostics under:
#   derivatives/sanity_checks/sub-XX/ses-01/dot_erf/
# - Group ERF arrays, contrasts, cluster stats, figures, and summaries under:
#   derivatives/sanity_checks/group/ses-01/dot_erf/
# - Linked and embedded reports under:
#   derivatives/reports/group/ses-01/dot_erf/
#
# Environment variables
# - DAD_PROJECT_ROOT: cluster DAD root.
# - DAD_MAMBA_ENV, MICROMAMBA_BIN: Python environment selection.
# - DAD_DOT_ERF_GROUP_SUBJECTS: comma-separated subjects; default P04-P14/P16-P39.
# - DAD_DOT_ERF_SENSOR_TYPE: mag or grad; default mag.
# - DAD_DOT_ERF_SENSOR_GROUPS: comma-separated Stage-07 groups; default is
#   occipital,parietal,frontal,temporal. SLURM_ARRAY_TASK_ID selects the Nth
#   group. When the list length changes, pass a matching --array range.
# - DAD_DOT_ERF_BASELINE: comma-separated baseline start,end in seconds;
#   default -0.300,-0.100.
# - DAD_DOT_ERF_TMIN, DAD_DOT_ERF_TMAX: analysis window in seconds. Defaults
#   -0.5..2.667 to match the existing Stage-04 epoch derivatives.
# - DAD_DOT_ERF_ANALYSIS_SFREQ: ERF/statistics sampling rate, default 200.
# - DAD_DOT_ERF_TOPOMAP_TIMES: optional comma-separated topomap latencies.
# - DAD_DOT_ERF_GROUP_OVERLAY_BAND: sem, std, or none; default sem.
# - DAD_CLUSTER_FORMING_P: pointwise p for cluster formation; default 0.05.
# - DAD_CLUSTER_ALPHA: cluster-level significance alpha; default 0.05.
# - DAD_CLUSTER_TAIL: two-sided, greater, or less; default two-sided.
# - DAD_N_PERMUTATIONS: all or positive integer; default 10000 (user approved).
# - Export DAD_EPOCH_PREPARATION and DAD_PROCESSING_RUN as in epoching to
#   select the matching isolated epochs/ERF outputs.
# - DAD_RANDOM_SEED: deterministic permutation/equalization seed; default 42.
# - DAD_DOT_ERF_INCLUDE_CLUSTER_TABLE: 1 includes verbose cluster tables.
# - DAD_DOT_ERF_LOCAL_EPOCH_CACHE_DIR: optional node-visible local epoch cache.
# - DAD_VALIDATE_ONLY: 1 validates inputs and exits before writing outputs.
# - DAD_ALLOW_MISSING: 1 builds from available subjects and lists skipped ones.
# - DAD_REBUILD: 1 overwrites existing ERF outputs; 0 resumes/reuses outputs.

#SBATCH --job-name=DAD_dot_erf
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --export=ALL
#SBATCH --array=1-4%4

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve paths and ERF analysis settings
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required; submit this as an SLURM array job.}"
REBUILD="${DAD_REBUILD:-0}"
VALIDATE_ONLY="${DAD_VALIDATE_ONLY:-0}"
ALLOW_MISSING="${DAD_ALLOW_MISSING:-0}"
SUBJECTS_CSV="${DAD_DOT_ERF_GROUP_SUBJECTS:-04,05,06,07,08,09,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39}"
SENSOR_TYPE="${DAD_DOT_ERF_SENSOR_TYPE:-mag}"
SENSOR_GROUPS_CSV="${DAD_DOT_ERF_SENSOR_GROUPS:-occipital,parietal,frontal,temporal}"
BASELINE_CSV="${DAD_DOT_ERF_BASELINE:--0.300,-0.100}"
TMIN="${DAD_DOT_ERF_TMIN:--0.5}"
TMAX="${DAD_DOT_ERF_TMAX:-2.667}"
ANALYSIS_SFREQ="${DAD_DOT_ERF_ANALYSIS_SFREQ:-200}"
TOPOMAP_TIMES_CSV="${DAD_DOT_ERF_TOPOMAP_TIMES:-}"
GROUP_OVERLAY_BAND="${DAD_DOT_ERF_GROUP_OVERLAY_BAND:-sem}"
CLUSTER_FORMING_P="${DAD_CLUSTER_FORMING_P:-0.05}"
CLUSTER_ALPHA="${DAD_CLUSTER_ALPHA:-0.05}"
CLUSTER_TAIL="${DAD_CLUSTER_TAIL:-two-sided}"
N_PERMUTATIONS="${DAD_N_PERMUTATIONS:-10000}"
RANDOM_SEED="${DAD_RANDOM_SEED:-42}"
INCLUDE_CLUSTER_TABLE="${DAD_DOT_ERF_INCLUDE_CLUSTER_TABLE:-0}"
LOCAL_EPOCH_CACHE_DIR="${DAD_DOT_ERF_LOCAL_EPOCH_CACHE_DIR:-}"

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
if [[ "$INCLUDE_CLUSTER_TABLE" != "0" && "$INCLUDE_CLUSTER_TABLE" != "1" ]]; then
    echo "ERROR: DAD_DOT_ERF_INCLUDE_CLUSTER_TABLE must be 0 or 1; got ${INCLUDE_CLUSTER_TABLE}." >&2
    exit 2
fi
case "$SENSOR_TYPE" in
    mag|grad)
        ;;
    *)
        echo "ERROR: DAD_DOT_ERF_SENSOR_TYPE must be mag or grad; got ${SENSOR_TYPE}." >&2
        exit 2
        ;;
esac
case "$GROUP_OVERLAY_BAND" in
    sem|std|none)
        ;;
    *)
        echo "ERROR: DAD_DOT_ERF_GROUP_OVERLAY_BAND must be sem, std, or none; got ${GROUP_OVERLAY_BAND}." >&2
        exit 2
        ;;
esac
case "$CLUSTER_TAIL" in
    two-sided|greater|less)
        ;;
    *)
        echo "ERROR: DAD_CLUSTER_TAIL must be two-sided, greater, or less; got ${CLUSTER_TAIL}." >&2
        exit 2
        ;;
esac

DEFAULT_CONFIG="$PROJECT_ROOT/meg/params/pipeline_default.yaml"
SITE_CONFIG="$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
SCRIPT="$PROJECT_ROOT/meg/code/sanity_checks/05_dot_erf_analysis.py"
# Multi-value settings accept comma-separated defaults inside this script and
# colon-separated sbatch overrides, because sbatch itself treats commas as the
# boundary between separate environment assignments.
if [[ "$SUBJECTS_CSV" == *:* ]]; then
    IFS=':' read -r -a SUBJECTS <<< "$SUBJECTS_CSV"
else
    IFS=',' read -r -a SUBJECTS <<< "$SUBJECTS_CSV"
fi
if [[ "$SENSOR_GROUPS_CSV" == *:* ]]; then
    IFS=':' read -r -a SENSOR_GROUPS <<< "$SENSOR_GROUPS_CSV"
else
    IFS=',' read -r -a SENSOR_GROUPS <<< "$SENSOR_GROUPS_CSV"
fi
if [[ "$BASELINE_CSV" == *:* ]]; then
    IFS=':' read -r -a BASELINE_VALUES <<< "$BASELINE_CSV"
else
    IFS=',' read -r -a BASELINE_VALUES <<< "$BASELINE_CSV"
fi
if [[ "$TOPOMAP_TIMES_CSV" == *:* ]]; then
    IFS=':' read -r -a TOPOMAP_TIMES <<< "$TOPOMAP_TIMES_CSV"
else
    IFS=',' read -r -a TOPOMAP_TIMES <<< "$TOPOMAP_TIMES_CSV"
fi

if (( ${#SENSOR_GROUPS[@]} == 0 )); then
    echo "ERROR: DAD_DOT_ERF_SENSOR_GROUPS must select at least one sensor group." >&2
    exit 2
fi
if (( TASK_ID < 1 || TASK_ID > ${#SENSOR_GROUPS[@]} )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID ${TASK_ID} does not map to DAD_DOT_ERF_SENSOR_GROUPS=${SENSOR_GROUPS_CSV}." >&2
    echo "Submit with --array=1-${#SENSOR_GROUPS[@]} when overriding the sensor-group list." >&2
    exit 2
fi
if [[ "${#BASELINE_VALUES[@]}" -ne 2 ]]; then
    echo "ERROR: DAD_DOT_ERF_BASELINE must contain exactly two values separated by ':' or ','; got ${BASELINE_CSV}." >&2
    exit 2
fi

SENSOR_GROUP="${SENSOR_GROUPS[$((TASK_ID - 1))]}"


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
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

JOB_TMP="${TMPDIR:-/tmp}/dad_dot_erf_${SLURM_JOB_ID:-manual}_${TASK_ID}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib"
export XDG_CACHE_HOME="$JOB_TMP/cache"


# -----------------------------------------------------------------------------
# Validate deployment and run the ERF analysis
# -----------------------------------------------------------------------------

cd "$PROJECT_ROOT"
mkdir -p cluster_logs

for required_path in \
    "$DEFAULT_CONFIG" \
    "$SITE_CONFIG" \
    "$SCRIPT" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/drsa_result_utils.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/pipeline_config.py"; do
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
    --group-subjects "${SUBJECTS[@]}"
    --sensor-type "$SENSOR_TYPE"
    --sensor-groups "$SENSOR_GROUP"
    --baseline "${BASELINE_VALUES[@]}"
    --tmin "$TMIN"
    --tmax "$TMAX"
    --analysis-sfreq "$ANALYSIS_SFREQ"
    --group-overlay-band "$GROUP_OVERLAY_BAND"
    --cluster-forming-p "$CLUSTER_FORMING_P"
    --cluster-alpha "$CLUSTER_ALPHA"
    --tail "$CLUSTER_TAIL"
    --n-permutations "$N_PERMUTATIONS"
    --random-seed "$RANDOM_SEED"
    --n-jobs "${SLURM_CPUS_PER_TASK:-4}"
    --print-summary
)

if [[ "${#TOPOMAP_TIMES[@]}" -gt 0 && -n "${TOPOMAP_TIMES[0]}" ]]; then
    COMMAND+=(--topomap-times "${TOPOMAP_TIMES[@]}")
fi
if [[ "$ALLOW_MISSING" == "1" ]]; then
    COMMAND+=(--allow-missing)
fi
if [[ "$INCLUDE_CLUSTER_TABLE" == "1" ]]; then
    COMMAND+=(--include-cluster-table)
fi
if [[ -n "$LOCAL_EPOCH_CACHE_DIR" ]]; then
    COMMAND+=(--local-epoch-cache-dir "$LOCAL_EPOCH_CACHE_DIR")
fi

if [[ "$VALIDATE_ONLY" == "1" ]]; then
    COMMAND+=(--validate-only)
    echo ">>> Validate-only mode: no ERF outputs will be written."
elif [[ "$REBUILD" == "1" ]]; then
    COMMAND+=(--overwrite)
    echo ">>> Rebuild mode: dot ERF outputs may be replaced."
else
    COMMAND+=(--resume)
    echo ">>> Resume mode: completed dot ERF outputs are reused when valid."
fi

echo ">>> Subjects: ${SUBJECTS[*]}"
echo ">>> Sensor type: $SENSOR_TYPE"
echo ">>> Sensor group task: $SENSOR_GROUP (${TASK_ID}/${#SENSOR_GROUPS[@]})"
echo ">>> Sensor group list: ${SENSOR_GROUPS[*]}"
echo ">>> Baseline: ${BASELINE_VALUES[*]}"
echo ">>> Time window: $TMIN..$TMAX s"
echo ">>> Analysis sfreq: $ANALYSIS_SFREQ Hz"
echo ">>> Group overlay band: $GROUP_OVERLAY_BAND"
echo ">>> Cluster-forming p: $CLUSTER_FORMING_P"
echo ">>> Cluster alpha: $CLUSTER_ALPHA"
echo ">>> Tail: $CLUSTER_TAIL"
echo ">>> N permutations: $N_PERMUTATIONS"
echo ">>> Include cluster table: $INCLUDE_CLUSTER_TABLE"
echo ">>> Local epoch cache: ${LOCAL_EPOCH_CACHE_DIR:-none}"
echo ">>> Host: $(hostname)"
echo ">>> Allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '>>> Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ">>> Group ERF outputs: $DERIVATIVES_ROOT/sanity_checks/group/ses-01/dot_erf"
echo ">>> Group ERF reports: $DERIVATIVES_ROOT/reports/group/ses-01/dot_erf"
