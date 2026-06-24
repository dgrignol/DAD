#!/bin/bash
#
# submit_dot_side_time_generalization.sh
#
# Purpose
# Run the optional dot-side time-generalization sanity check on the CIMeC
# cluster as one Sun Grid Engine array task per participant. Each task runs
# the four standard neural-condition/model-source analyses sequentially. This
# wrapper consumes completed Stage-04 epochs and the recorded Stage-04 source
# manifest; it never reruns epoching, dRSA export, MATLAB dRSA, or
# preprocessing.
#
# Submission examples (run from the cluster head node)
#   cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
#   mkdir -p cluster_logs
#   qsub -t 1-7 meg/cluster/submit_dot_side_time_generalization.sh
#
# Rebuild existing time-generalization outputs intentionally:
#   qsub -v DAD_REBUILD=1 -t 1-7 \
#     meg/cluster/submit_dot_side_time_generalization.sh
#
# Run only selected analyses and use identity-level pseudo-trials:
#   export DAD_DOT_SIDE_TG_ANALYSES="occluded_deviant_observed,occluded_deviant_predicted"
#   qsub -v DAD_DOT_SIDE_TG_ANALYSES,DAD_TRIALS_PER_AVERAGE=all_trials_per_identity \
#     -t 1-7 meg/cluster/submit_dot_side_time_generalization.sh
#
# Validate inputs only, without decoding or writing outputs:
#   qsub -v DAD_VALIDATE_ONLY=1 -t 1,4 \
#     meg/cluster/submit_dot_side_time_generalization.sh
#
# Inputs
# - Completed Stage-04 epochs FIF under derivatives/analysis/sub-XX/ses-01/meg.
# - source_manifest.json with completed Stage-04 epoching lineage.
# - Subject-specific trajectory MAT file resolved by configuration.
# - Shared/site/subject configuration files.
#
# Outputs
# - Heatmap PNG, matrix NPZ, long matrix CSV, summary JSON, and MNE HTML report
#   under derivatives/sanity_checks/sub-XX/ses-01/
#   dot_side_time_generalization/<condition>/<model-source>/.
#
# Environment variables
# - DAD_PROJECT_ROOT: cluster DAD root.
# - DAD_MAMBA_ENV, MICROMAMBA_BIN: Python environment selection.
# - DAD_REBUILD: 0 resumes complete outputs; 1 overwrites outputs.
# - DAD_VALIDATE_ONLY: 1 runs validation and exits before decoding.
# - DAD_DOT_SIDE_TG_ANALYSES: comma-separated analysis keys; default is all.
#   DAD_DOT_SIDE_ANALYSES is accepted as a fallback for parity with the 1D job.
# - DAD_DECODE_TMIN, DAD_DECODE_TMAX, DAD_DECODE_SFREQ: decoding time grid.
# - DAD_DECODE_WINDOW_SAMPLES: centered sliding-window width, default 10.
# - DAD_TRIALS_PER_AVERAGE: positive integer or all_trials_per_identity.
# - DAD_RANDOM_SEED: deterministic averaging/CV seed.
# - DAD_N_SPLITS: stratified CV folds.
# - DAD_DOT_SIDE_TG_OUT_DIR: optional explicit output root.
# Resource override example:
#   qsub -pe smp 8 -l h_vmem=48G -t 1-7 \
#     meg/cluster/submit_dot_side_time_generalization.sh

#$ -N DAD_dot_side_tg
#$ -S /bin/bash
#$ -cwd
#$ -o cluster_logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -j y
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -t 1-7

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve subject, paths, and runtime settings
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
DERIVATIVES_ROOT="${DAD_DERIVATIVES_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives}"
MAMBA_ENV="${DAD_MAMBA_ENV:-dad-meg}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
TASK_ID="${SGE_TASK_ID:?SGE_TASK_ID is required; submit this as an SGE array job.}"
REBUILD="${DAD_REBUILD:-0}"
VALIDATE_ONLY="${DAD_VALIDATE_ONLY:-0}"
DECODE_TMIN="${DAD_DECODE_TMIN:-0.0}"
DECODE_TMAX="${DAD_DECODE_TMAX:-2.667}"
DECODE_SFREQ="${DAD_DECODE_SFREQ:-120.0}"
DECODE_WINDOW_SAMPLES="${DAD_DECODE_WINDOW_SAMPLES:-10}"
TRIALS_PER_AVERAGE="${DAD_TRIALS_PER_AVERAGE:-1}"
RANDOM_SEED="${DAD_RANDOM_SEED:-42}"
N_SPLITS="${DAD_N_SPLITS:-5}"

if [[ -n "${DAD_DOT_SIDE_TG_ANALYSES:-}" ]]; then
    ANALYSES_RAW="$DAD_DOT_SIDE_TG_ANALYSES"
elif [[ -n "${DAD_DOT_SIDE_ANALYSES:-}" ]]; then
    ANALYSES_RAW="$DAD_DOT_SIDE_ANALYSES"
else
    ANALYSES_RAW="all"
fi

if [[ "$ANALYSES_RAW" == "all" ]]; then
    ANALYSES=(
        "always_visible_observed"
        "occluded_nondeviant_observed"
        "occluded_deviant_observed"
        "occluded_deviant_predicted"
    )
else
    IFS=',' read -r -a ANALYSES <<< "$ANALYSES_RAW"
fi

if (( ${#ANALYSES[@]} == 0 )); then
    echo "ERROR: no dot-side time-generalization analyses were selected." >&2
    exit 2
fi
for analysis_spec in "${ANALYSES[@]}"; do
    case "$analysis_spec" in
        always_visible_observed|occluded_nondeviant_observed|occluded_deviant_observed|occluded_deviant_predicted)
            ;;
        *)
            echo "ERROR: unsupported dot-side analysis: $analysis_spec" >&2
            exit 2
            ;;
    esac
done

if (( TASK_ID < 1 || TASK_ID > 7 )); then
    echo "ERROR: SGE_TASK_ID must be between 1 and 7; got ${TASK_ID}." >&2
    exit 2
fi
if [[ "$REBUILD" != "0" && "$REBUILD" != "1" ]]; then
    echo "ERROR: DAD_REBUILD must be 0 or 1; got ${REBUILD}." >&2
    exit 2
fi
if [[ "$VALIDATE_ONLY" != "0" && "$VALIDATE_ONLY" != "1" ]]; then
    echo "ERROR: DAD_VALIDATE_ONLY must be 0 or 1; got ${VALIDATE_ONLY}." >&2
    exit 2
fi
if [[ "$TRIALS_PER_AVERAGE" != "all_trials_per_identity" && ! "$TRIALS_PER_AVERAGE" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: DAD_TRIALS_PER_AVERAGE must be a positive integer or all_trials_per_identity; got ${TRIALS_PER_AVERAGE}." >&2
    exit 2
fi

SUBJECT_ID="$(printf "%02d" "$TASK_ID")"
SUBJECT_LABEL="sub-${SUBJECT_ID}"
DEFAULT_CONFIG="$PROJECT_ROOT/meg/params/pipeline_default.yaml"
SITE_CONFIG="$PROJECT_ROOT/meg/params/pipeline_cluster.yaml"
SUBJECT_CONFIG="$PROJECT_ROOT/meg/params/${SUBJECT_LABEL}.yaml"
REPORT_DIR="$DERIVATIVES_ROOT/reports/${SUBJECT_LABEL}/ses-01"
SOURCE_MANIFEST="$REPORT_DIR/source_manifest.json"
SCRIPT="$PROJECT_ROOT/meg/code/sanity_checks/03_dot_side_time_generalization.py"


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

JOB_TMP="${TMPDIR:-/tmp}/dad_dot_side_tg_${JOB_ID:-manual}_${TASK_ID}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib"
export XDG_CACHE_HOME="$JOB_TMP/cache"


# -----------------------------------------------------------------------------
# Validate deployment and launch the optional sanity check
# -----------------------------------------------------------------------------

cd "$PROJECT_ROOT"
mkdir -p cluster_logs

for required_path in \
    "$DEFAULT_CONFIG" \
    "$SITE_CONFIG" \
    "$SUBJECT_CONFIG" \
    "$SOURCE_MANIFEST" \
    "$SCRIPT" \
    "$PROJECT_ROOT/meg/code/sanity_checks/01_dot_side_decoding.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/downstream_stage_utils.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/drsa_result_utils.py" \
    "$PROJECT_ROOT/meg/code/preprocessing_basic/pipeline_config.py"; do
    if [[ ! -f "$required_path" ]]; then
        echo "ERROR: required file is missing: $required_path" >&2
        exit 2
    fi
done

echo ">>> Subject: $SUBJECT_LABEL"
echo ">>> Analyses: ${ANALYSES[*]}"
echo ">>> Decode grid: ${DECODE_TMIN}..${DECODE_TMAX} s @ ${DECODE_SFREQ} Hz"
echo ">>> Decode window samples: $DECODE_WINDOW_SAMPLES"
echo ">>> Trials per average: $TRIALS_PER_AVERAGE"
echo ">>> Random seed: $RANDOM_SEED"
echo ">>> Host: $(hostname)"
echo ">>> Slots: ${NSLOTS:-unknown}"
echo ">>> Started UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"

for analysis_spec in "${ANALYSES[@]}"; do
    MODEL_SOURCE="${analysis_spec##*_}"
    CONDITION="${analysis_spec%_*}"
    COMMAND=(
        python
        "$SCRIPT"
        --config "$DEFAULT_CONFIG"
        --site-config "$SITE_CONFIG"
        --subject-config "$SUBJECT_CONFIG"
        --source-manifest "$SOURCE_MANIFEST"
        --subject "$SUBJECT_ID"
        --condition "$CONDITION"
        --model-source "$MODEL_SOURCE"
        --decode-tmin "$DECODE_TMIN"
        --decode-tmax "$DECODE_TMAX"
        --decode-sfreq "$DECODE_SFREQ"
        --decode-window-samples "$DECODE_WINDOW_SAMPLES"
        --trials-per-average "$TRIALS_PER_AVERAGE"
        --random-seed "$RANDOM_SEED"
        --n-splits "$N_SPLITS"
    )

    if [[ -n "${DAD_DOT_SIDE_TG_OUT_DIR:-}" ]]; then
        COMMAND+=(--out-dir "$DAD_DOT_SIDE_TG_OUT_DIR/$CONDITION/$MODEL_SOURCE")
    fi

    if [[ "$VALIDATE_ONLY" == "1" ]]; then
        COMMAND+=(--validate-only --print-summary)
    elif [[ "$REBUILD" == "1" ]]; then
        COMMAND+=(--overwrite)
    else
        COMMAND+=(--resume)
    fi

    echo ">>> Analysis: condition=$CONDITION model_source=$MODEL_SOURCE"
    printf '>>> Command:'
    printf ' %q' "${COMMAND[@]}"
    printf '\n'
    "${COMMAND[@]}"
done

echo ">>> Completed UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ">>> Reports: $REPORT_DIR"
echo ">>> Time-generalization outputs: $DERIVATIVES_ROOT/sanity_checks/$SUBJECT_LABEL/ses-01/dot_side_time_generalization"
