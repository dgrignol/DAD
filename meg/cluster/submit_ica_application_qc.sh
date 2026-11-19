#!/bin/bash
# Audit the whole selected cohort after ICA application, with disjoint workers.
#
# Usage from the cluster DAD root (the user submits jobs):
#   DAD_QC_LABEL=ica-application-qc-NEW DAD_QC_SUBJECTS=04,05,06 \
#     sbatch --cpus-per-task=3 --mem=72G meg/cluster/submit_ica_application_qc.sh
#   DAD_QC_PHASE=epochs DAD_QC_TASK_SIGNAL=0 DAD_QC_LABEL=epoch-audit-NEW \
#     sbatch --cpus-per-task=8 --mem=24G meg/cluster/submit_ica_application_qc.sh
# Keep DAD_EPOCH_PREPARATION/DAD_PROCESSING_RUN exported from epoching in
# epoch mode; the same audits apply to in-memory no-ICA/ICA preparation.
# DAD_QC_PHASE=inputs performs read-only accepted input/Stage 2/3 preflight.
# Inputs: accepted receipts, unchanged prefit models and completed applied FIFs.
# Outputs: new reports/ica_review_updates/LABEL participant audits plus joined
# summary. Default application mode includes full auxiliary/native preservation
# and matched native/ICA-only/final task means; allocate 24 GB per worker.
# DAD_QC_PHASE=inventory only checks prerequisite presence. Epoch mode reconciles
# exact final trial identities and extra MNE drops. No scientific acceptance or
# signal/choice/annotation write is performed. Existing labels are immutable.
# Default subjects explicitly preserve P04-P14/P16-P39; P15 is absent.
#SBATCH --job-name=DAD_ica_qc
#SBATCH --output=cluster_logs/%x.%j.log
#SBATCH --error=cluster_logs/%x.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL

set -euo pipefail

# Resolve the allocated resources and user-selected immutable review namespace.
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
MAMBA_BIN="${MICROMAMBA_BIN:-$HOME/bin/micromamba}"
LABEL="${DAD_QC_LABEL:?Set a new DAD_QC_LABEL before submission.}"
PHASE="${DAD_QC_PHASE:-application}"
TASK_SIGNAL="${DAD_QC_TASK_SIGNAL:-1}"
SUBJECTS="${DAD_QC_SUBJECTS:-04,05,06,07,08,09,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39}"
JOBS="${DAD_QC_JOBS:-${SLURM_CPUS_PER_TASK:-8}}"
[[ "$TASK_SIGNAL" == 0 || "$TASK_SIGNAL" == 1 ]] || { echo "Invalid DAD_QC_TASK_SIGNAL" >&2; exit 2; }
[[ -x "$MAMBA_BIN" ]] || { echo "Missing micromamba" >&2; exit 2; }

# Activate the existing model-compatible environment and bound nested threads.
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "${DAD_MAMBA_ENV:-dad-meg}"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
JOB_TMP="${TMPDIR:-/tmp}/dad_ica_qc_${SLURM_JOB_ID:-manual}"
mkdir -p "$JOB_TMP/matplotlib" "$JOB_TMP/cache"
export MPLCONFIGDIR="$JOB_TMP/matplotlib" XDG_CACHE_HOME="$JOB_TMP/cache"

# Submit all selected participants to the audit's process pool; it joins before
# publishing the shared summary. The audit enforces allocation/cohort/path guards.
cd "$PROJECT_ROOT"
IFS=',' read -r -a SELECTED_SUBJECTS <<< "$SUBJECTS"
COMMAND=(python -u meg/code/preprocessing_basic/audit_ica_application.py
    --site-config meg/params/pipeline_cluster.yaml --subjects "${SELECTED_SUBJECTS[@]}"
    --jobs "$JOBS" --phase "$PHASE" --label "$LABEL")
if [[ "$TASK_SIGNAL" == 1 && "$PHASE" == application ]]; then COMMAND+=(--task-signal); fi
printf '>>> Command:'; printf ' %q' "${COMMAND[@]}"; printf '\n'
"${COMMAND[@]}"
