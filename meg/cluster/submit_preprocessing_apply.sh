#!/bin/bash
#
# submit_preprocessing_apply.sh
#
# Purpose
# Apply the exact reviewed final ICA model per P04-P39 participant, reusing
# accepted 4A/4B outputs. Never rebuild those upstream inputs after ICA review.
#
# Usage from remote DAD root after accepted component/filter decisions
#   sbatch --array=4-14,16-39 meg/cluster/submit_preprocessing_apply.sh
#   sbatch --array=4,7,24 meg/cluster/submit_preprocessing_apply.sh
# Requires 03_tuning_choices.yaml, accepted 4A/4B/4C evidence, model-bound
# exclusions (explicit [] is valid), and reviewed final spectrum-fit evidence.
# No FIR exposure margin is needed with the selected null high-pass. P38 run 1
# also receives spectrum-fit; its historical exception protects a future HP only.
# Data flow: accepted native annotated runs + exact prefit model -> native
# ICA correction -> ICA-only comparisons -> final temporal filtering -> guarded
# continuous outputs/QC. Raw data and Stage 1-3 evidence remain unchanged.
# Existing selected final outputs are explicitly overwritten. For interruption
# recovery use run_pipeline.py --stages ica_manual --processing-pass apply
# --resume instead; upstream outputs must remain identical.
#
#SBATCH --job-name=DAD_meg_apply
#SBATCH --output=cluster_logs/%x.%A_%a.log
#SBATCH --error=cluster_logs/%x.%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --export=ALL
#SBATCH --array=4-14,16-39

set -euo pipefail


# -----------------------------------------------------------------------------
# Delegate environment activation and execution to the reviewed shared wrapper
# -----------------------------------------------------------------------------

PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
REVIEW_WRAPPER="$PROJECT_ROOT/meg/cluster/submit_preprocessing_review.sh"

if [[ ! -x "$REVIEW_WRAPPER" ]]; then
    echo "ERROR: preprocessing review wrapper is missing or not executable:" >&2
    echo "  $REVIEW_WRAPPER" >&2
    exit 2
fi

# Final application reuses accepted SSS/annotation outputs and the exact reviewed
# model. It must never rebuild 4A/4B after component indices have been chosen.
export DAD_STAGES="${DAD_STAGES:-ica_manual}"
if [[ "$DAD_STAGES" != "ica_manual" ]]; then
    echo "ERROR: final apply selects ica_manual only; rerun upstream review explicitly." >&2
    exit 2
fi
export DAD_PROCESSING_PASS="apply"
export DAD_REBUILD="1"

exec "$REVIEW_WRAPPER"
