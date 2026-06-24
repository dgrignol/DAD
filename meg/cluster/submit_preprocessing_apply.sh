#!/bin/bash
#
# submit_preprocessing_apply.sh
#
# Purpose
# Submit the reviewed DAD stage-03 apply pass as one Sun Grid Engine array task
# per participant. This wrapper requires each subject's
# ``03_tuning_choices.yaml`` and intentionally overwrites review scientific
# derivatives with common-destination Maxwell outputs, reviewed artifact
# thresholds, and manual-ICA-cleaned FIF files.
#
# Usage
#   cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
#   mkdir -p cluster_logs
#   qsub -t 1-7 meg/cluster/submit_preprocessing_apply.sh
#
# Select only specific participants:
#   qsub -t 2,4,7 meg/cluster/submit_preprocessing_apply.sh
#
# Inputs
# - The same raw/configuration inputs as submit_preprocessing_review.sh.
# - A reviewed 03_tuning_choices.yaml for every submitted participant.
# - ``ica.manual_exclude`` must be a list; ``[]`` is valid and means reviewed
#   with no component rejection.
#
# Outputs
# - Rebuilt *_raw_sss.fif and *_ann.fif derivatives.
# - Final manual-ICA-cleaned *_ica.fif derivatives and apply reports.
#
# Safety
# Raw data remain read-only. The delegated runner uses the shared raw-tree
# guard. Overwrite is explicit because apply must replace review scientific
# derivatives rather than resume outputs created with different decisions.
#
# Grid Engine execution detail
# SGE copies the submitted script into a per-job spool directory before
# execution, so BASH_SOURCE points to that temporary copy. The shared review
# wrapper must therefore be resolved from DAD_PROJECT_ROOT rather than from
# the apply script's runtime directory.

#$ -N DAD_meg_apply
#$ -S /bin/bash
#$ -cwd
#$ -o cluster_logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -j y
#$ -pe smp 4
#$ -l h_vmem=24G
#$ -t 1-7

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

export DAD_PROCESSING_PASS="apply"
export DAD_REBUILD="1"

exec "$REVIEW_WRAPPER"
