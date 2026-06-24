#!/bin/bash
#
# deploy_review_to_cluster.sh
#
# Purpose
# Deploy the files required to run DAD preprocessing review/apply, downstream
# epoching/export, MATLAB dRSA, and optional sanity checks on the CIMeC
# cluster. This script runs on the
# local computer, creates the remote DAD directory outside the raw tree, and
# preserves the required ``meg/...`` layout with rsync.
#
# Usage
#   cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD
#   meg/cluster/deploy_review_to_cluster.sh
#
# Optional overrides
#   DAD_CLUSTER_HOST=user@cluster.example \
#   DAD_CLUSTER_ROOT=/remote/path/DAD \
#     meg/cluster/deploy_review_to_cluster.sh
#
# Inputs
# - Production/infrastructure Python files through dRSA export and MATLAB dRSA.
# - Optional sanity-check Python files.
# - Shared and subject parameter YAMLs.
# - Review/apply cluster submission, documentation, and environment files.
# - Maxwell cross-talk and calibration resources.
#
# Outputs
# - A minimal working copy at the configured cluster destination.
#
# Safety
# The destination defaults to ``Attention/DAD``, which is outside the
# read-only ``Attention/data`` tree. The script refuses a destination that is
# equal to or nested below a path component named ``data``. It does not use
# rsync deletion flags and therefore cannot remove unrelated remote files.

set -euo pipefail


# -----------------------------------------------------------------------------
# Resolve local and remote roots
# -----------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
REMOTE_HOST="${DAD_CLUSTER_HOST:-damiano.grignolio@cluster.cimec.unitn.it}"
REMOTE_ROOT="${DAD_CLUSTER_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"

case "/$REMOTE_ROOT/" in
    */data/*)
        echo "ERROR: refusing to deploy inside a path component named data:" >&2
        echo "  $REMOTE_ROOT" >&2
        exit 2
        ;;
esac


# -----------------------------------------------------------------------------
# Enumerate the minimal production-pipeline payload
# -----------------------------------------------------------------------------

FILES=(
    "meg/code/analysis/"
    "meg/code/preprocessing_basic/02_trigger_extraction_custom.py"
    "meg/code/preprocessing_basic/03_A_MaxFilt.py"
    "meg/code/preprocessing_basic/03_B_ArtefactAnnotation.py"
    "meg/code/preprocessing_basic/03_C_ICA.py"
    "meg/code/preprocessing_basic/04_epoching.py"
    "meg/code/preprocessing_basic/06_export_neural_for_drsa.py"
    "meg/code/preprocessing_basic/07_run_matlab_drsa.py"
    "meg/code/preprocessing_basic/08_plot_drsa_line_summaries.py"
    "meg/code/preprocessing_basic/dataset_inventory.py"
    "meg/code/preprocessing_basic/downstream_stage_utils.py"
    "meg/code/preprocessing_basic/drsa_canonical_meg_channels.txt"
    "meg/code/preprocessing_basic/drsa_result_utils.py"
    "meg/code/preprocessing_basic/pipeline_config.py"
    "meg/code/preprocessing_basic/preprocessing_stage_utils.py"
    "meg/code/preprocessing_basic/run_pipeline.py"
    "meg/code/preprocessing_basic/trigger_extraction.py"
    "meg/code/sanity_checks/"
    "meg/params/pipeline_default.yaml"
    "meg/params/pipeline_cluster.yaml"
    "meg/params/sub-01.yaml"
    "meg/params/sub-02.yaml"
    "meg/params/sub-03.yaml"
    "meg/params/sub-04.yaml"
    "meg/params/sub-05.yaml"
    "meg/params/sub-06.yaml"
    "meg/params/sub-07.yaml"
    "meg/cluster/README.md"
    "meg/cluster/environment.yml"
    "meg/cluster/requirements.txt"
    "meg/cluster/submit_preprocessing_review.sh"
    "meg/cluster/submit_preprocessing_apply.sh"
    "meg/cluster/submit_epoch_export.sh"
    "meg/cluster/submit_drsa_analysis.sh"
    "meg/cluster/submit_drsa_line_summaries.sh"
    "meg/cluster/submit_dot_side_decoding.sh"
    "meg/cluster/submit_dot_side_group_report.sh"
    "meg/cluster/deploy_review_to_cluster.sh"
    "meg/data/fwmegmaxfilter/ct_sparse_triux2.fif"
    "meg/data/fwmegmaxfilter/sss_cal_3182_250410.dat"
    "simulations/functions/"
)

for relative_path in "${FILES[@]}"; do
    if [[ ! -e "$LOCAL_ROOT/$relative_path" ]]; then
        echo "ERROR: required local file is missing: $relative_path" >&2
        exit 2
    fi
done


# -----------------------------------------------------------------------------
# Create the destination and transfer while preserving relative paths
# -----------------------------------------------------------------------------

echo ">>> Local root: $LOCAL_ROOT"
echo ">>> Remote: $REMOTE_HOST:$REMOTE_ROOT"
echo ">>> Files: ${#FILES[@]}"

ssh "$REMOTE_HOST" "mkdir -p '$REMOTE_ROOT'"
(
    cd "$LOCAL_ROOT"
    rsync -avR \
        --exclude '.DS_Store' \
        --exclude '__pycache__/' \
        "${FILES[@]}" \
        "$REMOTE_HOST:$REMOTE_ROOT/"
)

echo ">>> Minimal cluster production working copy deployed."
