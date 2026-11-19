#!/bin/bash
#
# deploy_review_to_cluster.sh
#
# Purpose
# Deploy the files required to run Stage 2 trial accounting, Stage 3 timing, preprocessing
# review/apply, downstream
# epoching/export, MATLAB dRSA, optional dRSA post-processing, and optional
# sanity checks on the DIPSEN HPC cluster. This script runs on the
# local computer, creates the remote DAD directory outside the raw tree, and
# preserves the required ``meg/...`` layout with rsync.
#
# Usage
#   cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD
#   meg/cluster/deploy_review_to_cluster.sh
#
# Optional overrides
#   DAD_CLUSTER_HOST=damiano.grignolio@hpc.dipsen.unitn.it \
#   DAD_CLUSTER_ROOT=/remote/path/DAD \
#     meg/cluster/deploy_review_to_cluster.sh
#
# Inputs
# - Production/infrastructure Python files through dRSA export, MATLAB dRSA,
#   and dRSA post-processing.
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
REMOTE_HOST="${DAD_CLUSTER_HOST:-damiano.grignolio@hpc.dipsen.unitn.it}"
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
    "meg/code/preprocessing_basic/09_plot_drsa_sliding_lag_surfaces.py"
    "meg/code/preprocessing_basic/10_drsa_stats.py"
    "meg/code/preprocessing_basic/dataset_inventory.py"
    "meg/code/preprocessing_basic/deep_trigger_report.py"
    "meg/code/preprocessing_basic/downstream_stage_utils.py"
    "meg/code/preprocessing_basic/trial_preparation.py"
    "meg/code/preprocessing_basic/visual_timing.py"
    "meg/code/preprocessing_basic/visual_timing_utils.py"
    "meg/code/preprocessing_basic/visual_timing_report.py"
    "meg/code/preprocessing_basic/visual_timing_application.py"
    "meg/code/preprocessing_basic/visual_timing_task2.py"
    "meg/code/preprocessing_basic/drsa_canonical_meg_channels.txt"
    "meg/code/preprocessing_basic/drsa_result_utils.py"
    "meg/code/preprocessing_basic/pipeline_config.py"
    "meg/code/preprocessing_basic/p24_acquisition_adapter.py"
    "meg/code/preprocessing_basic/p24_auxiliary_reconstruction.py"
    "meg/code/preprocessing_basic/preprocessing_stage_utils.py"
    "meg/code/preprocessing_basic/continuous_cleaning.py"
    "meg/code/preprocessing_basic/epoch_input_preparation.py"
    "meg/code/preprocessing_basic/epoch_annotation_trims.py"
    "meg/code/preprocessing_basic/manual_eye_annotations.py"
    "meg/code/preprocessing_basic/blink_annotation.py"
    "meg/code/preprocessing_basic/annotation_triage.py"
    "meg/code/preprocessing_basic/stage4b_adoption.py"
    "meg/code/preprocessing_basic/final_ica_review_book.py"
    "meg/code/preprocessing_basic/ica_review_plots.py"
    "meg/code/preprocessing_basic/final_ica_native_diagnostics.py"
    "meg/code/preprocessing_basic/reconcile_ica_convergence.py"
    "meg/code/preprocessing_basic/publish_ica_convergence_candidate.py"
    "meg/code/preprocessing_basic/prepare_ica_restart_candidate.py"
    "meg/code/preprocessing_basic/audit_final_ica_readiness.py"
    "meg/code/preprocessing_basic/record_provisional_ica_choices.py"
    "meg/code/preprocessing_basic/audit_ica_application.py"
    "meg/code/preprocessing_basic/verified_read_cache.py"
    "meg/code/preprocessing_basic/ica_optimizer_checkpoint.py"
    "meg/code/preprocessing_basic/diagnose_ica_convergence.py"
    "meg/code/preprocessing_basic/probe_ica_curvature_restart.py"
    "meg/code/preprocessing_basic/audit_ica_task_diagnostics.py"
    "meg/code/preprocessing_basic/review_final_ica.py"
    "meg/code/preprocessing_basic/participant_parallel.py"
    "meg/code/preprocessing_basic/review_ocular_components.py"
    "meg/code/preprocessing_basic/build_ocular_component_reports.py"
    "meg/code/preprocessing_basic/import_ocular_component_choices.py"
    "meg/code/preprocessing_basic/render_annotation_triage_report.py"
    "meg/code/preprocessing_basic/ocular_component_report.js"
    "meg/code/preprocessing_basic/ocular_component_report.html"
    "meg/code/preprocessing_basic/edf_blink_candidates.py"
    "meg/code/preprocessing_basic/review_cleaning.py"
    "meg/code/preprocessing_basic/render_stage03_mne_reports.py"
    "meg/code/preprocessing_basic/maxwell_qc.py"
    "meg/code/preprocessing_basic/squid_jumps.py"
    "meg/code/preprocessing_basic/review_squid_jumps.py"
    "meg/code/preprocessing_basic/run_pipeline.py"
    "meg/code/preprocessing_basic/trigger_extraction.py"
    "meg/code/sanity_checks/"
    "meg/params/pipeline_default.yaml"
    "meg/params/pipeline_cluster.yaml"
    "meg/params/stage3_visual_timing.yaml"
    "meg/params/sub-24-edf-aux.yaml"
    "meg/params/epoch_recovery_20261117/"
    "meg/cluster/README.md"
    "meg/cluster/environment.yml"
    "meg/cluster/requirements.txt"
    "meg/cluster/submit_preprocessing_review.sh"
    "meg/cluster/submit_ocular_component_reports.sh"
    "meg/cluster/submit_annotation_triage.sh"
    "meg/cluster/submit_stage4b_adoption.sh"
    "meg/cluster/submit_final_ica_reports.sh"
    "meg/cluster/submit_annotation_retention_report.sh"
    "meg/cluster/submit_preprocessing_stage1.sh"
    "meg/cluster/submit_preprocessing_stage2.sh"
    "meg/cluster/submit_preprocessing_stage2_summary.sh"
    "meg/cluster/submit_preprocessing_stage3.sh"
    "meg/cluster/submit_preprocessing_stage3_summary.sh"
    "meg/cluster/submit_preprocessing_apply.sh"
    "meg/cluster/submit_ica_application_qc.sh"
    "meg/cluster/submit_epoch_export.sh"
    "meg/cluster/submit_drsa_analysis.sh"
    "meg/cluster/submit_drsa_group_average.sh"
    "meg/cluster/submit_drsa_line_summaries.sh"
    "meg/cluster/submit_drsa_stats.sh"
    "meg/cluster/submit_dot_side_decoding.sh"
    "meg/cluster/submit_dot_side_group_report.sh"
    "meg/cluster/submit_dot_side_time_generalization.sh"
    "meg/cluster/submit_dot_side_time_generalization_group_report.sh"
    "meg/cluster/submit_dot_erf_analysis.sh"
    "meg/cluster/deploy_review_to_cluster.sh"
    "meg/docs/DRSA_SLIDING_LAG_SURFACES.md"
    "meg/docs/COHORT_PREPROCESSING_REVIEW_GUIDE.md"
    "meg/docs/STAGE4_CURRENT_MASK_HANDOFF.md"
    "meg/docs/FINAL_ICA_COMPONENT_REVIEW_GUIDE.md"
    "meg/docs/ICA_APPLICATION_HANDOFF_20261112.md"
    "meg/docs/INLINE_EPOCHING_HANDOFF_20261112.md"
    "meg/docs/SHORT_FRAGMENT_RETRY_20261117.md"
    "meg/docs/P19_ICA_NUMERICAL_RECONCILIATION.md"
    "meg/docs/PIPELINE_SPEC_basic.md"
    "meg/docs/meg_preprocessing_pipeline.mmd"
    "meg/docs/STAGE3_VISUAL_TIMING.md"
    "meg/docs/STAGE3_REVIEW_APPLICATION.md"
    "meg/tests/test_drsa_sliding_lag_surfaces.py"
    "meg/tests/test_drsa_stats.py"
    "meg/meg_calibration_files/fwmegmaxfilter/ct_sparse_triux2.fif"
    "meg/meg_calibration_files/fwmegmaxfilter/sss_cal_3182_250410.dat"
    "simulations/functions/"
)

# Subject selectors are explicit files in the deployment manifest so every
# array task resolves its own profile without copying anything from raw data.
for subject_number in $(seq -w 1 39); do
    FILES+=("meg/params/sub-${subject_number}.yaml")
done

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
