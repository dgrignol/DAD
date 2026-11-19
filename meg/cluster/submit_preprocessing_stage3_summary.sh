#!/bin/bash
# Purpose: aggregate 36 intact Stage 3 proposals without accepting their science.
# Usage after a successful participant array, from the deployed DAD root:
#   DAD_STAGE3_VERSION=s3-v001 sbatch --dependency="afterok:${stage3_job}" \
#     meg/cluster/submit_preprocessing_stage3_summary.sh
# Inputs: full P04-P39 proposals, current accepted Stage 2 evidence and method.
# Outputs: immutable derivatives/visual_timing/VERSION/group cohort report.
# Data flow: single job -> shared Stage 3 environment wrapper -> aggregation.
# Pilot/missing/stale/edited bundles stop this job. No signals/epochs are read.
# This wrapper never submits another job; the user controls execution.

#SBATCH --job-name=DAD_meg_stage3_summary
#SBATCH --output=cluster_logs/%x.%j.log
#SBATCH --error=cluster_logs/%x.%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --export=ALL

set -euo pipefail
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
export DAD_STAGE3_MODE=summary
exec bash "$PROJECT_ROOT/meg/cluster/submit_preprocessing_stage3.sh"
