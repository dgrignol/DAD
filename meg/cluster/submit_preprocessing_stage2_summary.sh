#!/bin/bash
#
# Purpose: summarize all 36 completed Stage 2 participant bundles.
# Usage from the deployed DAD root, after submitting the participant array:
#   sbatch --dependency="afterok:${stage2_job}" meg/cluster/submit_preprocessing_stage2_summary.sh
# Inputs: P04-P39 Stage 2 bundles and their unchanged Stage 1 evidence.
# Outputs: reports/group/trial_preparation cohort MD/JSON and count CSVs.
# Data flow: single SLURM job -> shared environment wrapper -> cohort aggregator.
# Missing, stale, or altered bundles stop aggregation. This job never submits
# another job. Explicit DAD_STAGE2_OVERWRITE=1 affects only summary outputs.

#SBATCH --job-name=DAD_meg_stage2_summary
#SBATCH --output=cluster_logs/%x.%j.log
#SBATCH --error=cluster_logs/%x.%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --export=ALL

set -euo pipefail
PROJECT_ROOT="${DAD_PROJECT_ROOT:-/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD}"
export DAD_STAGE2_MODE=summary
exec bash "$PROJECT_ROOT/meg/cluster/submit_preprocessing_stage2.sh"
