#!/usr/bin/env bash
#
# render_meg_preprocessing_pipeline.sh
#
# Purpose
# Render the editable Mermaid source for the MEG preprocessing pipeline into
# vector SVG and one-page PDF outputs using Mermaid CLI.
#
# Usage examples
# 1) From the project root:
#    ./docs/render_meg_preprocessing_pipeline.sh
#
# 2) From any directory:
#    /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/docs/render_meg_preprocessing_pipeline.sh
#
# Inputs
# - docs/meg_preprocessing_pipeline.mmd
#
# Outputs
# - docs/meg_preprocessing_pipeline.svg
# - docs/meg_preprocessing_pipeline.pdf
#
# Requirements
# - Node.js and npx must be available.
# - The script downloads/uses @mermaid-js/mermaid-cli through npx.

set -euo pipefail

# Resolve the project root from this script location so the command works
# regardless of the caller's current working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT_MMD="${PROJECT_ROOT}/docs/meg_preprocessing_pipeline.mmd"
OUTPUT_SVG="${PROJECT_ROOT}/docs/meg_preprocessing_pipeline.svg"
OUTPUT_PDF="${PROJECT_ROOT}/docs/meg_preprocessing_pipeline.pdf"

if ! command -v npx >/dev/null 2>&1; then
  echo "ERROR: npx was not found. Install Node.js/npm before rendering." >&2
  exit 1
fi

if [[ ! -f "${INPUT_MMD}" ]]; then
  echo "ERROR: Mermaid source not found: ${INPUT_MMD}" >&2
  exit 1
fi

cd "${PROJECT_ROOT}"

echo "Rendering SVG: ${OUTPUT_SVG}"
npx -y @mermaid-js/mermaid-cli \
  -i "${INPUT_MMD}" \
  -o "${OUTPUT_SVG}" \
  -b white

echo "Rendering PDF: ${OUTPUT_PDF}"
npx -y @mermaid-js/mermaid-cli \
  -i "${INPUT_MMD}" \
  -o "${OUTPUT_PDF}" \
  -b white \
  --pdfFit

echo "Rendered pipeline chart:"
echo "  ${OUTPUT_SVG}"
echo "  ${OUTPUT_PDF}"
