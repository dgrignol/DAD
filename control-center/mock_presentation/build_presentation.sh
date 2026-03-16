#!/usr/bin/env bash
# build_presentation.sh
#
# Purpose:
#   Build the complete MoveDot1 MEG meeting deck as a PowerPoint file.
#   The script first regenerates figure assets (real + synthetic), then
#   renders the markdown slides into a .pptx using pandoc.
#
# Usage examples:
#   bash control-center/presentation/build_presentation.sh
#   bash control-center/presentation/build_presentation.sh --subject 80
#
# Inputs:
#   - meg_meeting_slides.md
#   - generate_expected_figures.py
#   - required MAT files loaded by generate_expected_figures.py
#
# Outputs:
#   - control-center/presentation/figures/*.png
#   - control-center/presentation/MoveDot1_MEG_meeting_SubXX.pptx
#
# Key assumptions:
#   - python, pandoc, numpy, scipy, matplotlib are installed.
#   - script is run from any location (it resolves its own directory).

set -euo pipefail

# ---------------------------------------------------------------------------
# Parse CLI arguments
# Data flow: optional subject arg -> output filename + figure generation call.
# ---------------------------------------------------------------------------
SUBJECT="80"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --subject)
            SUBJECT="$2"
            shift 2
            ;;
        -h|--help)
            cat <<'EOF'
Usage: build_presentation.sh [--subject XX]

Build the MoveDot1 MEG presentation by regenerating figures and rendering
meg_meeting_slides.md to PowerPoint.
EOF
            exit 0
            ;;
        *)
            echo "[error] Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Resolve paths and tool binaries
# Data flow: script dir -> source files -> deterministic output locations.
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
PANDOC_BIN="${PANDOC_BIN:-pandoc}"

FIGURE_SCRIPT="$SCRIPT_DIR/generate_expected_figures.py"
SLIDES_MD="$SCRIPT_DIR/meg_meeting_slides.md"
OUTPUT_PPTX="$SCRIPT_DIR/MoveDot1_MEG_meeting_Sub$(printf "%02d" "$SUBJECT").pptx"

if [[ ! -f "$FIGURE_SCRIPT" ]]; then
    echo "[error] Missing figure script: $FIGURE_SCRIPT" >&2
    exit 1
fi
if [[ ! -f "$SLIDES_MD" ]]; then
    echo "[error] Missing slides markdown: $SLIDES_MD" >&2
    exit 1
fi
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    echo "[error] Python binary not found: $PYTHON_BIN" >&2
    exit 1
fi
if ! command -v "$PANDOC_BIN" >/dev/null 2>&1; then
    echo "[error] Pandoc binary not found: $PANDOC_BIN" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Generate figure assets
# Data flow: subject + outdir -> figures/*.png referenced by markdown slides.
# ---------------------------------------------------------------------------
echo "[build] Generating figure assets for subject ${SUBJECT}..."
"$PYTHON_BIN" "$FIGURE_SCRIPT" --subject "$SUBJECT" --outdir "$SCRIPT_DIR"

# ---------------------------------------------------------------------------
# Render markdown slides to PowerPoint
# Data flow: markdown + local figure references -> output .pptx file.
# ---------------------------------------------------------------------------
echo "[build] Rendering PowerPoint deck..."
(
    cd "$SCRIPT_DIR"
    "$PANDOC_BIN" "$SLIDES_MD" -t pptx -o "$OUTPUT_PPTX"
)

if [[ ! -s "$OUTPUT_PPTX" ]]; then
    echo "[error] PowerPoint build failed or produced empty file: $OUTPUT_PPTX" >&2
    exit 1
fi

echo "[ok] Presentation built: $OUTPUT_PPTX"
