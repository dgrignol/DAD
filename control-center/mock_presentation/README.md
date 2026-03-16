# MoveDot1 MEG meeting presentation

This folder contains a reproducible pipeline to generate figure assets and build
an editable PowerPoint presentation for the MoveDot1 MEG meeting.

## Files

- `meg_meeting_slides.md`: 12-slide markdown source for the talk.
- `generate_expected_figures.py`: figure generator (real Sub80 + synthetic expected patterns).
- `build_presentation.sh`: one-command build script.
- `figures/`: generated PNG assets referenced by slides.
- `MoveDot1_MEG_meeting_Sub80.pptx`: final deck output for subject 80.

## Requirements

- `python`
- `pandoc`
- Python packages: `numpy`, `scipy`, `matplotlib`

## Usage

Build default deck (subject 80):

```bash
bash control-center/presentation/build_presentation.sh
```

Build for a different subject id (if corresponding files exist):

```bash
bash control-center/presentation/build_presentation.sh --subject 80
```

Generate figures only:

```bash
python control-center/presentation/generate_expected_figures.py --subject 80 --outdir control-center/presentation
```

## Expected outputs

Figure files in `control-center/presentation/figures`:

- `fig01_paths_sub80_realistic.png`
- `fig02_deviant_predicted_divergence.png`
- `fig03_mock_predrift_lag.png`
- `fig04_mock_postdeviant_switch.png`
- `fig05_mock_drsa_matrix.png`

PowerPoint deck:

- `control-center/presentation/MoveDot1_MEG_meeting_Sub80.pptx`

## Acceptance checks

Run these checks after building:

```bash
# 1) Build success and non-empty pptx
bash control-center/presentation/build_presentation.sh
ls -lh control-center/presentation/MoveDot1_MEG_meeting_Sub80.pptx

# 2) Asset completeness
ls control-center/presentation/figures/fig0*.png

# 3) Markdown references resolve
rg "figures/fig" control-center/presentation/meg_meeting_slides.md
```

## Notes on figure interpretation

- `fig01` and `fig02` are derived from real Sub80 trajectory files.
- `fig03`, `fig04`, and `fig05` are deterministic synthetic dRSA illustrations
  and are labeled on-figure as synthetic.
- In synthetic dRSA figures, axis semantics are explicit:
  - x-axis = neural time
  - y-axis = model time
