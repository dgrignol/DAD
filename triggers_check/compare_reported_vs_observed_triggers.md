# compare_reported_vs_observed_triggers.py

## What this script does
`compare_reported_vs_observed_triggers.py` compares:
1. Reported triggers from experiment debug CSV files.
2. Observed trigger windows reconstructed from a `.fif` trigger channel.

It produces:
- A full alignment CSV (`MATCH`, `MISMATCH`, `REPORTED_ONLY`, `OBSERVED_ONLY`).
- An anomaly CSV.
- A Markdown report with summary tables and embedded plots.
- Per-anomaly plots and transition-pulse warning plots.

## Does it accept explicit input file paths?
Yes.
You can pass both required inputs explicitly with:
- `--reported-csv /absolute/path/to/debug_actual_triggers....csv`
- `--fif /absolute/path/to/file.fif`

Example:
```bash
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/compare_reported_vs_observed_triggers.py \
  --reported-csv /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/experiment/output_files/debug_actual_triggers_occlusion_v18_rescueTraject_sub01_block01_run02.csv \
  --fif /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/data/Subj01_4blocks.fif
```

## Does it accept an output folder?
There is no single `--output-dir` flag, but yes, you can route outputs to one folder by setting all output-path flags:
- `--out-csv`
- `--out-report`
- `--anomaly-csv`
- `--plots-dir`

Example:
```bash
OUT_DIR=/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/reports/sub01_block01_run02_custom
mkdir -p "$OUT_DIR"

python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/compare_reported_vs_observed_triggers.py \
  --reported-csv /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/experiment/output_files/debug_actual_triggers_occlusion_v18_rescueTraject_sub01_block01_run02.csv \
  --fif /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/data/Subj01_4blocks.fif \
  --out-csv "$OUT_DIR/confrontation.csv" \
  --out-report "$OUT_DIR/confrontation.md" \
  --anomaly-csv "$OUT_DIR/confrontation_anomalies.csv" \
  --plots-dir "$OUT_DIR/plots"
```

## Auto-resolution mode
If you do not pass `--reported-csv` and/or `--fif`, the script can auto-resolve from subject/block/run:

```bash
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/compare_reported_vs_observed_triggers.py \
  --subject 1 --block 1 --run 2
```

Reported CSV auto-resolution pattern:
- `experiment/output_files/debug_actual_triggers*_subXX_blockYY_runZZ.csv`

Optional filter:
- `--variant rescueTraject` keeps only names containing `_rescueTraject_`.

If multiple CSVs match, the script picks the best candidate by:
1. Highest parsed `_vNN` version number.
2. Newest file modification time.
3. Lexical filename order.

## Input expectations
### Reported CSV required columns
- `trigger`
- `trial`
- `frame`
- `seconds`
- `label`

### FIF requirements
- Must contain the selected trigger channel (default `STI101`).

## Flag reference (all flags)
Each row includes one concrete usage example.

| Flag | Purpose | Example |
| --- | --- | --- |
| `-h, --help` | Show CLI help. | `python .../compare_reported_vs_observed_triggers.py --help` |
| `--project-root` | Root used for auto-resolving input/output defaults. | `python .../compare_reported_vs_observed_triggers.py --project-root /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check --subject 1 --block 1 --run 2` |
| `--subject` | Subject number for auto-resolution. | `python .../compare_reported_vs_observed_triggers.py --subject 1 --block 1 --run 2` |
| `--block` | Block number for auto-resolution. | `python .../compare_reported_vs_observed_triggers.py --subject 1 --block 2 --run 1` |
| `--run` | Run number for auto-resolution. | `python .../compare_reported_vs_observed_triggers.py --subject 1 --block 1 --run 3` |
| `--variant` | Optional token filter for auto-resolved CSV names. | `python .../compare_reported_vs_observed_triggers.py --subject 1 --block 1 --run 2 --variant rescueTraject` |
| `--reported-csv` | Explicit debug CSV path. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug_actual_triggers_occlusion_v18_rescueTraject_sub01_block01_run02.csv --fif /abs/path/Subj01_4blocks.fif` |
| `--fif` | Explicit FIF path. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif` |
| `--trigger-channel` | Trigger channel name in FIF. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --trigger-channel STI101` |
| `--min-duration` | Minimum duration (s) used in event merge logic. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --min-duration 0.001` |
| `--transition-max-samples` | Max width (samples) treated as transient transition pulse. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --transition-max-samples 2` |
| `--out-csv` | Output path for full alignment CSV. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --out-csv /abs/out/confrontation.csv` |
| `--out-report` | Output path for Markdown report. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --out-report /abs/out/confrontation.md` |
| `--anomaly-csv` | Output path for anomaly CSV. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --anomaly-csv /abs/out/confrontation_anomalies.csv` |
| `--plots-dir` | Folder for generated plot PNGs. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --plots-dir /abs/out/plots` |
| `--plot-window-ms` | Half-window (ms) for anomaly plots. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --plot-window-ms 750` |
| `--transition-plot-window-ms` | Half-window (ms) for transition warning plots. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --transition-plot-window-ms 3` |
| `--disable-plot-transition-pulses` | Skip generating transition warning plots (warnings still in summary/report). | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --disable-plot-transition-pulses` |
| `--include-width-outliers` | Add WIDTH_OUTLIER anomalies. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --include-width-outliers` |
| `--open-scroller` | Open interactive trigger timeline scroller after generation. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --open-scroller` |
| `--scroller-duration` | Visible duration (s) for scroller window. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --open-scroller --scroller-duration 20` |
| `--scroller-start` | Initial scroller window start time (s). | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --open-scroller --scroller-start 120.5` |
| `--scroller-step` | Keyboard navigation step (s) in scroller. | `python .../compare_reported_vs_observed_triggers.py --reported-csv /abs/path/debug.csv --fif /abs/path/Subj01_4blocks.fif --open-scroller --scroller-step 2.5` |

## Practical command templates
### Template 1: explicit files + explicit outputs (recommended for reproducibility)
```bash
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/compare_reported_vs_observed_triggers.py \
  --reported-csv /ABS/PATH/debug_actual_triggers_...csv \
  --fif /ABS/PATH/file.fif \
  --out-csv /ABS/OUT/confrontation.csv \
  --out-report /ABS/OUT/confrontation.md \
  --anomaly-csv /ABS/OUT/confrontation_anomalies.csv \
  --plots-dir /ABS/OUT/plots
```

### Template 2: auto-resolve by subject/block/run
```bash
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/compare_reported_vs_observed_triggers.py \
  --subject XX --block YY --run ZZ
```

### Template 3: auto-resolve + variant filter
```bash
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/compare_reported_vs_observed_triggers.py \
  --subject XX --block YY --run ZZ --variant rescueTraject
```
