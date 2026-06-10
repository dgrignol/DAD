# MEG Preprocessing Basic Pipeline Specification

## Scope
This document tracks the implemented state of the local FLUX-aligned pipeline in:

- Project root: `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG`
- Code root: `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/code/preprocessing_basic`
- Current dataset focus: the existing scripts remain primarily validated on
  `sub-01`, `ses-01`; multisubject refactoring is still in progress.
- Current configured raw-data root: `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/data`

This is a living spec for what is implemented now, not a future plan.

## FLUX Reference Standard
Primary references used for adaptation:

- FLUX overview: <https://www.neuosc.com/flux>
- FLUX repository: <https://github.com/Neuronal-Oscillations/FLUX>
- FLUX MNE-Python notebooks: <https://github.com/Neuronal-Oscillations/FLUX/tree/main/MNEPython>

## Implemented Pipeline So Far

### Script inventory and status
| Order | Local script | Status | FLUX anchor | Current role |
| --- | --- | --- | --- | --- |
| Infra | `pipeline_config.py` | Implemented | Local infrastructure | Importable config merge/validation, authoritative derivative paths, and raw-tree output guards. |
| Infra | `dataset_inventory.py` | Implemented | FLUX first-look spirit | Reusable optional preflight; inventories files and FIF headers without producing a required pipeline derivative. |
| Runner | `run_pipeline.py` | Partially implemented | Local infrastructure | Resolves config, writes one config snapshot, runs preflight, and currently orchestrates trigger extraction. Later stages are added only after their CLI refactors. |
| 02 | `02_trigger_extraction_custom.py` | Implemented | Local trigger logic | Extract events from trigger channel and write event derivatives + reports. |
| 03_A | `03_A_MaxFilt.py` | Implemented | `ApplyingMaxFilter.ipynb` | First look, Maxwell bad-channel detection, SSS/Maxwell filter on all runs. |
| 03_B | `03_B_ArtefactAnnotation.py` | Implemented | `ArtefactAnnotation.ipynb` | Blink/muscle artifact annotations, save `*_ann.fif` and `*_ann.csv`. |
| 03_C (manual) | `03_C_ICA.py` | Implemented | `ICA.ipynb` | Manual ICA branch (FastICA + subject-specific manual exclude). |
| 03_C (auto) | `03_C_ICA_auto.py` | Implemented | `ICA.ipynb` + MEGNet extension | Automatic ICA branch using `mne-icalabel` MEGNet labels. |

### Infrastructure flow
1. `pipeline_config.py` merges and validates default/subject YAML in memory.
2. `dataset_inventory.py` performs reusable read-only discovery and header checks.
3. `run_pipeline.py` writes one resolved-config snapshot, runs preflight, and
   invokes only stages with a stable noninteractive CLI.
4. Inventory JSON/Markdown is optional QC output. No scientific stage consumes
   it as an input.

### Data flow (current preprocessing scripts)
1. Raw runs (`SUB01_block*_meg.fif`) -> `03_A_MaxFilt.py`
2. Maxwell-filtered runs (`*_raw_sss.fif`) -> `03_B_ArtefactAnnotation.py`
3. Annotated runs (`*_ann.fif`) -> ICA branch.
4. Manual branch outputs (`*_ica.fif`) from `03_C_ICA.py`.
5. Auto branch outputs (`*_ica_auto.fif`) from `03_C_ICA_auto.py`.

### Pipeline flow chart
```mermaid
flowchart TB
    config[/"Input:<br/>default + subject YAML"/]
    config --> runner["Pipeline runner<br/>run_pipeline.py"]
    config --> config_lib["Configuration library<br/>pipeline_config.py"]
    config_lib --> runner
    runner --> inventory["Read-only preflight<br/>dataset_inventory.py"]

    raw[/"Input:<br/>raw *_meg.fif"/]
    raw --> inventory
    inventory -.-> inventory_report>"Optional QC:<br/>dataset_inventory.json/md"
    runner --> trigger_extract

    raw --> trigger_extract["Trigger extraction<br/>02_trigger_extraction_custom.py<br/><br/>Key steps/params:<br/>&bull; read STI101 trigger channel<br/>&bull; min_duration=0.0 s<br/>&bull; transition_max=1 sample<br/>&bull; expected_lookahead=6<br/>&bull; identity codes=1-20<br/>&bull; condition codes=31/41/51"]
    expected_triggers[/"Input:<br/>reported/expected trigger CSVs<br/>+ config YAML"/] -.-> trigger_extract
    trigger_extract --> trigger_csv[/"Output:<br/>02_trigger_events.csv"/]
    trigger_extract --> trigger_json[/"Output:<br/>02_trigger_events.json"/]
    trigger_extract --> transition_warnings[/"Output:<br/>02_transition_warnings.csv"/]
    trigger_extract -.-> trigger_report>"Reports:<br/>02_trigger_extraction_custom_report.md/json<br/>+ deep trigger QA bundle"]

    raw --> maxfilt["Maxwell filtering / SSS<br/>03_A_MaxFilt.py<br/><br/>Key steps/params:<br/>&bull; find_bad_channels_maxwell<br/>&bull; cross-talk/calibration if present<br/>&bull; reuse block01 bads across runs<br/>&bull; fix_mag_coil_types<br/>&bull; maxwell_filter<br/>&bull; reattach bads into info['bads']<br/>&bull; PSD QC: fmax=60, n_fft=1000"]
    maxfilt --> sss[/"Output:<br/>*_raw_sss.fif<br/>"/]
    maxfilt -.-> max_report>"Report:<br/>03_A_MaxFilt_report.html"]

    sss --> ann_script["Artefact annotation<br/>03_B_ArtefactAnnotation.py<br/><br/>Key steps/params:<br/>&bull; sub-01: skip EOG blink here<br/>&bull; muscle detection on mag channels<br/>&bull; filter_freq=110-140 Hz<br/>&bull; threshold=10<br/>&bull; min_length_good=0.2 s<br/>&bull; combine existing + muscle annotations"]
    ann_script --> ann_fif[/"Output:<br/>*_ann.fif<br/>"/]
    ann_script --> ann_csv[/"Output:<br/>*_ann.csv"/]
    ann_script -.-> ann_report>"Report:<br/>03_B_ArtefactAnnotation_report.html"]

    ann_fif --> branch{"ICA step<br/>choose one branch"}

    branch -->|default now| manual_ica["Manual ICA branch<br/>03_C_ICA.py<br/><br/>Key steps/params:<br/>&bull; preserve input info['bads']<br/>&bull; concatenate annotated runs for ICA fit copy (on_mismatch=warn)<br/>&bull; fit copy: MEG only<br/>&bull; resample=200 Hz<br/>&bull; band-pass=1-40 Hz<br/>&bull; FastICA, n_components=30<br/>&bull; random_state=97<br/>&bull; manual exclude sub-01=[1,5]<br/>&bull; EOG diagnostics threshold=0.6"]
    manual_ica --> manual_prefit[("Output:<br/>ICA model file<br/>*_ica_prefit.fif")]
    manual_ica --> manual_postreject[("Output:<br/>ICA model file<br/>*_ica_postreject.fif")]
    manual_ica --> manual_out[/"Output:<br/>*_ica.fif<br/>with propagated info['bads']"/]
    manual_ica -.-> manual_report>"Report:<br/>03_C_ICA_report.html"]

    branch -->|implemented - not default - to be improved| auto_ica["Auto ICA branch<br/>03_C_ICA_auto.py<br/><br/>Key steps/params:<br/>&bull; preserve input info['bads']<br/>&bull; concatenate annotated runs for ICA fit copy (on_mismatch=warn)<br/>&bull; fit copy: MEG only<br/>&bull; band-pass=1-100 Hz<br/>&bull; notch=50/60 Hz<br/>&bull; resample=250 Hz<br/>&bull; Infomax, n_components=20<br/>&bull; random_state=97, max_iter=auto<br/>&bull; MEGNet auto-exclude labels: eye blink/eye movement/heart beat if p>=0.80"]
    auto_ica --> auto_prefit[("Output:<br/>ICA model file<br/>*_ica_auto_prefit.fif")]
    auto_ica --> auto_postreject[("Output:<br/>ICA model file<br/>*_ica_auto_postreject.fif")]
    auto_ica --> auto_labels[/"Output:<br/>MEGNet labels CSV<br/>*_ica_auto_megnet_labels.csv"/]
    auto_ica --> auto_out[/"Output:<br/>*_ica_auto.fif<br/>with propagated info['bads']"/]
    auto_ica -.-> auto_report>"Report:<br/>03_C_ICA_auto_report.html"]

    manual_out --> cleaned_choice{"Selected cleaned ICA input"}
    auto_out --> cleaned_choice

    trigger_csv --> epoching
    cleaned_choice --> epoching["Epoching + metadata propagation<br/>04_epoching.py<br/><br/>Key steps/params:<br/>&bull; current ica_variant=manual<br/>&bull; events: 31, 41, 51<br/>&bull; map trial_sequence_identity from 02 events<br/>&bull; concatenate cleaned runs + per-run events (on_mismatch=warn)<br/>&bull; propagate bads into epochs.info['bads']<br/>&bull; tmin=-1.0, tmax=2.667 s<br/>&bull; baseline=None<br/>&bull; linear detrend=1<br/>&bull; reject_by_annotation=True<br/>&bull; reject: grad=5000e-13, mag=4e-12<br/>&bull; drop gaze-break code 150 trials<br/>&bull; drop catch code 116 for condition 31<br/>&bull; blink-center rejection=OFF"]
    epoching --> epo[/"Output:<br/>sub-*_ses-*_task-Attention_run-all_epo.fif<br/>contains trial_sequence_identity metadata<br/>and propagated epochs.info['bads']"/]
    epoching --> epoch_metadata[/"Output:<br/>*_epoch_metadata.csv"/]
    epoching -.-> epoch_figs(["QC figures:<br/>04_epoching_*.png"])
    epoching -.-> epoch_report>"Report:<br/>04_epoching_report.html"]

    epo --> drsa["Export neural tensors for dRSA<br/>06_export_neural_for_drsa.py<br/><br/>Key steps/params:<br/>&bull; input = step-04 epochs FIF only<br/>&bull; requires trial_sequence_identity metadata<br/>&bull; exclude bad MEG channels from export<br/>&bull; guard channel order before MAT export<br/>&bull; identity values=1-20<br/>&bull; crop=0-2.667 s<br/>&bull; resample=120 Hz<br/>&bull; average by condition x identity"]
    drsa --> mat[/"Output:<br/>condition-specific MAT files<br/>*_condition-*_desc-identityAveraged_neural.mat<br/><br/>Variables:<br/>neural = 20 x good MEG channels x timepoints<br/>channel_labels = exported channel order"/]
    drsa -.-> drsa_figs(["QC figures:<br/>06_export_neural_for_drsa_*.png"])
    drsa -.-> drsa_report>"Report:<br/>06_export_neural_for_drsa_report.html"]

    classDef inactive fill:#f6f6f6,stroke:#9a9a9a,color:#666,stroke-dasharray: 4 2;
    class auto_ica,auto_prefit,auto_postreject,auto_labels,auto_out,auto_report inactive;

```

## Step Details (As Implemented)

### 03_A MaxFilt
- FLUX-aligned operations:
- `find_bad_channels_maxwell` on block01 with cross-talk/calibration sidecars when present.
- First-look QC raw + PSD figures and report logging.
- `fix_mag_coil_types()` before Maxwell filtering.
- `maxwell_filter()` for each run.
- Local adaptation currently in script:
- Detect bad channels for each additional run as a consistency check, but still reuse block01 bads for SSS application (FLUX-style reuse policy).
- Inputs:
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/data/sub01/SUB01_block*_meg.fif`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/data/fwmegmaxfilter/ct_sparse_triux2.fif`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/data/fwmegmaxfilter/sss_cal_3182_250410.dat`
- Outputs:
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/sss/*_raw_sss.fif`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/reports/sub-01/ses-01/03_A_MaxFilt_report.html`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/figures/03_A_MaxFilt_*.png`

### 03_B ArtefactAnnotation
- FLUX-aligned operations:
- Load `*_raw_sss.fif`.
- Blink detection from `EOG001` (FLUX strategy).
- Muscle annotation via `annotate_muscle_zscore` with magnetometers (`110-140 Hz`, threshold `10`).
- Merge annotations and save annotated files.
- Save QC figures and report.
- Local adaptation currently in script:
- Subject-specific branch for `sub-01`: EOG blink annotation is skipped in 03_B (deferred downstream).
- Inputs:
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/sss/*_raw_sss.fif`
- Outputs:
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/sss/*_ann.fif`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/sss/*_ann.csv`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/reports/sub-01/ses-01/03_B_ArtefactAnnotation_report.html`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/figures/03_B_ArtefactAnnotation_*.png`

### 03_C ICA (manual branch)
- FLUX-aligned operations:
- Concatenate all annotated runs for ICA fit copy.
- MEG-only fit copy, resample `200 Hz`, band-pass `1-40 Hz`.
- FastICA fit (`n_components=30`, `random_state=97`).
- Inspect sources/components and apply ICA to original annotated raws.
- Local adaptations currently in script:
- `on_mismatch="warn"` during run concatenation because `dev_head_t` differs across runs.
- ICA cache split into prefit and postreject model files.
- Manual exclude list for `sub-01`.
- Subject-specific blink workaround in this branch uses ICA000-derived blink annotations.
- Inputs:
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/sss/*_ann.fif`
- Outputs:
- Cleaned runs: `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/sss/*_ica.fif`
- ICA models:
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/ica_models/sub-01_ses-01_ica_prefit.fif`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/ica_models/sub-01_ses-01_ica_postreject.fif`
- Report: `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/reports/sub-01/ses-01/03_C_ICA_report.html`

### 03_C ICA auto (MEGNet branch)
- Purpose:
- Alternative to manual component rejection; keeps same input/output stage as 03_C but auto-labels components with MEGNet.
- Fit-copy settings:
- Concatenate all annotated runs for ICA fit copy (`on_mismatch="warn"` for `dev_head_t` differences).
- MEG-only, band-pass `1-100 Hz`, notch `50/60 Hz`, resample `250 Hz`.
- ICA method `infomax`, current default `n_components=20`.
- Auto exclusion policy:
- Target labels: `eye blink`, `eye movement`, `heart beat`.
- Exclude when predicted label probability >= `0.80`.
- `brain/other` is not auto-excluded by this rule.
- Current subject-specific behavior for `sub-01`:
- Select MEGNet-labeled blink component and inject blink windows before saving cleaned runs.
- Inputs:
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/sss/*_ann.fif`
- Outputs:
- Cleaned runs: `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/sss/*_ica_auto.fif`
- ICA models:
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/ica_models/sub-01_ses-01_ica_auto_prefit.fif`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/ica_models/sub-01_ses-01_ica_auto_postreject.fif`
- Labels CSV:
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/preprocessing/sub-01/ses-01/ica_models/sub-01_ses-01_ica_auto_megnet_labels.csv`
- Report: `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/reports/sub-01/ses-01/03_C_ICA_auto_report.html`

## Current Derivative Footprint (Observed)
Current `sub-01/ses-01` derivatives include:

- SSS files: `*_raw_sss.fif`
- Annotation files: `*_ann.fif`, `*_ann.csv`
- Manual ICA cleaned files: `*_ica.fif`
- Auto ICA cleaned files: `*_ica_auto.fif`
- ICA model caches in `ica_models/` for both branches
- HTML reports:
- `03_A_MaxFilt_report.html`
- `03_B_ArtefactAnnotation_report.html`
- `03_C_ICA_report.html`
- `03_C_ICA_auto_report.html`

## Logged FLUX Deviations
Deviation notes currently tracked in:

- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/docs/deviations/DEVIATION_03B_SUB01_SKIP_EOG_BLINK.md`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/docs/deviations/DEVIATION_03C_ICA_DEV_HEAD_T.md`

Required rule for future updates:

- Any new FLUX deviation must be documented in `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/docs/deviations/` with the same template fields already used in existing deviation notes.

## Next Planned Stages (Not Yet Implemented in `preprocessing_basic`)
- Epoching / condition extraction integration after ICA branch decision.
- Downstream export for analysis-specific models.

## Keep-in-Sync Checklist
When behavior changes, update this file sections:

- Script inventory and status
- Infrastructure flow and runner-supported stages
- Step details and parameter defaults
- Derivative filenames and paths
- Deviation log index
