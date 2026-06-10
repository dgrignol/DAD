# RUN_LOG.md

## Purpose
Track concrete pipeline/test executions with commands, outcomes, and notes.

## 2026-06-10 - Configuration/Inventory Infrastructure Refactor
- Replaced numbered infrastructure scripts `00_config.py` and
  `01_inventory_and_load.py` with:
  - `pipeline_config.py`, an importable configuration and path-safety library;
  - `dataset_inventory.py`, an optional reusable preflight;
  - `run_pipeline.py`, the orchestration entrypoint.
- `project.derivatives_root` is now authoritative for derivative templates.
- Trigger extraction imports `pipeline_config.py` directly and checks every
  output against the read-only raw-data root.
- Inventory output is QC/provenance only; no preprocessing stage consumes it.
- Verification:
  - `python -m compileall -q meg/code/preprocessing_basic meg/tests`: passed.
  - `pytest -q meg/tests`: `6 passed`.
  - `run_pipeline.py --dry-run`: passed and produced the expected resolved
    config, optional inventory, and run-log files in a temporary directory.

## 2026-06-10 - Canonical Derivative Root Change
- Active configured derivative root:
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives`
- Former derivative root:
  `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives`
- The new derivative root is a sibling of the read-only
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/data` tree and passes the shared
  raw-tree output guard.

## 2026-04-21 - Full Test Suite (Current Baseline)
- Scope: `meg/tests` full suite after adding trigger parser and no-EDF behavior tests.
- Command:
  `python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests -q`
- Result:
  `3 passed in 0.72s`
- Covered checks:
  - expected-trigger-first symmetric transition handling (`A->B` and `B->A`);
  - transition fallback behavior (`width_fallback` and unresolved deterministic fallback);
  - no-EDF contract (`eye_model` skip metadata emitted, pipeline validity retained).
- Notes:
  - This run confirms current test baseline is green in the local environment.

## 2026-04-21 - Step 02 End-to-End Trigger Extraction + Deep Report
- Scope: `02_trigger_extraction_custom.py` full run on `sub-01/ses-01` with integrated deep confrontation report enabled.
- Command:
  `python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml --print-summary`
- Generated UTC (from report): `2026-04-21T13:12:04.878997Z`
- Result status:
  - step status: `warning`
  - resolved events: `3624`
  - transition warnings: `13`
  - required trigger coverage passed: `False`
  - missing required codes: `[113, 152]`
  - deep trigger report status: `ok`
- Key outputs:
  - events CSV:
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/preprocessing/sub-01/ses-01/events/02_trigger_events.csv`
  - step companion JSON:
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/reports/sub-01/ses-01/02_trigger_extraction_custom_report.json`
  - step companion Markdown:
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/reports/sub-01/ses-01/02_trigger_extraction_custom_report.md`
  - deep confrontation report Markdown:
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/reports/sub-01/ses-01/02_triggers_extraction_custum_report_deeper/sub01_all_blocks_reported_vs_observed.md`
  - deep confrontation alignment CSV:
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/reports/sub-01/ses-01/02_triggers_extraction_custum_report_deeper/sub01_all_blocks_reported_vs_observed.csv`
