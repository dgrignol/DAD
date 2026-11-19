#!/usr/bin/env python3
"""
02_trigger_extraction_custom.py

Purpose
Run the configuration-driven trigger extraction stage for all 39 DAD
participants. The implementation lives in ``trigger_extraction.py`` so parser
and alignment behavior can be tested without invoking the command-line entry
point.
The current analysis cohort is P04-P39; older profiles remain for provenance.
P37's reviewed recording/block boundaries contain its unavailable first-run
tail, and the full deep report is enabled for the cohort rerun.

Observed STI101 pulses are the scientific event source. Saved
``debug_actual_triggers`` CSV files are parsed as behavioral block/run
segments and used only to resolve transition pulses and attach alignment/QC
metadata. Missing reported rows never create MEG events.
Exact human-reviewed transition overrides can preserve both recorded pulses
when a pair is confirmed to represent two consecutive logged triggers.

Usage examples
1) Extract triggers for a v18 participant:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --print-summary

2) Extract triggers for a v20 participant:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-04.yaml \
     --print-summary

3) Extract triggers for a v25 participant after the cohort preflight:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-24.yaml \
     --print-summary

4) Run without a pre-existing source manifest:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-02.yaml \
     --source-manifest ""

Outputs
- ``02_trigger_events.csv/json``: observed events with raw sample frame,
  source FIF/acquisition run, aligned behavioral block/run/trial, replay
  status, gaze-break/ESC flags, and transition provenance.
- ``02_trigger_alignment.csv``: reported-vs-observed QC, including
  reported-only and observed-only rows.
- ``02_transition_warnings.csv``: collapsed non-zero transition pulses.
- ``02_trigger_source_associations.json``: evidence-based acquisition-run to
  behavioral-segment associations.
- JSON/Markdown companion reports.
- ``02_triggers_extraction_custum_report_deeper/``: legacy-format full
  alignment CSV, anomaly CSV, block/run Markdown report, anomaly plots, and
  transition-pulse plots generated from the canonical step-02 alignment. Its
  non-MATCH tables retain raw statuses and include a separate adjudication for
  next-frame matched gaze breaks that exclude their trial by policy.

Safety
Every output path is checked by
``pipeline_config.assert_output_path_outside_raw`` before writing. Raw FIF and
CSV inputs are opened read-only.

Configuration
``--site-config`` overrides machine-specific project paths between the default
scientific YAML and the subject selector.
"""

from trigger_extraction import *  # noqa: F403
from trigger_extraction import main


if __name__ == "__main__":
    raise SystemExit(main())
