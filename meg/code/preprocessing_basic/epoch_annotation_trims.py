#!/usr/bin/env python3
"""Apply exact user-approved annotation end trims in a new epoch variant.

Usage through Stage 04::

    DAD_EPOCH_PREPARATION=inline_no_ica \
    DAD_PROCESSING_RUN=noica-recovery-20261117-v1 python 04_epoching.py \
      --config ../../params/pipeline_default.yaml \
      --subject-config ../../params/epoch_recovery_20261117/sub-21.yaml \
      --tuning-config /derivatives/reports/sub-21/ses-01/03_tuning_choices.yaml \
      --source-manifest /derivatives/reports/sub-21/ses-01/source_manifest.json \
      --timing-review /derivatives/visual_timing/s3-v002/reviews/sub-21_review.json

Input rules identify the exact source, reviewed trial key/onset, original whole
annotation hash and one exact annotation row. The approved end is one native
sample before the full epoch starts. Only that annotation duration is shortened;
its onset, label, channel scope and every other annotation stay unchanged.

Validation binds the requested trials to the accepted Stage-2/3 metadata before
processing. Operative trims are applied in memory before continuous line removal
and epoching, so both operators see the requested annotation. Source FIFs,
original mask files and ICA models remain unchanged. The changed setting and
implementation bind a new processing-run contract; earlier epochs stay intact.
This narrowly scoped adapter supports only explicitly labelled P21 no-ICA runs.
It neither disables amplitude rejection nor introduces general trial overrides.
"""

from copy import deepcopy

import numpy as np

from continuous_cleaning import annotation_rows, digest


# -----------------------------------------------------------------------------
# Exact annotation identity and native sample arithmetic; no fuzzy trial repair.
# -----------------------------------------------------------------------------

def plan_trims(raw, cfg, source_filename):
    """Validate exact native annotation rows and return bound shortening plans."""
    rules=[r for r in cfg["epoching"].get("annotation_end_trims",[]) if r["source_fif"]==source_filename]
    if not rules:return []
    if (str(cfg["subject"]["id"]).zfill(2)!="21" or cfg["epoching"].get("preparation_mode")!="inline_no_ica"
            or not cfg["epoching"].get("processing_run")):
        raise RuntimeError("Annotation end trims require a new labelled P21 no-ICA variant")
    original=digest(annotation_rows(raw));plans=[];used=set();sfreq=float(raw.info["sfreq"])
    for rule in rules:
        if original!=rule["original_annotations_sha256"]:
            raise RuntimeError("Source annotations changed since the approved P21 review")
        anchor=rule["epoch_event_sample_raw"]
        start=anchor+round(cfg["epoching"]["tmin_sec"]*sfreq)
        if rule["new_end_sample_raw"]!=start-1:
            raise RuntimeError("Approved annotation must end one native sample before the full epoch")
        indices=[i for i,a in enumerate(raw.annotations)
            if str(a["description"])==rule["annotation_description"]
            and np.isclose(a["onset"],rule["annotation_onset_sec"],rtol=0,atol=1e-9)
            and np.isclose(a["duration"],rule["annotation_duration_sec"],rtol=0,atol=1e-9)]
        if len(indices)!=1 or indices[0] in used:
            raise RuntimeError("Approved P21 annotation row is absent, duplicated or reused")
        index=indices[0];used.add(index);ann=raw.annotations[index]
        end=rule["new_end_sample_raw"]/sfreq;old_end=float(ann["onset"]+ann["duration"])
        if not (raw.first_samp<=start<raw.first_samp+raw.n_times
                and ann["onset"]<end<start/sfreq<old_end<anchor/sfreq
                and str(ann["description"]).lower().startswith("bad")):
            raise RuntimeError("Approved trim must shorten a pre-onset BAD overlap only")
        plans.append({**deepcopy(rule),"annotation_index":index,"new_duration_sec":end-float(ann["onset"]),
            "epoch_first_sample_raw":start})
    return plans


def apply_trims(raw, cfg, source_filename):
    """Shorten only the approved durations on the disposable epoch input Raw."""
    plans=plan_trims(raw,cfg,source_filename)
    if not plans:return []
    before=annotation_rows(raw)
    for plan in plans:raw.annotations.duration[plan["annotation_index"]]=plan["new_duration_sec"]
    expected=deepcopy(before)
    for plan in plans:expected["rows"][plan["annotation_index"]][1]=plan["new_duration_sec"]
    if annotation_rows(raw)!=expected:raise RuntimeError("An unapproved annotation property changed")
    return plans


# -----------------------------------------------------------------------------
# Bind requested trial identities to accepted eligibility and visual timing.
# -----------------------------------------------------------------------------

def validate_trial_bindings(cfg, trials):
    """Reject stale/mistyped source, behavioral trial or corrected-onset keys."""
    rules=cfg["epoching"].get("annotation_end_trims",[])
    if not rules:return
    if trials is None:raise RuntimeError("P21 annotation trims require accepted trial/timing metadata")
    for rule in rules:
        matched=[r for r in trials if r["trial_key"]==rule["trial_key"]]
        if len(matched)!=1:raise RuntimeError("Approved P21 trial key is absent or duplicated")
        row=matched[0]
        checks={"source_fif":rule["source_fif"],"epoch_event_sample_raw":rule["epoch_event_sample_raw"],
            "behavioral_run":rule["behavioral_run"],"behavioral_trial":rule["behavioral_trial"]}
        for key,value in checks.items():
            if str(row[key]).removesuffix(".0")!=str(value):raise RuntimeError(f"Approved P21 trial binding differs: {key}")
        if (str(row.get("final_analysis_eligible")).lower()!="true" or int(float(row["trial_sequence_identity"]))!=31
                or row["condition_name"]!="occluded_deviant"):
            raise RuntimeError("Annotation trim cannot change trial eligibility, condition or identity")
