#!/usr/bin/env python3
"""Stage a hash-bound post-Stage-3 eligibility and restricted-donor application.

Usage from DAD (metadata and saved measurements only; no waveform rereads)::

    python meg/code/preprocessing_basic/visual_timing_application.py \
      --config meg/params/pipeline_default.yaml --source-version s3-v001 \
      --target-version s3-v002 --workers 4 --validate-only

    python meg/code/preprocessing_basic/visual_timing_application.py \
      --config meg/params/pipeline_default.yaml --source-version s3-v001 \
      --target-version s3-v002 --workers 4

Data flow: verified original Stage 2 ledger and s3-v001 measurements -> exact
saved MAT block/run/trial/replay/condition join -> separate final-eligibility
overlay -> eight explicit whole-trial exclusions -> P24-only donor removal and
dependent estimates/proposals -> clearly marked Task-1 staging bundles/reports.
Original eligibility, keys, sources, recorded events and measured edges remain
unchanged. Saved catch_type_code=2 qualifies only motion onsets 81/91; no code
118 inference, synthetic onset, missing trial or generic exclusion waiver is
permitted. Participant-analysis eligibility never implicitly filters donors.

The final s3-v002 directory is deliberately NOT published here. Task 2 reads
the verified staging bundles with load_task1_stage(), adds its separately
supported onset estimates, and publishes one coherent final version. Staging
mode cannot pass the existing accepted-timing gate. This script creates no
accepted review.json, changes no formal decision/working log, and launches no
cleaning, epoching, export, analysis or scheduler job.
"""

from __future__ import annotations

import argparse
import copy
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import html
import io
import json
import multiprocessing
from pathlib import Path
import re
import sys

import numpy as np
from scipy.io import loadmat

from downstream_stage_utils import csv_bool
from pipeline_config import ANALYSIS_SUBJECT_IDS,load_and_resolve_config
from preprocessing_stage_utils import guarded_path,utc_now
from trial_preparation import csv_text,digest,file_hash,load_prepared_decisions,read_csv,stage2_paths
from visual_timing_utils import donor_allowed,estimate_lag,round_half_away,verify_bundle,window_coverage


APPLICATION_SCHEMA = "1.0"
POLICY_ID = "saved_occlusion_type2_and_exact_sequence_exclusions_v1"
# These are exact human decisions, not a numerical sequence-error threshold.
# Source/block identities remain distinct, including P29/P38 short recordings.
EXCLUSIONS = (
    ("24","sub24_run05_100926.fif",582257,81,5,5,2,51,True),
    ("15","sub15_run01_030926.fif",1074390,81,1,1,3,57,True),
    ("28","sub28_run05_150926.fif",352984,91,5,5,2,13,False),
    ("29","sub29_run01_150926.fif",1500707,81,1,1,3,10,False),
    ("29","sub29_run01_150926.fif",1519036,81,1,1,3,12,False),
    ("29","sub29_run01_150926.fif",1523905,81,1,1,3,13,False),
    ("33","sub33_run04_220926.fif",811230,81,4,4,3,39,False),
    ("38","sub38_run02_051026.fif",622800,81,2,2,2,24,False),
)
P24_DONOR_KEYS = tuple(f"sub-24:sub24_run05_100926.fif:{sample}:{code}" for code,sample in
                       ((81,582257),(111,583362),(112,583496),(114,584004),(115,584138)))
IDENTITY_FIELDS = ("source_fif","raw_sample_frame","condition_code","acquisition_run",
                   "behavioral_block","behavioral_run","behavioral_trial","is_replay_trial")


# -----------------------------------------------------------------------------
# Exact saved-trial metadata and the narrow, independent eligibility revision
# -----------------------------------------------------------------------------


def policy_definition():
    """Return the complete supported policy for immutable application records."""
    return dict(policy_id=POLICY_ID,saved_catch_type_code=2,condition_codes=[81,91],
        waived_reason="catch_type_2",preserve_all_independent_exclusions=True,
        sequence_exclusions=[dict(zip(("subject","source_fif","raw_sample_frame","condition_code",
            "acquisition_run","behavioral_block","behavioral_run","behavioral_trial","is_replay_trial"),case)) for case in EXCLUSIONS],
        donor_excluded_event_keys=list(P24_DONOR_KEYS),donor_filter_is_independent_of_analysis_eligibility=True)


def context_key(row):
    """Canonical saved context; replay is an identity dimension, not row order."""
    replay=csv_bool(row["is_replay_trial"]) if "is_replay_trial" in row else row.get("replay_status")=="replay"
    values=[None if row.get(k) in (None,"") else int(row[k]) for k in ("behavioral_block","behavioral_run","behavioral_trial")]
    return (*values,replay)


def portable_basename(value):
    """Saved lab paths use Windows separators; mounted ledger paths use POSIX."""
    return Path(str(value).replace("\\","/")).name


def scalar_integer(obj,name,default=None):
    """Read only finite scalar MAT metadata; never interpret an array as a class."""
    value=getattr(obj,name,default)
    if value is None:return None
    value=np.asarray(value)
    if value.size!=1 or not np.isfinite(float(value.item())):
        raise RuntimeError(f"Saved {name} is not a finite scalar.")
    number=float(value.item())
    if number!=int(number):raise RuntimeError(f"Saved {name} is not an integer.")
    return int(number)


def read_saved_metadata(cfg,source_manifest,contexts):
    """Hash/parse only existing saved MAT runResults; return exact trial metadata.

    The recorded source manifest supplies filenames and expected hashes. MAT
    outputs without recorded onsets are retained as metadata, never new events.
    """
    folder=Path(cfg["resolved_paths"]["subject_data_dir"])/"output_files"
    result=[]
    for block,context in sorted(contexts.items(),key=lambda item:int(item[0])):
        name=context["source_mat"];path=folder/name
        expected=source_manifest["saved_mat_sha256"].get(name)
        if not expected or file_hash(path)!=expected:raise RuntimeError(f"Saved MAT evidence changed: {path}")
        saved=loadmat(path,variable_names=["runResults"],squeeze_me=True,struct_as_record=False)
        for run in np.atleast_1d(saved.get("runResults",[])):
            trials=[trial for trial in np.atleast_1d(getattr(run,"output",[])) if hasattr(trial,"trial_num")]
            # Interrupted/short blocks retain empty planned-run placeholders.
            # They are not saved trial outputs and carry no usable run number.
            if not trials:continue
            run_number=scalar_integer(run,"run_num")
            debug=str(getattr(run,"debugTriggerFile","") or "")
            for trial in trials:
                if scalar_integer(trial,"block_num",int(block))!=int(block):
                    raise RuntimeError("Saved trial block disagrees with its MAT file.")
                catch=scalar_integer(trial,"catch_type_code")
                condition=scalar_integer(trial,"condition_trigger")
                if catch not in (0,1,2) or condition not in (71,81,91):
                    raise RuntimeError(f"Unsupported/missing saved catch class or condition: {name}")
                result.append(dict(subject_bids_id=cfg["layout"]["subject_bids_id"],source_mat=name,
                    source_mat_sha256=expected,behavioral_block=int(block),behavioral_run=run_number,
                    behavioral_trial=scalar_integer(trial,"trial_num"),
                    is_replay_trial=bool(scalar_integer(trial,"is_replay_trial",0)),
                    saved_condition_code=condition,saved_catch_type_code=catch,
                    saved_sequence_identity=scalar_integer(trial,"sequence_trigger_code"),
                    saved_trial_completed=scalar_integer(trial,"trial_completed"),
                    saved_frame_missed_count=scalar_integer(trial,"frame_missed_count"),
                    debug_source_csv=portable_basename(debug) if debug else ""))
        if file_hash(path)!=expected:raise RuntimeError("Saved MAT changed during metadata parsing.")
    keys=[context_key(r) for r in result]
    if len(set(keys))!=len(keys):raise RuntimeError("Ambiguous saved block/run/trial/replay identities.")
    return result


def build_eligibility_overlay(ledger,saved_metadata,subject,*,require_exact_exclusions=True):
    """Apply only the saved-class waiver and explicit whole-trial exclusions.

    Return one overlay row for every original ledger key, preserving original
    eligibility/reasons. Catch inclusion removes only catch_type_2; gaze, ESC,
    unavailable tails, type-1 and every other reason remain independent. Exact
    sequence decisions override the waiver. No timing or donor row is altered.
    """
    saved={context_key(r):r for r in saved_metadata}
    if len(saved)!=len(saved_metadata):raise RuntimeError("Duplicate saved trial context.")
    cases={f"sub-{case[0]}:{case[1]}:{case[2]}":case for case in EXCLUSIONS if case[0]==subject}
    result=[];matched=set();seen_context=set()
    for original in ledger:
        key=context_key(original)
        metadata=saved.get(key)
        if metadata and key in seen_context:raise RuntimeError("Saved trial context joins multiple recorded onset keys.")
        if metadata:seen_context.add(key)
        # Accepted interrupted recordings can retain an onset absent from the
        # saved completed MAT output. Preserve their original disposition; an
        # absent class never authorizes a new inclusion or a fabricated trial.
        metadata=metadata or dict(source_mat="",saved_catch_type_code=None,saved_condition_code=None,debug_source_csv="")
        if metadata["saved_condition_code"] is not None and int(original["condition_code"])!=int(metadata["saved_condition_code"]):
            raise RuntimeError("Saved condition contradicts recorded trial identity.")
        debug=metadata.get("debug_source_csv","")
        if debug and debug!=portable_basename(original["reported_source_csv"]):
            raise RuntimeError("Saved debug CSV contradicts accepted source/trial context.")
        before=csv_bool(original["included_by_trial_policy"])
        reasons=[r for r in str(original.get("exclusion_reasons","")).split(",") if r]
        if before!=(not bool(reasons)):raise RuntimeError("Original eligibility/reason accounting is inconsistent.")
        qualifies=metadata["saved_catch_type_code"] is not None and int(metadata["saved_catch_type_code"])==2 and int(original["condition_code"]) in (81,91)
        remaining=[r for r in reasons if not (qualifies and r=="catch_type_2")]
        case=cases.get(original["trial_key"])
        if case:
            actual=(subject,str(original["source_fif"]),int(original["raw_sample_frame"]),int(original["condition_code"]),
                int(original["acquisition_run"]),int(original["behavioral_block"]),int(original["behavioral_run"]),
                int(original["behavioral_trial"]),csv_bool(original["is_replay_trial"]))
            if actual!=case:raise RuntimeError("Exact sequence-exclusion identity mismatch.")
            remaining.append("stage3_sequence_exclusion");matched.add(original["trial_key"])
        if metadata["saved_catch_type_code"] is not None and int(metadata["saved_catch_type_code"])==1 and not remaining:
            raise RuntimeError("A saved type-1 catch would become eligible.")
        result.append(dict(trial_key=original["trial_key"],**{f:original[f] for f in IDENTITY_FIELDS},
            source_mat=metadata["source_mat"],saved_catch_type_code=metadata["saved_catch_type_code"],
            saved_metadata_join="exact_block_run_trial_replay_condition" if metadata["source_mat"] else "unavailable_original_disposition_preserved",
            saved_condition_code=metadata["saved_condition_code"],reported_source_csv=portable_basename(original["reported_source_csv"]),
            original_included_by_trial_policy=before,original_exclusion_reasons=original.get("exclusion_reasons",""),
            original_accepted_exclusion_reason=original.get("accepted_exclusion_reason",""),
            qualifies_saved_occlusion_type2=qualifies,waived_catch_reason="catch_type_2" if qualifies and "catch_type_2" in reasons else "",
            sequence_exclusion=bool(case),sequence_decision_reference=f"review_notes.md#p{subject}" if case else "",
            final_analysis_eligible=not remaining,analysis_exclusion_reasons=",".join(remaining),
            newly_included=(not before and not remaining),newly_excluded=(before and bool(remaining))))
    if require_exact_exclusions and matched!=set(cases):raise RuntimeError("A confirmed whole-trial exclusion is missing.")
    return result


def attach_eligibility(maps,overlay):
    """Add final eligibility without changing the original included-policy field."""
    lookup={r["trial_key"]:r for r in overlay}
    if len(lookup)!=len(overlay) or set(lookup)!={r["trial_key"] for r in maps}:
        raise RuntimeError("Eligibility overlay does not cover exactly the original timing keys.")
    output=[]
    for original in maps:
        selection=lookup[original["trial_key"]]
        if csv_bool(original["included_by_trial_policy"])!=csv_bool(selection["original_included_by_trial_policy"]):
            raise RuntimeError("Application changed original Stage 2 eligibility.")
        row=dict(original)
        row.update({key:selection[key] for key in ("final_analysis_eligible","analysis_exclusion_reasons",
            "original_included_by_trial_policy","original_exclusion_reasons","saved_catch_type_code",
            "qualifies_saved_occlusion_type2","newly_included","newly_excluded","sequence_exclusion")})
        output.append(row)
    return output


# -----------------------------------------------------------------------------
# P24-only donor application; measured evidence is never removed or relabelled
# -----------------------------------------------------------------------------


def rebuild_p24_donors(maps,proposals,estimates,measurements,ledger,method,epoch_window):
    """Remove the five decided measured keys and rebuild only dependent estimates.

    P11 and all other participant donors are outside this function. The seven
    analysis-only sequence decisions never enter donor membership. Codes
    116/117 retain the existing insufficient-donor outcome and no numeric value.
    """
    by_key={r["event_key"]:r for r in measurements}
    if not set(P24_DONOR_KEYS).issubset(by_key):raise RuntimeError("Missing decided P24 donor evidence.")
    if any(by_key[k]["measurement_status"]!="measured" for k in P24_DONOR_KEYS):
        raise RuntimeError("P24 donor exclusion no longer identifies measured evidence.")
    result=copy.deepcopy(estimates);audit=[]
    for code in (81,111,112,114,115):
        eid=f"P24_code{code}";old=estimates[eid]
        donors=[r for r in measurements if donor_allowed("24",dict(condition_code=code,acquisition_run=1),r)
                and r["event_key"] not in P24_DONOR_KEYS]
        estimate=estimate_lag(donors,method)
        for name in ("target_context_limitation",):
            if name in old:estimate[name]=old[name]
        estimate.update(excluded_donor_event_keys=[k for k in P24_DONOR_KEYS if int(by_key[k]["trigger_code"])==code],
                        source_estimate_sha256=digest(old))
        if len(donors)!=old["n_donors"]-1:raise RuntimeError("P24 donor count did not remove exactly one event per affected code.")
        result[eid]=estimate
        audit.append(dict(estimate_id=eid,old_n_donors=old["n_donors"],new_n_donors=estimate["n_donors"],
            old_lag_ms=old.get("lag_estimate_ms"),new_lag_ms=estimate.get("lag_estimate_ms"),
            old_sd_ms=old.get("donor_sd_ms"),new_sd_ms=estimate.get("donor_sd_ms")))
    original_ledger={r["trial_key"]:r for r in ledger}
    new_maps=copy.deepcopy(maps);new_proposals=copy.deepcopy(proposals)
    for row in new_maps:
        eid=row.get("donor_estimate_id","")
        if row["assignment_basis"]!="estimated" or eid not in {f"P24_code{c}" for c in (81,111,112,114,115)}:continue
        estimate=result[eid];lag=estimate["lag_estimate_ms"];sf=float(row["sfreq_hz"])
        shift=round_half_away(lag*sf/1000);proposed=int(row["raw_sample_frame"])+shift
        row.update(lag_estimate_ms=lag,applied_shift_samples=shift,proposed_sample_raw=proposed,rounding_residual_samples=shift-lag*sf/1000)
        row.update(window_coverage(original_ledger[row["trial_key"]],proposed,*epoch_window,tuple(estimate["donor_observed_p05_p95_ms"])))
        row["proposed_recording_time_sec"]=(proposed-int(row["recording_first_sample_raw"]))/sf
    for row in new_proposals:
        eid=row.get("donor_estimate_id","")
        if row["assignment_basis"]!="estimated" or eid not in {f"P24_code{c}" for c in (81,111,112,114,115)}:continue
        lag=result[eid]["lag_estimate_ms"];shift=round_half_away(lag*float(row["sfreq_hz"])/1000)
        row.update(lag_estimate_ms=lag,applied_shift_samples=shift,proposed_event_sample_raw=int(row["digital_sample_raw"])+shift)
    for code in (116,117):
        estimate=result[f"P24_code{code}"]
        if estimate["status"]!="insufficient_donors" or estimate["n_donors"]!=6 or estimate.get("lag_estimate_ms") is not None:
            raise RuntimeError("P24 codes 116/117 must remain unresolved without a numeric estimate.")
    return new_maps,new_proposals,result,audit


# -----------------------------------------------------------------------------
# Hash-bound downstream validation; human acceptance remains a separate gate
# -----------------------------------------------------------------------------


def load_applied_eligibility(cfg,prepared,bundle_folder,bundle,review,application_path=None):
    """Verify accepted application hashes and recompute its narrow eligibility.

    Called only after the existing Stage 2 and accepted timing gates. The final
    bundle must bind application.json, saved metadata and overlay bytes. Review
    JSON must explicitly accept that application hash and its scientific scopes.
    No generated Task-1 artifact supplies this human acceptance automatically.
    """
    folder=Path(bundle_folder).resolve()
    path=Path(application_path).resolve() if application_path else folder/"application.json"
    if path!=folder/"application.json" or file_hash(path)!=bundle.get("review_application_sha256"):
        raise RuntimeError("Eligibility application is not bound to this timing bundle.")
    sha=file_hash(path)
    if review.get("accepted_eligibility_application_sha256")!=sha:
        raise RuntimeError("Timing review did not accept this exact eligibility application.")
    application=json.loads(path.read_text())
    if (application.get("schema_version")!=APPLICATION_SCHEMA or application.get("policy")!=policy_definition()
            or application.get("subject_bids_id")!=cfg["layout"]["subject_bids_id"]
            or application.get("target_version")!=bundle["version"]):
        raise RuntimeError("Unsupported or mismatched eligibility application.")
    # Derive the minimum scope from the recomputed overlay and participant.
    # A bundle hash binds what the application claims; it cannot make a missing
    # human-review requirement optional by omitting that claim.
    overlay=build_eligibility_overlay(prepared,read_csv(folder/"saved_trial_metadata.csv"),cfg["subject"]["id"])
    required=mandatory_review_scopes(overlay,cfg["subject"]["id"])
    declared=set(application.get("scientific_review_requirements",[]))
    if not required.issubset(declared):
        raise RuntimeError("Eligibility application omits a mandatory scientific review scope.")
    if not required.issubset(review.get("accepted_review_scopes",[])):
        raise RuntimeError("New-catch/final-eligibility or participant scientific scope is not explicitly accepted.")
    if application["stage2_report_sha256"]!=bundle["stage2_report_sha256"] or file_hash(stage2_paths(cfg)["ledger"])!=application["stage2_ledger_sha256"]:
        raise RuntimeError("Eligibility application changed its original Stage 2 accounting.")
    source=Path(cfg["resolved_paths"]["derivatives_root"])/"visual_timing"/application["source_version"]/cfg["layout"]["subject_bids_id"]/"ses-01"
    if file_hash(source/"manifest.json")!=application["source_manifest_sha256"]:
        raise RuntimeError("Eligibility application's original Stage 3 source changed.")
    source_manifest=json.loads((source/"manifest.json").read_text())
    if file_hash(folder/"event_measurements.csv")!=source_manifest["output_sha256"]["event_measurements.csv"]:
        raise RuntimeError("Application changed original measured timing evidence.")
    for name,expected in application["saved_mat_sha256"].items():
        if file_hash(Path(cfg["resolved_paths"]["subject_data_dir"])/"output_files"/name)!=expected:
            raise RuntimeError("Saved catch-class evidence changed after eligibility application.")
    for name,field in (("saved_trial_metadata.csv","saved_metadata_sha256"),("eligibility_overlay.csv","eligibility_overlay_sha256")):
        if file_hash(folder/name)!=application[field]:raise RuntimeError("Eligibility evidence hash mismatch.")
    expected=overlay
    if csv_text(expected).encode()!=(folder/"eligibility_overlay.csv").read_bytes():
        raise RuntimeError("Eligibility overlay contradicts original reasons, saved class or exact exclusions.")
    donor_keys=list(P24_DONOR_KEYS) if cfg["subject"]["id"]=="24" else []
    if application["donor_excluded_event_keys"]!=donor_keys:
        raise RuntimeError("Application changed the independent decided donor exclusions.")
    return {r["trial_key"]:r for r in expected},dict(application_path=str(path),application_sha256=sha,
        eligibility_overlay_sha256=application["eligibility_overlay_sha256"],policy_id=POLICY_ID,
        donor_excluded_event_keys=donor_keys,accepted_review_scopes=sorted(required))


def mandatory_review_scopes(overlay,subject):
    """Return review scopes required by recomputed facts, independent of claims."""
    scopes={"final_analysis_eligibility"}
    if any(bool(row["newly_included"]) for row in overlay):
        scopes.add("newly_included_catch_onsets")
    if str(subject)=="33":
        scopes.add("p33_remaining_original_onsets")
    return scopes


def analysis_counts(ledger,overlay,maps):
    """Count original/final eligibility by source, condition, identity and replay."""
    selections={r["trial_key"]:r for r in overlay};timing={r["trial_key"]:r for r in maps}
    groups=defaultdict(list)
    for row in ledger:
        for factor in ("source_fif","condition_code","trial_sequence_identity","is_replay_trial","behavioral_block"):
            groups[(factor,str(row[factor]))].append(row)
        groups[("condition_identity",f"{row['condition_code']}:{row['trial_sequence_identity']}")].append(row)
        groups[("source_condition",f"{row['source_fif']}:{row['condition_code']}")].append(row)
        groups[("source_block_condition_replay",f"{row['source_fif']}:{row['behavioral_block']}:{row['condition_code']}:{row['is_replay_trial']}")].append(row)
    output=[]
    for (factor,level),group in sorted(groups.items()):
        final=[r for r in group if selections[r["trial_key"]]["final_analysis_eligible"]]
        output.append(dict(factor=factor,level=level,n_recorded=len(group),
            n_original_eligible=sum(csv_bool(r["included_by_trial_policy"]) for r in group),n_final_eligible=len(final),
            n_newly_included=sum(selections[r["trial_key"]]["newly_included"] for r in group),
            n_newly_excluded=sum(selections[r["trial_key"]]["newly_excluded"] for r in group),
            n_final_unassigned=sum(timing[r["trial_key"]].get("proposed_sample_raw") in (None,"") for r in final),
            n_final_window_outside=sum(timing[r["trial_key"]].get("proposed_sample_raw") not in (None,"") and
                timing[r["trial_key"]].get("proposed_window_inside_recording") not in (None,"") and
                not csv_bool(timing[r["trial_key"]]["proposed_window_inside_recording"]) for r in final),
            n_final_coverage_unknown=sum(timing[r["trial_key"]].get("proposed_window_inside_recording") in (None,"") for r in final),
            n_final_trial_end_overlap=sum(csv_bool(timing[r["trial_key"]].get("epoch_overlaps_trial_interval_end",False)) for r in final)))
    return output


def short_recording_checks(subject,ledger,overlay):
    """Protect the accepted 37- and 16-trial short recordings and block identity."""
    specifications={"29":("sub29_run02_150926_halfBlock.fif",37,44),"38":("sub38_run01_051026.fif",16,20)}
    if subject not in specifications:return []
    source,eligible,total=specifications[subject]
    selected={r["trial_key"]:r for r in overlay};group=[r for r in ledger if r["source_fif"]==source]
    if len(group)!=total or sum(selected[r["trial_key"]]["final_analysis_eligible"] for r in group)!=eligible:
        raise RuntimeError("Accepted short-recording eligibility changed.")
    if 6 not in {int(r["behavioral_block"]) for r in ledger}:
        raise RuntimeError("Short-recording participant lost separate block-6 identity.")
    if subject=="29" and sum("catch_type_1" in r["exclusion_reasons"].split(",") for r in group)!=6:
        raise RuntimeError("P29 short recording lost its six type-1 exclusions.")
    if subject=="29":
        terminal=[r for r in group if int(r["raw_sample_frame"])==269558]
        if len(terminal)!=1 or selected[terminal[0]["trial_key"]]["final_analysis_eligible"] or not {
                "esc_termination","stage1_reviewed_interrupted_trial"}.issubset(terminal[0]["exclusion_reasons"].split(",")):
            raise RuntimeError("P29 short recording lost its exact interrupted final-trial exclusion.")
    return [dict(source_fif=source,n_recorded=total,n_final_eligible=eligible,n_independent_exclusions=total-eligible,
        separate_behavioral_blocks=sorted({int(r["behavioral_block"]) for r in group}),block6_identity_preserved=True)]


# -----------------------------------------------------------------------------
# Staging publication and explicit Task-2 API (never final s3-v002 publication)
# -----------------------------------------------------------------------------


def stage_subject(config_path,subject,source_version,target_version,review_log,staging_root,*,validate_only=False,site_config=None):
    """Build one complete Task-1 staging result from verified existing evidence."""
    implementation_files=[Path(__file__),Path(__file__).with_name("visual_timing_utils.py"),
        Path(__file__).with_name("visual_timing_report.py"),Path(__file__).with_name("04_epoching.py")]
    implementation_bytes={path.name:path.read_bytes() for path in implementation_files}
    implementation_hashes={name:hashlib.sha256(value).hexdigest() for name,value in implementation_bytes.items()}
    cfg,warnings,_,_=load_and_resolve_config(config_path,Path(config_path).parent/f"sub-{subject}.yaml",site_config_path=site_config)
    if warnings:raise RuntimeError(f"Configuration warnings: {warnings}")
    root=Path(cfg["resolved_paths"]["derivatives_root"])/"visual_timing"
    source=root/source_version/f"sub-{subject}"/"ses-01"
    target=Path(staging_root)/f"sub-{subject}"/"ses-01"
    if target.resolve().is_relative_to(source.resolve()):raise ValueError("Staging cannot overwrite its immutable source.")
    manifest=verify_bundle(source);source_sha=file_hash(source/"manifest.json")
    if manifest["mode"]!="cohort_proposal" or manifest["subject_bids_id"]!=f"sub-{subject}":
        raise RuntimeError("Task 1 requires the matching complete cohort source.")
    ledger=load_prepared_decisions(cfg);paths=stage2_paths(cfg)
    if file_hash(paths["report"])!=manifest["stage2_report_sha256"]:
        raise RuntimeError("Stage 2 completion record differs from the timing source.")
    contexts=json.loads((source/"contexts.json").read_text())
    saved=read_saved_metadata(cfg,manifest,contexts)
    overlay=build_eligibility_overlay(ledger,saved,subject)
    maps=read_csv(source/"trial_timing_map.csv");proposals=read_csv(source/"event_timing_proposals.csv")
    if {r["trial_key"] for r in maps}!={r["trial_key"] for r in ledger}:raise RuntimeError("Stage 2/timing key sets differ.")
    estimates=json.loads((source/"estimates.json").read_text());donor_audit=[]
    if subject=="24":
        from visual_timing_report import normalize
        maps,proposals,estimates,donor_audit=rebuild_p24_donors(maps,proposals,estimates,
            normalize(read_csv(source/"event_measurements.csv")),ledger,manifest["method"],manifest["epoch_window_sec"])
    maps=attach_eligibility(maps,overlay)
    counts=analysis_counts(ledger,overlay,maps);short=short_recording_checks(subject,ledger,overlay)
    recorded={context_key(r) for r in ledger};unrecorded=[]
    unavailable={context_key(r):r for r in read_csv(paths["unavailable"])}
    for metadata in saved:
        key=context_key(metadata)
        if metadata["saved_catch_type_code"]!=2 or metadata["saved_condition_code"] not in (81,91) or key in recorded:continue
        if key not in unavailable:raise RuntimeError("Unmatched saved type-2 trial is not an accepted unavailable tail.")
        unrecorded.append(dict(metadata,source_fif=unavailable[key]["constrained_source_fif"],recorded_onset=False,trial_key="",
                               reason="accepted_recording_tail_unavailable;no_onset_fabricated"))
    summary=dict(subject_bids_id=f"sub-{subject}",n_original_keys=len(ledger),
        n_original_eligible=sum(r["original_included_by_trial_policy"] for r in overlay),
        n_final_eligible=sum(r["final_analysis_eligible"] for r in overlay),
        n_newly_included=sum(r["newly_included"] for r in overlay),n_newly_excluded=sum(r["newly_excluded"] for r in overlay),
        n_saved_type2_recorded=sum(r["qualifies_saved_occlusion_type2"] for r in overlay),
        n_catch_only_before=sum(r["original_exclusion_reasons"]=="catch_type_2" and r["qualifies_saved_occlusion_type2"] for r in overlay),
        n_sequence_excluded=sum(r["sequence_exclusion"] for r in overlay),n_unrecorded_saved_type2=len(unrecorded),
        n_recorded_without_saved_metadata=sum(not r["source_mat"] for r in overlay),
        short_recording_checks=short,p24_estimate_changes=donor_audit)
    if validate_only:return summary
    saved_text=csv_text(saved);overlay_text=csv_text(overlay)
    scopes=["final_analysis_eligibility"]
    if summary["n_newly_included"]:scopes.append("newly_included_catch_onsets")
    if subject=="33":scopes.append("p33_remaining_original_onsets")
    log_bytes=Path(review_log).read_bytes()
    application=dict(schema_version=APPLICATION_SCHEMA,status="applied_pending_scientific_reconciliation",
        subject_bids_id=f"sub-{subject}",source_version=source_version,target_version=target_version,policy=policy_definition(),
        source_manifest_sha256=source_sha,stage2_report_sha256=manifest["stage2_report_sha256"],
        stage2_ledger_sha256=file_hash(paths["ledger"]),saved_mat_sha256=manifest["saved_mat_sha256"],
        saved_metadata_sha256=hashlib.sha256(saved_text.encode()).hexdigest(),
        eligibility_overlay_sha256=hashlib.sha256(overlay_text.encode()).hexdigest(),
        review_log_sha256=hashlib.sha256(log_bytes).hexdigest(),review_log_source=str(Path(review_log)),
        donor_excluded_event_keys=list(P24_DONOR_KEYS) if subject=="24" else [],
        scientific_review_requirements=scopes,summary=summary)
    application_text=json.dumps(application,indent=2,sort_keys=True,allow_nan=False)
    application_sha=hashlib.sha256(application_text.encode()).hexdigest()
    selection={r["trial_key"]:r for r in overlay};measurement_lookup={r["event_key"]:r for r in read_csv(source/"event_measurements.csv")}
    new_qc=[]
    for row in maps:
        if not selection[row["trial_key"]]["newly_included"]:continue
        measurement_key=row.get("measurement_event_key") or f"sub-{subject}:{row['source_fif']}:{row['raw_sample_frame']}:{row['condition_code']}"
        measurement=measurement_lookup.get(measurement_key,{})
        new_qc.append(dict(row,measurement_status=measurement.get("measurement_status","not_measured"),
            measurement_qc_reasons=measurement.get("qc_reasons",""),threshold_span_ms=measurement.get("threshold_span_ms",""),
            scientific_scope="new_catch_not_yet_accepted"))
    fixed=("event_measurements.csv","recording_summary.csv","physical_intervals.csv","contexts.json",
           "diagnostics.png","waveform_windows.npz","diagnostic_index.csv")
    artifacts={name:(source/name).read_bytes() for name in fixed}
    artifacts.update({"source_manifest.json":(source/"manifest.json").read_bytes(),
        "original_stage2_ledger.csv":paths["ledger"].read_bytes(),"review_notes_snapshot.md":log_bytes,
        "saved_trial_metadata.csv":saved_text,"eligibility_overlay.csv":overlay_text,"application.json":application_text,
        "trial_timing_map.csv":csv_text(maps),"event_timing_proposals.csv":csv_text(proposals) if subject=="24" else (source/"event_timing_proposals.csv").read_bytes(),
        "estimates.json":json.dumps(estimates,indent=2,sort_keys=True,allow_nan=False) if subject=="24" else (source/"estimates.json").read_bytes(),
        "donor_exclusions.csv":csv_text([dict(event_key=k,reason="P24_whole_trial_sequence_exclusion",analysis_eligibility_filter=False) for k in application["donor_excluded_event_keys"]],("event_key","reason","analysis_eligibility_filter")),
        "p24_estimate_changes.csv":csv_text(donor_audit,("estimate_id","old_n_donors","new_n_donors")),
        "analysis_counts.csv":csv_text(counts),"new_catch_timing_qc.csv":csv_text(new_qc,("trial_key","final_analysis_eligible","assignment_basis")),
        "unrecorded_saved_type2.csv":csv_text(unrecorded,("source_mat","behavioral_block","behavioral_run","behavioral_trial","recorded_onset","trial_key")),
        "application_summary.json":json.dumps(summary,indent=2,sort_keys=True),
        "code_visual_timing_application.py.txt":implementation_bytes[Path(__file__).name]})
    artifacts.update({f"code_{name}.txt":value for name,value in implementation_bytes.items()})
    text=(f"# Task 1 staging: P{subject}\n\nNot final s3-v002; scientific acceptance is pending.\n\n"
        f"Original keys: {len(ledger)}. Original eligible: {summary['n_original_eligible']}. Final analysis eligible: {summary['n_final_eligible']}. "
        f"Newly included: {summary['n_newly_included']}; newly excluded: {summary['n_newly_excluded']}; whole-trial sequence exclusions: {summary['n_sequence_excluded']}.\n\n"
        "Saved class-2 condition-81/91 catches use ordinary motion onset; only catch_type_2 is waived. Independent exclusions and exact sequence decisions remain. "
        "Original Stage 2 eligibility stays separate. Measured evidence is retained; only the five exact P24 keys change donor inputs. "
        "New catch timing/QC/coverage and P33's remaining original onset scope require scientific reconciliation. "
        "P04-P07/P17 timing estimates belong to Task 2. No accepted review.json is created.\n\n"
        "See application.json, eligibility_overlay.csv, analysis_counts.csv, new_catch_timing_qc.csv, donor_exclusions.csv and p24_estimate_changes.csv. "
        "Unavailable MAT-only outputs have no onset or trial key. No epoching, exports or analyses ran.\n")
    artifacts["report.md"]=text
    downloads="".join(f"<li><a download='{name}' href='data:text/csv;base64,{__import__('base64').b64encode(value.encode() if isinstance(value,str) else value).decode()}'>{name}</a></li>"
        for name,value in artifacts.items() if name.endswith('.csv'))
    artifacts["report.html"]="<!doctype html><html><meta charset='utf-8'><title>Task 1 staging P"+subject+"</title><body><pre style='white-space:pre-wrap;font:17px/1.6 system-ui'>"+html.escape(text)+"</pre><ul>"+downloads+"</ul></body></html>"
    new_manifest=dict(manifest,version=target_version,mode="task1_staging",task1_staging=True,
        generated_utc=utc_now(),scientific_acceptance="pending_independent_review",review_application_sha256=application_sha,
        source_version=source_version,source_manifest_sha256=source_sha,application_summary=summary,
        eligible_assignment_counts=dict(Counter(r["assignment_basis"] for r in maps if csv_bool(r["final_analysis_eligible"]))),
        source_n_signal_windows_read=manifest.get("n_signal_windows_read"),n_signal_windows_read=0,
        implementation_sha256=implementation_hashes)
    new_manifest["input_contract"]=dict(source_manifest_sha256=source_sha,
        source_input_contract_sha256=manifest.get("input_contract_sha256"),
        stage2_report_sha256=manifest["stage2_report_sha256"],review_application_sha256=application_sha,
        implementation_sha256=new_manifest["implementation_sha256"])
    new_manifest["input_contract_sha256"]=digest(new_manifest["input_contract"])
    # P24's report must reflect its actual filtered donor estimates and final
    # eligibility. Other Task-1 reports are compact handoffs; Task 2 produces the
    # final full-cohort presentations after its method is applied coherently.
    if subject=="24":
        from visual_timing_report import build_report_artifacts,load_cached_waveforms
        reporting=build_report_artifacts(read_csv(source/"event_measurements.csv"),maps,estimates,
            read_csv(source/"recording_summary.csv"),read_csv(source/"physical_intervals.csv"),contexts,new_manifest,
            load_cached_waveforms(source),source_identity=dict(source_manifest_sha256=source_sha,review_application_sha256=application_sha),proposals=proposals)
        for name,value in reporting.items():
            if name not in {"event_measurements.csv","recording_summary.csv","physical_intervals.csv","contexts.json","estimates.json","event_timing_proposals.csv","trial_timing_map.csv"}:
                artifacts[name]=value
    verify_bundle(source)
    if file_hash(source/"manifest.json")!=source_sha:raise RuntimeError("Source changed during Task-1 staging.")
    if any(file_hash(path)!=implementation_hashes[path.name] for path in implementation_files):
        raise RuntimeError("Application/consumer/report implementation changed during participant staging.")
    from visual_timing import publish
    publish(target,artifacts,new_manifest,cfg["resolved_paths"]["data_root"])
    return summary


def load_task1_stage(staging_root,subject):
    """Task-2 API: verify and return (manifest, filename->bytes) for one subject.

    Keep application.json, its overlay/metadata/source snapshots and original
    included_by_trial_policy intact. Task 2 may change only its supported timing
    method/estimates/proposals/maps and presentations. Final mode must become
    cohort_proposal with task1_staging=False, preserving application hash and
    all original source/Stage-2 checks. No acceptance record is implied.
    """
    root=Path(staging_root)
    handoff=json.loads((root/"task1_handoff.json").read_text())
    if (handoff.get("status")!="task1_staged_awaiting_task2" or
            {r["subject_bids_id"] for r in handoff["participants"]}!={f"sub-{s}" for s in ANALYSIS_SUBJECT_IDS}):
        raise RuntimeError("Task 2 requires the complete coordinated 36-participant Task-1 checkpoint.")
    folder=root/f"sub-{subject}"/"ses-01"
    if file_hash(folder/"manifest.json")!=handoff["staged_manifest_sha256"][f"sub-{subject}"]:
        raise RuntimeError("Task-1 participant changed after cohort checkpoint sealing.")
    manifest=verify_bundle(folder)
    if manifest["mode"]!="task1_staging" or not manifest.get("task1_staging"):
        raise RuntimeError("Not a complete Task-1 staging bundle.")
    if file_hash(folder/"application.json")!=manifest["review_application_sha256"]:
        raise RuntimeError("Task-1 application fingerprint mismatch.")
    return manifest,{name:(folder/name).read_bytes() for name in manifest["output_sha256"]}


def main():
    """Validate/stage P04-P39, reserving final publication for Task 2."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",default="meg/params/pipeline_default.yaml",type=Path)
    parser.add_argument("--site-config",type=Path)
    parser.add_argument("--source-version",default="s3-v001")
    parser.add_argument("--target-version",default="s3-v002")
    parser.add_argument("--staging-root",type=Path)
    parser.add_argument("--review-log",type=Path)
    parser.add_argument("--subjects",nargs="+",default=list(ANALYSIS_SUBJECT_IDS))
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--validate-only",action="store_true")
    args=parser.parse_args()
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*",value) for value in (args.source_version,args.target_version)) or args.source_version==args.target_version:
        raise ValueError("Source and target must be distinct safe version names.")
    cfg,*_=load_and_resolve_config(args.config,args.config.parent/"sub-04.yaml",site_config_path=args.site_config)
    root=Path(cfg["resolved_paths"]["derivatives_root"])/"visual_timing"
    staging=args.staging_root or root/f"{args.target_version}-task1-staging"
    review_log=args.review_log or root/args.source_version/"review_notes.md"
    if staging.resolve()==(root/args.target_version).resolve():raise ValueError("Task 1 cannot freeze the final coordinated version.")
    if not staging.resolve().is_relative_to(root.resolve()) or staging.resolve().is_relative_to((root/args.source_version).resolve()):
        raise ValueError("Choose new staging below derivatives/visual_timing, outside the source version.")
    guarded_path(staging,cfg["resolved_paths"]["data_root"])
    if staging.exists() and not args.validate_only:raise ValueError("Choose a fresh Task-1 staging directory; overwrite is not supported.")
    if args.workers<1 or args.workers>4:raise ValueError("Use one to four metadata/report workers.")
    results=[]
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context("spawn")) as pool:
        futures=[pool.submit(stage_subject,str(args.config.resolve()),subject,args.source_version,args.target_version,
            str(review_log.resolve()),str(staging.resolve()),validate_only=args.validate_only,site_config=args.site_config) for subject in args.subjects]
        for future in as_completed(futures):
            try:result=future.result()
            except BaseException:
                for pending in futures:pending.cancel()
                raise
            results.append(result)
            print(json.dumps(result),flush=True)
    results.sort(key=lambda r:r["subject_bids_id"])
    summary=dict(status="validated" if args.validate_only else "task1_staged_awaiting_task2",source_version=args.source_version,
        target_version=args.target_version,staging_root=str(staging),participants=results,
        totals={field:sum(r[field] for r in results) for field in ("n_original_keys","n_original_eligible","n_final_eligible",
            "n_newly_included","n_newly_excluded","n_saved_type2_recorded","n_catch_only_before","n_sequence_excluded","n_unrecorded_saved_type2","n_recorded_without_saved_metadata")},
        api="visual_timing_application.load_task1_stage(staging_root, subject)",
        original_eligibility_field="included_by_trial_policy",analysis_eligibility_field="final_analysis_eligible",
        application_file="application.json",donor_excluded_event_keys=list(P24_DONOR_KEYS),
        task2_targets=["04","05","06","07","17"],final_version_not_published=True,
        no_accepted_review_created=True,raw_signal_reads=0,scientific_acceptance="pending_independent_review")
    if not args.validate_only:
        staging.mkdir(parents=True,exist_ok=True)
        summary["staged_manifest_sha256"]={r["subject_bids_id"]:file_hash(staging/r["subject_bids_id"]/"ses-01/manifest.json") for r in results}
        logs={json.loads((staging/r["subject_bids_id"]/"ses-01/application.json").read_text())["review_log_sha256"] for r in results}
        if len(logs)!=1:raise RuntimeError("Task-1 participants used different working-log snapshots.")
        summary["review_log_snapshot_sha256"]=next(iter(logs))
        summary["policy_sha256"]=digest(policy_definition())
        (staging/"task1_handoff.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
        flat=[{k:v for k,v in row.items() if k not in ("short_recording_checks","p24_estimate_changes")} for row in results]
        (staging/"participant_counts.csv").write_text(csv_text(flat))
        changes="\n".join(f"| {r['subject_bids_id']} | {r['n_original_eligible']} | {r['n_final_eligible']} | {r['n_newly_included']} | {r['n_newly_excluded']} |" for r in results)
        note=("# Task 1 handoff: staged application for coordinated s3-v002\n\n"
            "Status: Task 1 staged; Task 2 must add its onset method and publish the one final s3-v002. No accepted review.json exists.\n\n"
            "| Participant | Original eligible | Final policy eligible | Added catches | Removed original trials |\n|---|---:|---:|---:|---:|\n"+changes+
            "\n\nOriginal Stage 2 eligibility remains in included_by_trial_policy; final_analysis_eligible controls reviewed analysis selection. "
            "All eight exact sequence decisions apply; only P24 removes donor events. P29's 37 short-recording trials and P38's 16 remain separately identified. "
            "Five P37 saved type-2 outputs remain unrecorded; absent saved classes never authorize inclusion.\n\n"
            "API: visual_timing_application.load_task1_stage(staging_root, subject) returns a verified manifest and filename-to-bytes dictionary. "
            "Keep application.json, eligibility_overlay.csv, saved_trial_metadata.csv and source snapshots unchanged. "
            "Preserve event_measurements.csv byte-for-byte; supply Task-2 cross-participant donor evidence as a separate artifact. "
            "Do not filter donors using final analysis eligibility; remove only the five P24 donor_excluded_event_keys listed in task1_handoff.json. "
            "Task 2 must extend the estimated-assignment validator for its five targets, rebuild final presentations/counts/coverage, update the coherent method/input contract, "
            "set mode=cohort_proposal and task1_staging=False, and publish all 36 into the fresh final version. "
            "Keep review_application_sha256 bound to the frozen application.\n\n"
            "Scientific holds: newly included catches need timing/QC/coverage acceptance; P33's remaining original onsets need explicit reconciliation; "
            "Task-2 numerical estimates/donor provenance still need review. Per-participant accepted timing records must name the bundle hash, "
            "accepted_eligibility_application_sha256 and accepted_review_scopes from application.json. No automatic acceptance, cleaning, epoching, export or analysis is authorized by this handoff.\n")
        (staging/"TASK1_HANDOFF.md").write_text(note)
        rows="".join(f"<tr><td><a href='{r['subject_bids_id']}/ses-01/report.html'>{r['subject_bids_id']}</a></td><td>{r['n_original_eligible']}</td><td>{r['n_final_eligible']}</td><td>{r['n_newly_included']}</td><td>{r['n_newly_excluded']}</td></tr>" for r in results)
        (staging/"index.html").write_text("<!doctype html><html lang='en'><meta charset='utf-8'><title>Task 1 Stage 3 staging</title><body style='font:17px/1.6 system-ui;max-width:1100px;margin:35px auto'><h1>Task 1: staged review application</h1><p>Awaiting Task 2 and scientific reconciliation. The final s3-v002 is not published here. Original Stage 2 flags and timing evidence remain separate and unchanged.</p><p><a href='TASK1_HANDOFF.md'>Task 2 handoff</a> | <a href='task1_handoff.json'>Exact checkpoint hashes</a> | <a href='participant_counts.csv'>Counts CSV</a></p><table><tr><th>Report</th><th>Original eligible</th><th>Final policy eligible</th><th>Added catches</th><th>Removed original trials</th></tr>"+rows+"</table></body></html>")
    print(json.dumps(summary),flush=True)
    return 0


if __name__=="__main__":raise SystemExit(main())
