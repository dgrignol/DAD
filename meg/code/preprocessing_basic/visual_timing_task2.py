#!/usr/bin/env python3
"""Build the coordinated Stage 3 cohort timing proposal from sealed Task 1.

Usage from the DAD root (reads saved Stage 3/Stage 2 artifacts, never raw FIF):

    python meg/code/preprocessing_basic/visual_timing_task2.py \
      --config meg/params/pipeline_default.yaml --workers 2 --validate-only
    python meg/code/preprocessing_basic/visual_timing_task2.py \
      --config meg/params/pipeline_default.yaml --workers 2

The script verifies all sealed Task 1 bundles, derives one ordinary median from
each trigger's unique usable measured cohort events, applies it only to the
five decided no-diode participants, rebuilds reports and publishes a complete
new version atomically. Donor selection is independent of trial eligibility;
only the exact five P24 keys in the handoff are excluded. Source measurements,
application artifacts, Stage 2 accounting and accepted source maps are retained.
No signal windows are reread and no review acceptance is created.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import visual_timing_application as application
import visual_timing as timing
from pipeline_config import ANALYSIS_SUBJECT_IDS,load_and_resolve_config
from preprocessing_stage_utils import utc_now
from trial_preparation import csv_text,digest,file_hash,read_csv
from visual_timing_utils import round_half_away,verify_bundle
from downstream_stage_utils import csv_bool

TRIGGERS=(71,81,91)
TARGETS={"04","05","06","07","17"}
METHOD_ID="ordinary_median_of_unique_qc_usable_measured_cohort_events_v1"


def rows_from_bytes(content):
    """Parse a saved CSV byte stream without changing its original bytes."""
    return list(csv.DictReader(io.StringIO(content.decode("utf-8-sig"))))


def pooled_donors(staging_root,handoff):
    """Collect each genuine measured onset once, bound to its source hashes."""
    excluded=set(handoff["donor_excluded_event_keys"])
    donors=[]; sources={}; seen=set()
    for participant in handoff["participants"]:
        subject=participant["subject_bids_id"].replace("sub-","")
        manifest,artifacts=application.load_task1_stage(staging_root,subject)
        rows=rows_from_bytes(artifacts["event_measurements.csv"])
        for r in rows:
            if int(r["trigger_code"]) not in TRIGGERS or r["measurement_status"]!="measured" or r.get("qc_reasons","").strip():
                continue
            key=r["event_key"]
            if key in excluded:continue
            if key in seen:raise RuntimeError(f"Duplicate cohort donor identity: {key}")
            seen.add(key)
            sf=float(r["sfreq_hz"]);edge=int(r["diode_edge_sample_raw"]);anchor=int(r["digital_sample_raw"])
            lag=(edge-anchor)*1000/sf
            if not math.isclose(lag,float(r["lag_ms"]),rel_tol=0,abs_tol=1e-8):
                raise RuntimeError(f"Donor sample/lag arithmetic failed: {key}")
            donors.append({"event_key":key,"subject_bids_id":r["subject_bids_id"],"source_fif":r["source_fif"],
                "trigger_code":int(r["trigger_code"]),"digital_sample_raw":anchor,"diode_edge_sample_raw":edge,
                "sfreq_hz":sf,"lag_ms":lag,"behavioral_block":r.get("behavioral_block",""),
                "behavioral_run":r.get("behavioral_run",""),"acquisition_run":r.get("acquisition_run",""),
                "runtime":r.get("runtime",""),"qc_reasons":""})
        source_manifest=json.loads(artifacts["source_manifest.json"])
        sources[participant["subject_bids_id"]]={"source_version":manifest["source_version"],
            "source_manifest_sha256":manifest["source_manifest_sha256"],
            "event_measurements_sha256":hashlib.sha256(artifacts["event_measurements.csv"]).hexdigest(),
            "source_event_measurements_sha256":source_manifest["output_sha256"]["event_measurements.csv"]}
        if sources[participant["subject_bids_id"]]["event_measurements_sha256"]!=sources[participant["subject_bids_id"]]["source_event_measurements_sha256"]:
            raise RuntimeError("Task 1 event measurements differ from the sealed original source.")
    donors.sort(key=lambda r:r["event_key"])
    summaries={}
    expected={71:7341,81:7422,91:7357}
    for code in TRIGGERS:
        values=np.asarray([r["lag_ms"] for r in donors if r["trigger_code"]==code],float)
        if len(values)!=expected[code]:raise RuntimeError(f"Trigger {code}: expected {expected[code]} measured donors, found {len(values)}")
        summaries[str(code)]={"n":len(values),"median_ms":float(np.median(values)),"mean_ms":float(np.mean(values)),
            "sd_ms":float(np.std(values,ddof=1)),"minimum_ms":float(np.min(values)),"maximum_ms":float(np.max(values)),
            "q25_ms":float(np.quantile(values,.25)),"q75_ms":float(np.quantile(values,.75)),
            "p05_ms":float(np.quantile(values,.05)),"p95_ms":float(np.quantile(values,.95)),
            "estimator":METHOD_ID,"transfer_uncertainty":"unknown;observed cohort spread is not target transfer error"}
        if summaries[str(code)]["median_ms"]!=19.0:raise RuntimeError(f"Unexpected cohort median for trigger {code}")
    return donors,summaries,sources


def load_pooled_cache(build,handoff):
    """Resume only from the hash-bound donor evidence produced by this run."""
    donor_path=build/"cohort_pooled_donor_measurements.csv"
    provenance_path=build/"cohort_pooled_donor_provenance.json"
    if not donor_path.is_file():
        donor_path=build/"sub-04"/"ses-01"/"cohort_pooled_donor_measurements.csv"
        provenance_path=build/"sub-04"/"ses-01"/"cohort_pooled_donor_provenance.json"
    if not donor_path.is_file() or not provenance_path.is_file():return None
    donor_bytes=donor_path.read_bytes();provenance=json.loads(provenance_path.read_text())
    if hashlib.sha256(donor_bytes).hexdigest()!=provenance.get("donor_csv_sha256"):
        raise RuntimeError("Temporary cohort donor cache hash mismatch.")
    if provenance.get("donor_excluded_event_keys")!=handoff["donor_excluded_event_keys"] or provenance.get("analysis_eligibility_filter") is not False:
        raise RuntimeError("Temporary cohort donor cache does not match the sealed handoff.")
    donors=rows_from_bytes(donor_bytes);summaries={}
    expected={71:7341,81:7422,91:7357}
    for code,n in expected.items():
        values=np.asarray([float(r["lag_ms"]) for r in donors if int(r["trigger_code"])==code],float)
        if len(values)!=n or float(np.median(values))!=19.:raise RuntimeError("Temporary donor cache summary mismatch.")
        summaries[str(code)]={"n":n,"median_ms":19.,"mean_ms":float(np.mean(values)),"sd_ms":float(np.std(values,ddof=1)),
            "minimum_ms":float(np.min(values)),"maximum_ms":float(np.max(values)),"q25_ms":float(np.quantile(values,.25)),
            "q75_ms":float(np.quantile(values,.75)),"p05_ms":float(np.quantile(values,.05)),"p95_ms":float(np.quantile(values,.95)),
            "estimator":METHOD_ID,"transfer_uncertainty":"unknown;observed cohort spread is not target transfer error"}
    return donors,summaries,provenance.get("donor_sources",{})


def apply_pool(subject,maps,proposals,estimates,summaries,tmin,tmax):
    """Apply native-sample-rounded medians only to the five target participants."""
    if subject not in TARGETS:return maps,proposals,estimates
    lookup={str(code):summary for code,summary in summaries.items()}
    for row in maps:
        code=int(row["condition_code"])
        if code not in TRIGGERS or row["assignment_basis"]!="digital_fallback":continue
        summary=lookup[str(code)];sf=float(row["sfreq_hz"]);lag=summary["median_ms"]
        shift=round_half_away(lag*sf/1000);proposed=int(row["raw_sample_frame"])+shift
        row.update(assignment_basis="estimated",assignment_status="proposed",donor_estimate_id=f"COHORT_code{code}",
            lag_estimate_ms=lag,applied_shift_samples=shift,proposed_sample_raw=proposed,
            uncertainty_kind="observed_cohort_dispersion;target_transfer_unknown",
            uncertainty_ms=None,donor_provenance="cohort_pooled_donor_provenance.json")
        from visual_timing_utils import window_coverage
        cov=window_coverage(row,proposed,tmin,tmax)
        row.update(cov)
        row["proposed_recording_time_sec"]=(proposed-int(row["recording_first_sample_raw"]))/sf
    maps_by={r["trial_key"]:r for r in maps}
    for row in proposals:
        if int(row["trigger_code"]) in TRIGGERS and row.get("assignment_basis")=="digital_fallback":
            # These target proposals correspond to onset events; map coverage is
            # exact-key and the map remains authoritative for trial epoching.
            lag=lookup[str(int(row["trigger_code"]))]["median_ms"]
            shift=round_half_away(lag*float(row["sfreq_hz"])/1000)
            row.update(assignment_basis="estimated",donor_estimate_id=f"COHORT_code{int(row['trigger_code'])}",
                lag_estimate_ms=lag,applied_shift_samples=shift,
                proposed_event_sample_raw=int(row["digital_sample_raw"])+shift,
                uncertainty_kind="observed_cohort_dispersion;target_transfer_unknown")
    for code,s in lookup.items():estimates[f"COHORT_code{code}"]=dict(status="estimated",estimator="cohort_event_median",
        method_id=METHOD_ID,
        lag_estimate_ms=s["median_ms"],n_donors=s["n"],donor_sd_ms=s["sd_ms"],
        donor_observed_min_max_ms=[s["minimum_ms"],s["maximum_ms"]],donor_observed_q25_q75_ms=[s["q25_ms"],s["q75_ms"]],
        donor_observed_p05_p95_ms=[s["p05_ms"],s["p95_ms"]],transfer_uncertainty=s["transfer_uncertainty"])
    return maps,proposals,estimates


def final_derived_tables(artifacts,maps,subject,summaries):
    """Refresh cohort accounting from the final map, preserving Task 1 evidence."""
    ledger=rows_from_bytes(artifacts["original_stage2_ledger.csv"])
    overlay=rows_from_bytes(artifacts["eligibility_overlay.csv"])
    for row in overlay:
        for key in ("final_analysis_eligible","original_included_by_trial_policy","newly_included","newly_excluded","sequence_exclusion"):
            if key in row:row[key]=csv_bool(row[key])
    counts=application.analysis_counts(ledger,overlay,maps)
    by_trial={row["trial_key"]:row for row in maps}
    new_catches=[]
    for item in overlay:
        if item.get("newly_included"):
            timing_row=by_trial[item["trial_key"]]
            new_catches.append({"trial_key":item["trial_key"],"final_analysis_eligible":True,
                "assignment_basis":timing_row["assignment_basis"],"assignment_status":timing_row.get("assignment_status",""),
                "proposed_sample_raw":timing_row.get("proposed_sample_raw",""),
                "lag_estimate_ms":timing_row.get("lag_estimate_ms",""),
                "proposed_window_inside_recording":timing_row.get("proposed_window_inside_recording","")})
    assignments=Counter(row["assignment_basis"] for row in maps)
    eligible=Counter(row["assignment_basis"] for row in maps if csv_bool(row.get("final_analysis_eligible",False)))
    old_summary=json.loads(artifacts["application_summary.json"])
    final_summary=dict(old_summary,task1_staging_summary=old_summary,
        final_assignment_counts=dict(assignments),final_eligible_assignment_counts=dict(eligible),
        final_analysis_eligible=sum(csv_bool(row.get("final_analysis_eligible",False)) for row in maps),
        final_timing_version="s3-v002",pooled_onset_estimation_method=METHOD_ID,
        onset_medians_ms={code:summaries[str(code)]["median_ms"] for code in TRIGGERS})
    refreshed=dict(artifacts)
    refreshed.update({"analysis_counts.csv":csv_text(counts).encode(),
        "new_catch_timing_qc.csv":csv_text(new_catches,("trial_key","final_analysis_eligible","assignment_basis","assignment_status","proposed_sample_raw","lag_estimate_ms","proposed_window_inside_recording")).encode(),
        "task1_application_summary.json":artifacts["application_summary.json"],
        "application_summary.json":json.dumps(final_summary,indent=2,sort_keys=True,allow_nan=False).encode()+b"\n",
        "final_timing_summary.json":json.dumps(dict(subject_bids_id=f"sub-{subject}",assignment_counts=dict(assignments),
            final_eligible_assignment_counts=dict(eligible),n_final_analysis_eligible=final_summary["final_analysis_eligible"],
            pooled_onset_estimation_method=METHOD_ID,onset_medians_ms=final_summary["onset_medians_ms"]),
            indent=2,sort_keys=True,allow_nan=False).encode()+b"\n"})
    return refreshed,assignments,eligible,final_summary


def refresh_candidate(folder,report,artifacts,maps,subject,summaries,implementation_hashes,provenance_sha):
    """Repair unpublished candidate accounting and manifest after resumption."""
    outputs,assignments,eligible,summary=final_derived_tables(artifacts,maps,subject,summaries)
    outputs["code_visual_timing_task2.py.txt"]=Path(__file__).read_bytes()
    report.update(assignment_counts=dict(assignments),eligible_assignment_counts=dict(eligible),
        n_final_analysis_eligible=summary["final_analysis_eligible"],
        implementation_sha256=implementation_hashes,
        input_contract=dict(report.get("input_contract",{}),task2_implementation_sha256=implementation_hashes,
            donor_provenance_sha256=provenance_sha),
        task1_application_summary=summary["task1_staging_summary"])
    report["input_contract_sha256"]=digest(report["input_contract"])
    for name,value in outputs.items():
        if name=="manifest.json":continue
        path=folder/name
        path.write_bytes(value if isinstance(value,bytes) else value.encode())
        report.setdefault("output_sha256",{})[name]=file_hash(path)
    report["output_sha256"].pop("manifest.json",None)
    manifest_path=folder/"manifest.json"
    manifest_path.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    verify_bundle(folder)
    return report


def main():
    """Verify, build and atomically publish all 36 final participant bundles."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=Path("meg/params/pipeline_default.yaml"))
    parser.add_argument("--site-config",type=Path);parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--validate-only",action="store_true");args=parser.parse_args()
    cfg,*_=load_and_resolve_config(args.config,args.config.parent/"sub-04.yaml",site_config_path=args.site_config)
    root=Path(cfg["resolved_paths"]["derivatives_root"])/"visual_timing"
    source=root/"s3-v002-task1-staging";final=root/"s3-v002";build=root/"s3-v002-task2-build"
    handoff=json.loads((source/"task1_handoff.json").read_text())
    cached=load_pooled_cache(build,handoff) if build.exists() else None
    donors,summaries,sources=cached if cached is not None else pooled_donors(source,handoff)
    print(json.dumps({"donors":len(donors),"summaries":summaries},sort_keys=True),flush=True)
    if args.validate_only:return
    if final.exists():raise RuntimeError("Final s3-v002 already exists; refusing overwrite.")
    if build.exists() and cached is None:raise RuntimeError("Temporary build lacks a hash-valid donor cache; refusing overwrite.")
    build.mkdir(parents=True,exist_ok=True)
    donor_bytes=csv_text(donors).encode()
    donor_sha=hashlib.sha256(donor_bytes).hexdigest()
    provenance={"schema_version":"1.0","method":METHOD_ID,"estimator":"ordinary event-weighted median once per unique measured event",
        "trigger_summaries":summaries,"donor_excluded_event_keys":handoff["donor_excluded_event_keys"],
        "analysis_eligibility_filter":False,"n_donors":len(donors),"donor_sources":sources,
        "donor_csv_sha256":donor_sha,"transfer_uncertainty":"unknown;cohort spread does not quantify cross-context transfer"}
    provenance_bytes=json.dumps(provenance,indent=2,sort_keys=True,allow_nan=False).encode()+b"\n"
    provenance_sha=hashlib.sha256(provenance_bytes).hexdigest()
    implementation_hashes={name:file_hash(HERE/name) for name in ("visual_timing_task2.py","visual_timing_application.py",
        "visual_timing_utils.py","visual_timing_report.py","visual_timing.py","04_epoching.py")}
    # Per-participant rebuild consumes only sealed stage artifacts and metadata.
    for participant in handoff["participants"]:
        subject=participant["subject_bids_id"].replace("sub-","")
        subject_folder=build/f"sub-{subject}"/"ses-01"
        if subject_folder.exists():
            existing=verify_bundle(subject_folder)
            candidate_artifacts={name:(subject_folder/name).read_bytes() for name in
                ("original_stage2_ledger.csv","eligibility_overlay.csv","application_summary.json")}
            candidate_maps=read_csv(subject_folder/"trial_timing_map.csv")
            refresh_candidate(subject_folder,existing,candidate_artifacts,candidate_maps,subject,summaries,
                implementation_hashes,provenance_sha)
            print(f"reconciled candidate sub-{subject}",flush=True)
            continue
        staged_manifest,artifacts=application.load_task1_stage(source,subject)
        folder=source/f"sub-{subject}"/"ses-01"
        maps=rows_from_bytes(artifacts["trial_timing_map.csv"]);proposals=rows_from_bytes(artifacts["event_timing_proposals.csv"])
        estimates=json.loads(artifacts["estimates.json"]);maps,proposals,estimates=apply_pool(subject,maps,proposals,estimates,summaries,
            cfg["epoching"]["tmin_sec"],cfg["epoching"]["tmax_sec"])
        output={name:value for name,value in artifacts.items()}
        output.update({"trial_timing_map.csv":csv_text(maps).encode(),"event_timing_proposals.csv":csv_text(proposals).encode(),
            "estimates.json":json.dumps(estimates,indent=2,sort_keys=True,allow_nan=False).encode()+b"\n",
            "cohort_pooled_donor_measurements.csv":donor_bytes,"cohort_pooled_donor_provenance.json":provenance_bytes})
        report=dict(staged_manifest,version="s3-v002",mode="cohort_proposal",task1_staging=False,
            generated_utc=utc_now(),scientific_acceptance="pending_independent_review",
            onset_estimation_method={"method_id":METHOD_ID,"trigger_codes":list(TRIGGERS),"summary":summaries,
                "donor_provenance_sha256":hashlib.sha256(provenance_bytes).hexdigest(),"donor_measurements_sha256":donor_sha},
            implementation_sha256=implementation_hashes,
            input_contract=dict(staged_manifest.get("input_contract",{}),task1_manifest_sha256=file_hash(folder/"manifest.json"),
                donor_provenance_sha256=hashlib.sha256(provenance_bytes).hexdigest()),
            review_application_sha256=staged_manifest["review_application_sha256"])
        report["input_contract_sha256"]=digest(report["input_contract"])
        report["cohort_pooled_donor_provenance"]=provenance
        # Full explained reports are regenerated from the saved native windows;
        # those caches contain only previously captured diagnostic snippets.
        from visual_timing_report import build_report_artifacts,load_cached_waveforms
        contexts=json.loads(artifacts["contexts.json"])
        generated=build_report_artifacts(rows_from_bytes(artifacts["event_measurements.csv"]),maps,estimates,
            rows_from_bytes(artifacts["recording_summary.csv"]),rows_from_bytes(artifacts["physical_intervals.csv"]),
            contexts,report,load_cached_waveforms(folder),source_identity={"task1_manifest_sha256":file_hash(folder/"manifest.json"),
            "cohort_donor_provenance_sha256":hashlib.sha256(provenance_bytes).hexdigest()},proposals=proposals,
            donor_measurements=donors,evidence_bytes={"cohort_pooled_donor_measurements.csv":donor_bytes,
                "cohort_pooled_donor_provenance.json":provenance_bytes})
        output.update(generated)
        # Frozen input evidence bytes override report-normalized companions.
        for name in ("event_measurements.csv","saved_trial_metadata.csv","eligibility_overlay.csv","application.json",
                     "source_manifest.json","original_stage2_ledger.csv","review_notes_snapshot.md"):
            if name in artifacts:output[name]=artifacts[name]
        output,assignments,eligible,final_summary=final_derived_tables(output,maps,subject,summaries)
        output["code_visual_timing_task2.py.txt"]=Path(__file__).read_bytes()
        report.update(assignment_counts=dict(assignments),eligible_assignment_counts=dict(eligible),
            n_final_analysis_eligible=final_summary["final_analysis_eligible"],
            task1_application_summary=final_summary["task1_staging_summary"])
        report["input_contract"]["task2_implementation_sha256"]=implementation_hashes
        report["input_contract_sha256"]=digest(report["input_contract"])
        timing.publish(subject_folder,output,report,cfg["resolved_paths"]["data_root"])
        verified=verify_bundle(subject_folder)
        if verified["mode"]!="cohort_proposal" or verified.get("task1_staging") is not False:
            raise RuntimeError("Final participant publication mode check failed.")
        print(f"published candidate sub-{subject}",flush=True)
    # The cohort provenance is independently available at the version root.
    (build/"cohort_pooled_donor_measurements.csv").write_bytes(donor_bytes)
    (build/"cohort_pooled_donor_provenance.json").write_bytes(provenance_bytes)
    summary={"version":"s3-v002","participants":len(handoff["participants"]),"donor_summaries":summaries,
        "donor_measurements_sha256":donor_sha,"provenance_sha256":hashlib.sha256(provenance_bytes).hexdigest(),
        "task1_handoff_sha256":file_hash(source/"task1_handoff.json"),"scientific_acceptance":"pending_independent_review"}
    (build/"cohort_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    if not final.exists():os.replace(build,final)
    print(json.dumps(summary,sort_keys=True),flush=True)

if __name__=="__main__":main()
