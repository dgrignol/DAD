#!/usr/bin/env python3
"""Bounded historical regression and independent cached participant-pilot checks.

Usage from the DAD root after the four bounded participant pilots exist::

    python meg/code/preprocessing_basic/visual_timing_validation.py \
      --project-root "$PWD" \
      --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
      --raw-root /Volumes/MORWUR/Projects/DAMIANO/Attention/data \
      --pilot-version s3-pilot-20261020-v001 --version s3-validation-20261020-v001

Reads 14 preselected 250-ms historical test windows: v22's three run starts plus
control/cadence, v23's late-arm event plus control/cadence, v24's three low-tail
onsets plus control/cadence/mid-trial. It does not rerun an historical audit.
Original v21/100-Hz recordings are unavailable at recorded paths; no substitute
calibration is invented. Raw FIF hashes are not recomputed by reading whole
recordings; recorded provenance hashes and current header/size/anchor evidence
are kept distinct. Participant comparison uses cached pilot windows only.

An independent maximum-positive-slope sample is a coarse waveform cross-check,
not a second accepted detector. Reports expose differences and flags without
automatically accepting the method, cohort estimates or participant inclusion.
All outputs are new guarded derivatives; original evidence is never modified.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from contextlib import nullcontext
from types import SimpleNamespace

import mne
import numpy as np

from trial_preparation import csv_text, file_hash, read_csv
from preprocessing_stage_utils import utc_now
from visual_timing import DEFAULT_METHOD_PATH, diagnostic_images, load_method, publish
from visual_timing_utils import CONDITION_CODES, measure_edge, verify_bundle
from downstream_stage_utils import csv_bool


# -----------------------------------------------------------------------------
# Prespecified historical windows and independent coarse edge comparison
# -----------------------------------------------------------------------------


def historical_cases(root):
    """Resolve exact catalog inputs and choose cases from existing measurements."""
    root=Path(root)
    definitions=[
        ("v22",root/"tmp/260731_sub98/sub98_testPhotodiodeTimingTriggers_260731.fif",
         root/"tmp/photodiode_timing_audit/v22_scheduled_onset_20260731/photodiode_events_detailed.csv"),
        ("v23",root/"tmp/260820_sub97/Sub97_testPhotodiodeTimingTriggers_260820.fif",
         root/"tmp/260820_sub97/analysis/photodiode_events_detailed.csv"),
        ("v24",root/"tmp/260820_sub96/Sub96_testPhotodiodeTimingTriggers_260820.fif",
         root/"tmp/photodiode_timing_audit/v24_all_trial_guard_resume_retest_20260820/photodiode_events_detailed.csv")]
    result=[]
    for label,fif,table in definitions:
        rows=read_csv(table)
        previous=None
        for row in sorted(rows,key=lambda r:int(r["sample"])):
            row["_previous_code"] = int(previous["trigger"]) if previous else None
            row["_previous_gap_samples"] = int(row["sample"])-int(previous["sample"]) if previous else None
            previous=row
        onsets=[r for r in rows if int(r["trigger"]) in CONDITION_CODES]
        if label=="v22":
            selected=[next(r for r in onsets if r["run"]==run) for run in sorted({r["run"] for r in onsets})]
        elif label=="v23":selected=[min(onsets,key=lambda r:float(r["delay_ms"]))]
        else:selected=[r for r in onsets if float(r["delay_ms"])<18]
        selected.append(sorted(onsets,key=lambda r:float(r["delay_ms"]))[len(onsets)//2])
        selected.append([r for r in rows if r["trigger"]=="153"][1])
        if label=="v24":selected.append(next(r for r in rows if r["trigger"]=="111"))
        result.append((label,fif,table,selected))
    if sum(len(x[3]) for x in result)!=14:
        raise RuntimeError("Historical case selection no longer matches the 14-window bound.")
    return result


def independent_slope(values,anchor,sfreq,polarity=1,*,warmup=False):
    """Locate the steepest positive sample increment without fractional logic."""
    a=max(1,anchor+round((0 if warmup else -.020)*sfreq))
    b=min(len(values),anchor+round(.040*sfreq))
    gradient=np.diff(values*polarity)
    index=int(a+np.argmax(gradient[a-1:b-1]))
    return (index-anchor)*1000/sfreq


# -----------------------------------------------------------------------------
# Native historical reads and participant cache-only comparisons
# -----------------------------------------------------------------------------


def run_validation(project_root,derivatives_root,raw_root,pilot_version,version,*,historical_cache_from=None):
    """Publish bounded validation evidence; no participant FIFs are reopened."""
    root=Path(project_root)
    output=Path(derivatives_root)/"visual_timing"/version/"checks"
    if output.exists():raise FileExistsError(output)
    method=load_method(DEFAULT_METHOD_PATH)
    implementation={p.name:file_hash(p) for p in [Path(__file__),Path(__file__).with_name("visual_timing.py"),Path(__file__).with_name("visual_timing_utils.py")]}
    historical_cache, cached_provenance, cache_reference = None, {}, None
    if historical_cache_from:
        folder=Path(derivatives_root)/"visual_timing"/historical_cache_from/"checks"
        manifest=json.loads((folder/"manifest.json").read_text())
        if manifest.get("status")!="completed" or manifest.get("mode")!="bounded_validation":
            raise RuntimeError("Historical cache is not a completed bounded-validation bundle.")
        for name,expected in manifest["output_sha256"].items():
            if Path(name).name!=name or file_hash(folder/name)!=expected:
                raise RuntimeError("Historical cache output fingerprint mismatch.")
        historical_cache={}
        with np.load(folder/"historical_all_windows.npz",allow_pickle=False) as cache:
            for entry in read_csv(folder/"historical_all_windows_index.csv"):
                historical_cache[entry["event_key"]]=cache[entry["array"]].copy()
        if len(historical_cache)!=14:raise RuntimeError("Historical cache must contain exactly 14 native windows.")
        cached_provenance={Path(r["source"]).name:r for r in json.loads((folder/"source_provenance.json").read_text())}
        cache_reference={"version":historical_cache_from,"manifest_sha256":file_hash(folder/"manifest.json")}
    measurements,waveforms,comparisons,provenance=[],{},[],[]
    for label,fif,table,selected in historical_cases(root):
        if historical_cache is None:
            reader=mne.io.read_raw_fif(fif,preload=False,verbose="ERROR")
        else:
            previous=cached_provenance[fif.name]
            if file_hash(table)!=previous["archived_event_table_sha256"]:
                raise RuntimeError("Archived event table changed after historical-window caching.")
            time_axis=next(v[0] for k,v in historical_cache.items() if k.startswith(label+":"))
            reader=nullcontext(SimpleNamespace(info={"sfreq":1000/(time_axis[1]-time_axis[0])},
                first_samp=previous["first_samp"],n_times=previous["n_times"],filenames=previous["split_files"]))
        with reader as raw:
            sf=float(raw.info["sfreq"])
            inventory=table.parent/"source_inventory.csv"
            if inventory.is_file():
                matches=[r for r in read_csv(inventory) if Path(r["source_path"]).name==fif.name]
                if len(matches)!=1 or int(matches[0]["size_bytes"])!=fif.stat().st_size:
                    raise RuntimeError("Historical raw source identity/size differs from recorded provenance.")
            for i,original in enumerate(selected):
                sample=int(original["sample"])
                center=sample-int(raw.first_samp)
                start,stop=center-round(.075*sf),center+round(.175*sf)
                if start<0 or stop>raw.n_times:raise RuntimeError("Historical window crosses recording boundary.")
                anchor=center-start
                code=int(original["trigger"])
                key=f"{label}:{sample}:{code}"
                if historical_cache is None:
                    values=raw.get_data(picks=["STI101","MISC008"],start=start,stop=stop)
                else:
                    cached=historical_cache[key]
                    values=cached[1:]
                    if values.shape!=(2,stop-start) or not np.isclose(cached[0,anchor],0):
                        raise RuntimeError("Historical cached native coordinates differ.")
                if int(round(values[0,anchor]))!=code:
                    raise RuntimeError("Historical recorded anchor/code differs from archived measurement table.")
                post_warmup=bool(code in CONDITION_CODES and original["_previous_code"]==153
                    and 0<original["_previous_gap_samples"]*1000/sf<=method["post_warmup_max_gap_ms"])
                measured=measure_edge(values[1],sf,anchor,1,int(original["expected_pulse_frames"]),120,method,
                                      warmup=code==153,post_warmup=post_warmup)
                measured.update(event_key=key,source_fif=fif.name,digital_sample_raw=sample,
                                first_samp=int(raw.first_samp),sfreq_hz=sf,trigger_code=code,
                                archived_lag_ms=float(original["delay_ms"]),historical_source=label)
                measured["difference_from_archived_ms"]=(measured["lag_ms"]-float(original["delay_ms"])) if measured["lag_ms"] is not None else None
                measured["independent_slope_lag_ms"]=independent_slope(values[1],anchor,sf,warmup=code==153 or post_warmup)
                measurements.append(measured);waveforms[key]=(values,anchor,sf)
            provenance.append({"source":str(fif),"size_bytes":fif.stat().st_size,
                               "first_samp":int(raw.first_samp),"n_times":int(raw.n_times),
                               "split_files":[str(p) for p in raw.filenames],
                               "archived_event_table":str(table),"archived_event_table_sha256":file_hash(table),
                               "raw_sha256_recomputed":False})
    for subject in ("08","09","11","24"):
        folder=Path(derivatives_root)/"visual_timing"/pilot_version/f"sub-{subject}/ses-01"
        manifest=verify_bundle(folder)
        if manifest["mode"]!="pilot" or manifest["n_signal_windows_read"]>12:
            raise RuntimeError("Expected an intact bounded pilot bundle.")
        rows={r["event_key"]:r for r in read_csv(folder/"event_measurements.csv")}
        with np.load(folder/"waveform_windows.npz",allow_pickle=False) as cache:
            for entry in read_csv(folder/"diagnostic_index.csv"):
                t,sti,diode=cache[entry["array"]]
                row=rows[entry["event_key"]]
                anchor=int(np.flatnonzero(np.isclose(t,0))[0])
                sf=float(row["sfreq_hz"])
                polarity=int(row.get("polarity") or 1)
                primary=float(row["lag_ms"]) if row.get("lag_ms") else None
                independent=independent_slope(diode,anchor,sf,polarity,
                    warmup=row["trigger_code"]=="153" or csv_bool(row.get("post_warmup_context",False)))
                comparisons.append({"event_key":entry["event_key"],"measurement_status":row["measurement_status"],
                                    "candidate_50_lag_ms":primary,"independent_slope_lag_ms":independent,
                                    "difference_ms":primary-independent if primary is not None else None,
                                    "digital_anchor_value":int(round(sti[anchor])),"trigger_code":row["trigger_code"]})
    lines=["# Stage 3 bounded validation","","Scientific acceptance remains pending.","",
           f"Historical windows assessed: 14 x 250 ms; new historical signal reads: {0 if historical_cache is not None else 14}. Participant comparison: cached pilot windows only.","",
           "The maximum-slope comparison is independent of fractional-amplitude extraction, but is a coarse check rather than an accepted timing definition.","",
           "Inspect archived differences, candidate flags and raw waveforms. Early/low-tail observations remain visible even when the baseline/pulse QC flags them.","",
           "The original v21 FIF/MAT/CSVs and previous 100-Hz FIF are unavailable at recorded locations; those historical results were not independently remeasured.","",
           "![Historical native waveform comparisons](diagnostics.png)",""]
    contents={"historical_measurements.csv":csv_text(measurements),"participant_cache_checks.csv":csv_text(comparisons),
              "source_provenance.json":json.dumps(provenance,indent=2,sort_keys=True),
              "report.md":"\n".join(lines),**diagnostic_images(measurements,waveforms)}
    # Store all historical reads, not only the plotted representative subset.
    import io
    cache=io.BytesIO();np.savez_compressed(cache,**{f"window_{i:03d}":np.vstack([(np.arange(v.shape[1])-a)*1000/sf,v])
        for i,(v,a,sf) in enumerate(waveforms.values())})
    contents["historical_all_windows.npz"]=cache.getvalue()
    contents["historical_all_windows_index.csv"]=csv_text([{"array":f"window_{i:03d}","event_key":key}
        for i,key in enumerate(waveforms)])
    report={"schema_version":"1.0","status":"completed","mode":"bounded_validation",
            "scientific_acceptance":"pending_independent_review","version":version,"pilot_version":pilot_version,
            "generated_utc":utc_now(),"historical_windows_read":0 if historical_cache is not None else 14,
            "historical_windows_assessed":14,"historical_cache":cache_reference,"participant_windows_reopened":0,
            "participant_cached_windows_checked":len(comparisons),
            "implementation_sha256":implementation}
    for name,expected in implementation.items():
        p=Path(__file__).with_name(name)
        if file_hash(p)!=expected:raise RuntimeError("Implementation changed during validation.")
        contents[f"code_{name}.txt"]=p.read_bytes()
    publish(output,contents,report,raw_root)
    return report


def main():
    """Run the fixed validation scope and write a new guarded evidence bundle."""
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("project-root","derivatives-root","raw-root","pilot-version","version"):
        p.add_argument(f"--{key}",required=True)
    p.add_argument("--historical-cache-from",help="Replay an intact 14-window historical cache without raw signal rereads.")
    a=p.parse_args()
    # Use the same non-traversing identifier rule as production versions.
    import re
    for version in (a.version,a.pilot_version,a.historical_cache_from):
        if version is None:continue
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}",version) or ".." in version:
            raise ValueError("Invalid immutable version identifier.")
    result=run_validation(a.project_root,a.derivatives_root,a.raw_root,a.pilot_version,a.version,
                          historical_cache_from=a.historical_cache_from)
    print(json.dumps({k:result[k] for k in ("status","version","historical_windows_read","historical_windows_assessed","participant_cached_windows_checked")},sort_keys=True))


if __name__=="__main__":
    main()
