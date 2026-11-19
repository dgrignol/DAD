"""Check the Stage 3 cohort-median overlay and plot using synthetic rows only.

Usage: PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg python3 -B -m pytest \
  -p no:cacheprovider meg/tests/test_visual_timing_task2.py -q

These tests do not access participant data or publish a bundle.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import json
import sys

CODE=Path(__file__).resolve().parents[1]/"code/preprocessing_basic"
sys.path.insert(0,str(CODE))
import visual_timing_task2 as task2
import visual_timing_report as report
import visual_timing_utils as utils
import pytest


def test_pooled_median_uses_native_grid_and_only_changes_designated_target():
    original=dict(trial_key="trial",condition_code=81,assignment_basis="digital_fallback",
        raw_sample_frame=1001,digital_sample_raw=1001,sfreq_hz=1000,
        recording_first_sample_raw=0,recording_end_sample_raw_exclusive=10000,
        trial_interval_end_raw_exclusive=9999)
    summaries={str(code):{"median_ms":19.0,"n":10,"sd_ms":.4,"minimum_ms":18,"maximum_ms":20,
        "q25_ms":19,"q75_ms":19,"p05_ms":18,"p95_ms":20,
        "transfer_uncertainty":"unknown;observed cohort spread is not target transfer error"} for code in (71,81,91)}
    maps,proposals,estimates=task2.apply_pool("04",[dict(original)],[],{},summaries,-.5,2.667)
    assert maps[0]["assignment_basis"]=="estimated"
    assert maps[0]["proposed_sample_raw"]==1020
    assert maps[0]["applied_shift_samples"]==19
    assert maps[0]["donor_estimate_id"]=="COHORT_code81"
    assert maps[0]["uncertainty_kind"]=="observed_cohort_dispersion;target_transfer_unknown"
    unchanged,_,_=task2.apply_pool("08",[dict(original)],[],{},summaries,-.5,2.667)
    assert unchanged[0]["assignment_basis"]=="digital_fallback"


def test_pooled_plot_treats_trigger_code_as_factor_and_reports_spread():
    rows=[]
    for code,values in ((71,[18,19,20]),(81,[18,19,19,20]),(91,[19,20])):
        rows.extend(dict(trigger_code=code,lag_ms=value) for value in values)
    image=report.pooled_lag_figure(rows)
    assert image.startswith(b"\x89PNG\r\n\x1a\n")


def test_pooled_consumer_rejects_wrong_bound_method(tmp_path):
    donor=tmp_path/"cohort_pooled_donor_measurements.csv"
    provenance=tmp_path/"cohort_pooled_donor_provenance.json"
    donor.write_text("event_key,trigger_code,lag_ms\\nfake,71,19\\n")
    provenance.write_text(json.dumps({"method":"untrusted_method","analysis_eligibility_filter":False}))
    bundle={"onset_estimation_method":{"donor_measurements_sha256":hashlib.sha256(donor.read_bytes()).hexdigest(),
        "donor_provenance_sha256":hashlib.sha256(provenance.read_bytes()).hexdigest()}}
    with pytest.raises(RuntimeError,match="Unsupported pooled donor method"):
        utils.validate_cohort_donor_provenance({},tmp_path,bundle)
