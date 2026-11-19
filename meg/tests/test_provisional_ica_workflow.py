"""Check provisional ICA import boundaries without accessing DAD recordings.

Usage: python -m pytest meg/tests/test_provisional_ica_workflow.py -q

Temporary authorization tables exercise null versus explicit empty choices,
complete current-model index coverage, changed proposals, and the whole-cohort
validation barrier. No production tuning, native data or scheduler is touched.
"""
import copy
import json
from pathlib import Path
import sys
import os
import subprocess

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
import annotation_triage
import record_provisional_ica_choices as workflow
import audit_ica_application as application_audit
import mne
import numpy as np


def frozen_authorization(directory):
    """Make a complete working cohort with one unresolved participant."""
    rows = [{"subject": int(subject), "component_count": 2, "model_sha256": "model",
             "exclusions": None if int(subject) == 32 else [], "deferred_components_retained": [1]}
            for subject in workflow.WORKING_SUBJECT_IDS]
    auth = {"participants": rows, "working_cohort": [r["subject"] for r in rows],
        "authorized_participant_count": 34, "authorized_exclusion_count": 0,
        "blocked_participants": [32], "status": "provisional_component_processing_authorization",
        "user_visual_review_completed": False, "final_scientific_acceptance": False, "readiness_sha256": "sha"}
    directory.mkdir(parents=True)
    path = directory / "authorization.json"
    path.write_text(json.dumps(auth))
    path.with_name("choices.json").write_text(json.dumps({"authorization_file": path.name, "participants": rows}))
    return path, auth


def test_null_remains_pending_and_completed_empty_is_eligible(tmp_path):
    """An unresolved table cannot acquire processing permission from an empty list."""
    path, _ = frozen_authorization(tmp_path / "authorization")
    _, eligible = workflow.load_authorization(path)
    assert len(eligible) == 34 and all(row["subject"] != 32 for row in eligible)
    auth = json.loads(path.read_text())
    next(row for row in auth["participants"] if row["subject"] == 32)["exclusions"] = []
    path.write_text(json.dumps(auth))
    with pytest.raises(RuntimeError, match="disagree"):
        workflow.load_authorization(path)


@pytest.mark.parametrize("change", ["missing", "duplicate", "changed", "foreign", "unknown"])
def test_invalid_component_table_is_rejected(change):
    """Missing/duplicated indices or new exclusions never inherit authorization."""
    row = {"subject": 4, "model_sha256": "model", "component_count": 2,
        "exclusions": [], "deferred_components_retained": [1]}
    proposal = {"subject": "P04", "current_model_sha256": "model",
        "rows": [{"component": 0, "proposed": "keep"}, {"component": 1, "proposed": "defer"}]}
    workflow.validate_proposal(row, proposal)
    if change == "missing":
        proposal["rows"].pop()
    elif change == "duplicate":
        proposal["rows"][1] = copy.deepcopy(proposal["rows"][0])
    elif change == "changed":
        proposal["rows"][0]["proposed"] = "exclude"
    elif change == "foreign":
        proposal["current_model_sha256"] = "other"
    else:
        proposal["rows"][0]["proposed"] = "maybe"
    with pytest.raises(RuntimeError):
        workflow.validate_proposal(row, proposal)


def test_one_failed_participant_prevents_all_recording(tmp_path, monkeypatch):
    """The full table must join successfully before the first production mutation."""
    path, auth = frozen_authorization(tmp_path / "reports/ica_review_updates/auth")
    ready = tmp_path / "reports/ica_review_updates/final-ica-readiness-20261111-v1/summary.json"
    ready.parent.mkdir()
    ready.write_text(json.dumps({"participants": auth["participants"]}))
    monkeypatch.setattr(workflow, "file_hash", lambda path: "sha")
    monkeypatch.setattr(annotation_triage, "subject_config", lambda *args: {"resolved_paths": {"derivatives_root": str(tmp_path)}})
    calls = []
    def failed_validation(worker, tasks, workers):
        calls.append(worker)
        assert len(tasks) == 34
        return [{"subject": 4, "ok": False, "error": "changed source"}]
    monkeypatch.setattr(workflow, "map_participants", failed_validation)
    output = tmp_path / "reports/ica_review_updates/attempt"
    monkeypatch.setattr(sys, "argv", ["workflow", "--authorization", str(path), "--output", str(output), "--record"])
    with pytest.raises(RuntimeError, match="NO bulk decision writes"):
        workflow.main()
    assert calls == [workflow.validate_participant]
    assert not (output / "before").exists() and not (output / "recording.json").exists()


def test_task_preservation_uses_exact_native_trial_bounds():
    """All three signal views use identical accepted absolute sample windows."""
    raw = mne.io.RawArray(np.arange(40.).reshape(2, 20),
        mne.create_info(["MEG0111", "MEG0112"], 250., ["mag", "grad"]), first_samp=123, verbose=False)
    corrected = raw.copy(); corrected._data *= .5
    final = corrected.copy(); final._data += 1.
    rows = [{"condition_code": 71, "proposed_epoch_first_sample_raw": 126,
        "proposed_epoch_last_sample_raw": 130}]
    cfg = {"epoching": {"baseline_sec": None, "detrend": None}}
    original = raw.get_data().copy()
    sums, counts = application_audit.trial_signal_sums((raw, corrected, final), rows, cfg, [0, 1])
    assert counts == {"71": 1}
    np.testing.assert_array_equal(sums["condition_71_native"], original[:, 3:8])
    np.testing.assert_array_equal(sums["condition_71_ica_only"], original[:, 3:8] * .5)
    np.testing.assert_array_equal(sums["condition_71_final"], original[:, 3:8] * .5 + 1.)
    np.testing.assert_array_equal(raw.get_data(), original)
    rows[0]["proposed_epoch_first_sample_raw"] = 122
    with pytest.raises(RuntimeError, match="boundary"):
        application_audit.trial_signal_sums((raw, corrected, final), rows, cfg, [0, 1])


@pytest.mark.parametrize("fail_validation", [False, True])
@pytest.mark.parametrize("inline,validate_only", [(False, False), (True, False), (True, True)])
def test_epoch_wrapper_passes_accepted_timing_and_stops_before_export(tmp_path, fail_validation, inline, validate_only):
    """A missing/failed timing gate must never dispatch epoch creation or export."""
    project = tmp_path / "DAD"; derivatives = tmp_path / "derivatives"
    for name in ["meg/params/pipeline_default.yaml", "meg/params/pipeline_cluster.yaml", "meg/params/sub-04.yaml",
        "meg/code/preprocessing_basic/04_epoching.py", "meg/code/preprocessing_basic/06_export_neural_for_drsa.py",
        "meg/code/preprocessing_basic/downstream_stage_utils.py", "meg/code/preprocessing_basic/drsa_canonical_meg_channels.txt"]:
        path = project / name; path.parent.mkdir(parents=True, exist_ok=True); path.touch()
    for name in ["reports/sub-04/ses-01/source_manifest.json", "reports/sub-04/ses-01/03_tuning_choices.yaml",
                 "visual_timing/s3-v002/reviews/sub-04_review.json"]:
        path = derivatives / name; path.parent.mkdir(parents=True, exist_ok=True); path.touch()
    binary = tmp_path / "bin"; binary.mkdir()
    mamba = binary / "micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n"); mamba.chmod(0o755)
    python = binary / "python"
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$DAD_TEST_OUTPUT"\n'
        'if [[ "$DAD_TEST_FAIL_VALIDATION" == 1 && "$*" == *04_epoching.py* && "$*" == *--validate-only* ]]; then exit 7; fi\n')
    python.chmod(0o755)
    capture = tmp_path / "calls.txt"
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
        DAD_PROJECT_ROOT=str(project), DAD_DERIVATIVES_ROOT=str(derivatives), MICROMAMBA_BIN=str(mamba),
        SLURM_ARRAY_TASK_ID="4", SLURM_CPUS_PER_TASK="1", DAD_REBUILD="0",
        DAD_TEST_FAIL_VALIDATION=str(int(fail_validation)), DAD_TEST_OUTPUT=str(capture), TMPDIR=str(tmp_path / "tmp"))
    env.update(DAD_VALIDATE_ONLY=str(int(validate_only)))
    if inline:
        env.update(DAD_EPOCH_PREPARATION="inline_no_ica", DAD_PROCESSING_RUN="noica-fixture")
    else:
        env.pop("DAD_EPOCH_PREPARATION", None); env.pop("DAD_PROCESSING_RUN", None)
    wrapper = Path(__file__).resolve().parents[1] / "cluster/submit_epoch_export.sh"
    completed = subprocess.run(["bash", str(wrapper)], env=env, capture_output=True, text=True)
    calls = capture.read_text().splitlines()
    assert not any("run_pipeline.py" in line for line in calls)
    if fail_validation:
        assert completed.returncode == 7 and len(calls) == 2
        return
    assert completed.returncode == 0, completed.stderr
    if validate_only:
        assert len(calls) == 2 and not any("06_export" in line for line in calls)
        return
    assert len(calls) == 4
    if inline:
        assert "--source-manifest " + str(derivatives / "processing_runs/noica-fixture/reports/sub-04/ses-01/source_manifest.json") in calls[3]
        assert "--source-manifest " + str(derivatives / "reports/sub-04/ses-01/source_manifest.json") in calls[2]
    assert "--timing-review " + str(derivatives / "visual_timing/s3-v002/reviews/sub-04_review.json") in calls[1]
    assert "--timing-review" in calls[2] and "--resume" in calls[2]
    assert "06_export_neural_for_drsa.py" in calls[3] and "--resume" in calls[3]
    assert all("--tuning-config" in line for line in calls[1:])
