"""Verify the three-participant epoch recovery without participant writes/jobs.

Usage: OPENBLAS_NUM_THREADS=1 python -m pytest meg/tests/test_epoch_recovery.py -q

Synthetic data verify exact pre-epoch annotation shortening, unchanged samples
and unrelated annotations, stale-approval rejection and actual MNE retention.
The real Stage-04 implementation verifies that disabling amplitude rejection
still rejects BAD intervals and trial-policy exclusions. Layered real selector
YAMLs verify participant-specific scope. Wrapper tests use disposable fake CLIs
to check that epoching and export receive the same recovery selector and label.
"""

from copy import deepcopy
from pathlib import Path
import os
import subprocess
import sys

import mne
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code/preprocessing_basic"))
from continuous_cleaning import annotation_rows,digest
from epoch_annotation_trims import apply_trims,plan_trims,validate_trial_bindings
from pipeline_config import read_yaml_mapping,deep_merge_dicts,apply_subject_profile
from test_inline_epoch_preparation import inputs
from test_epoch_export_smoke import EPOCHING


# -----------------------------------------------------------------------------
# Approved exact annotation changes and retained native samples.
# -----------------------------------------------------------------------------

def trim_fixture():
    """Build a pre-onset blink tail and a second protected annotation."""
    raw=mne.io.RawArray(np.zeros((1,6000)),mne.create_info(["MEG0111"],1000.,"mag"),first_samp=100000,verbose="ERROR")
    raw.set_annotations(mne.Annotations([1.3,5.],[.5,.3],["BAD_blink_auto","BAD_muscle"]))
    a=raw.annotations[0]
    rule={"source_fif":"sub21_run03_090926.fif","trial_key":"reviewed-trial", "behavioral_run":3,"behavioral_trial":30,
        "epoch_event_sample_raw":102000,"new_end_sample_raw":101499,
        "original_annotations_sha256":digest(annotation_rows(raw)),"annotation_description":str(a["description"]),
        "annotation_onset_sec":float(a["onset"]),"annotation_duration_sec":float(a["duration"])}
    cfg={"subject":{"id":"21"},"epoching":{"tmin_sec":-.5,"preparation_mode":"inline_no_ica","processing_run":"new",
        "annotation_end_trims":[rule]}}
    return raw,cfg,rule


def test_trim_changes_only_approved_duration_and_retains_exact_mne_epoch():
    raw,cfg,rule=trim_fixture();before=annotation_rows(raw);samples=raw.get_data().copy()
    def epochs():return mne.Epochs(raw,np.array([[102000,0,91]]),event_id={"deviant":91},tmin=-.5,tmax=2.,
        baseline=None,detrend=1,reject={"mag":4e-12},reject_by_annotation=True,preload=True,verbose="ERROR")
    assert len(epochs())==0
    plans=apply_trims(raw,cfg,rule["source_fif"])
    assert len(plans)==1 and len(epochs())==1
    assert np.isclose(raw.annotations.onset[0]+raw.annotations.duration[0],101.499,rtol=0,atol=1e-12)
    after=annotation_rows(raw);before["rows"][0][1]=after["rows"][0][1]
    assert after==before
    np.testing.assert_array_equal(raw.get_data(),samples)


def test_trim_refuses_stale_masks_wrong_boundaries_and_wrong_trial_keys():
    raw,cfg,rule=trim_fixture()
    wrong=deepcopy(cfg);wrong["epoching"]["annotation_end_trims"][0]["new_end_sample_raw"]+=1
    with pytest.raises(RuntimeError,match="one native sample"):plan_trims(raw,wrong,rule["source_fif"])
    raw.annotations.duration[1]+=.001
    with pytest.raises(RuntimeError,match="annotations changed"):apply_trims(raw,cfg,rule["source_fif"])
    with pytest.raises(RuntimeError,match="trial key"):validate_trial_bindings(cfg,[])


# -----------------------------------------------------------------------------
# Actual epoch-stage amplitude switch and participant-specific configuration.
# -----------------------------------------------------------------------------

def test_disabled_amplitude_keeps_bad_and_trial_policy_rejection(inputs):
    cfg,manifest,ann,calls=inputs
    cfg["epoching"]["amplitude_reject"]={"mag":1e-30}
    cfg["epoching"]["amplitude_rejection_enabled"]=False
    report=EPOCHING.run_epoching(cfg,source_manifest=manifest,overwrite=False)
    assert report["amplitude_rejection_enabled"] is False and report["amplitude_reject"]=={}
    assert report["counts"]["n_final_epochs"]==2
    assert report["counts"]["mne_drop_reasons"]=={"BAD_blink":1}
    assert report["counts"]["trial_policy_exclusions"]=={"gaze_break":1}


def test_recovery_selectors_disable_only_p23_amplitude_screen():
    root=Path(__file__).resolve().parents[1]/"params"
    base=read_yaml_mapping(root/"pipeline_default.yaml")
    for subject in (21,23,24):
        cfg=apply_subject_profile(deep_merge_dicts(base,read_yaml_mapping(root/f"epoch_recovery_20261117/sub-{subject}.yaml")))
        assert cfg["epoching"].get("amplitude_rejection_enabled",True)==(subject!=23)
        assert cfg["epoching"]["reject_by_annotation"] is True
        assert len(cfg["epoching"].get("annotation_end_trims",[]))==(2 if subject==21 else 0)


# -----------------------------------------------------------------------------
# Real shell wrapper dispatch, with disposable CLIs instead of processing/jobs.
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("subject",[21,23,24])
def test_epoch_and_export_wrapper_share_selected_recovery_config(tmp_path,subject):
    project=tmp_path/"project";derivatives=tmp_path/"derivatives";custom=project/"meg/params/recovery"
    required=[project/"meg/params/pipeline_default.yaml",project/"meg/params/pipeline_cluster.yaml",custom/f"sub-{subject}.yaml",
        derivatives/f"reports/sub-{subject}/ses-01/source_manifest.json",derivatives/f"reports/sub-{subject}/ses-01/03_tuning_choices.yaml",
        derivatives/f"visual_timing/s3-v002/reviews/sub-{subject}_review.json"]
    required += [project/"meg/code/preprocessing_basic"/n for n in ["04_epoching.py","06_export_neural_for_drsa.py","downstream_stage_utils.py","drsa_canonical_meg_channels.txt"]]
    for path in required:path.parent.mkdir(parents=True,exist_ok=True);path.write_text("fixture")
    binary=tmp_path/"bin";binary.mkdir();mamba=binary/"micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n");mamba.chmod(0o755)
    python=binary/"python";capture=tmp_path/"calls"
    python.write_text('#!/bin/bash\nprintf "%s\\t" "$@" >> "$DAD_TEST_CAPTURE"\nprintf "\\n" >> "$DAD_TEST_CAPTURE"\n');python.chmod(0o755)
    env=dict(os.environ,PATH=str(binary)+os.pathsep+os.environ["PATH"],DAD_PROJECT_ROOT=str(project),DAD_DERIVATIVES_ROOT=str(derivatives),
        MICROMAMBA_BIN=str(mamba),DAD_SUBJECT_CONFIG_DIR=str(custom),DAD_EPOCH_PREPARATION="inline_no_ica",
        DAD_PROCESSING_RUN="recovery",DAD_REBUILD="0",DAD_VALIDATE_ONLY="0",DAD_TEST_CAPTURE=str(capture),
        SLURM_ARRAY_TASK_ID=str(subject),SLURM_CPUS_PER_TASK="1",TMPDIR=str(tmp_path/"tmp"))
    wrapper=Path(__file__).resolve().parents[1]/"cluster/submit_epoch_export.sh"
    run=subprocess.run(["bash",str(wrapper)],env=env,capture_output=True,text=True)
    assert run.returncode==0,run.stdout+run.stderr
    calls=[line.split("\t") for line in capture.read_text().splitlines() if "epoching.py" in line or "export_neural" in line]
    assert len(calls)==3
    for call in calls:assert call[call.index("--subject-config")+1]==str(custom/f"sub-{subject}.yaml")
