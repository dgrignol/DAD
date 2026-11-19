"""Verify P24 EDF reconstruction without touching participant data or jobs.

Usage: OPENBLAS_NUM_THREADS=1 python -m pytest \
    meg/tests/test_p24_auxiliary_reconstruction.py -q

Synthetic original-clock EDF fields exercise recording gaps, missing eyes,
unmatched intervals and the six-channel pupil/x/y order. Small native Raw
fixtures verify the final 325-channel ordering contract independently of array
generation, retained-sample equality, source immutability, run-specific high-pass
metadata and unchanged annotation rejection. Array tampering must fail closed.
"""

from copy import deepcopy
import json
from pathlib import Path
import sys

import mne
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code/preprocessing_basic"))
import p24_auxiliary_reconstruction as aux
from continuous_cleaning import annotation_rows, file_hash


# -----------------------------------------------------------------------------
# Original clock mapping: missing data must never become a connected waveform.
# -----------------------------------------------------------------------------

def test_mapping_preserves_gaps_missing_eyes_and_measured_zero():
    """Zero fill stays distinct from measured zero in the validity mask."""
    samples = np.array([[t, 10, 0, 20, 30, 40, 50] for t in [1., 1.001, 1.003, 2., 2.001]])
    samples[1, 4:] = np.nan
    segment = {"status":"matched", "start":1., "stop":1.004,
        "edf_anchor_sec":1., "meg_anchor_sample":1100, "sfreq":1000.}
    data, valid = aux.map_samples(samples, [segment], first_samp=1000, n_times=1200, sfreq=1000.)
    assert data.shape == (6,1200)
    assert valid[1,100] and data[1,100] == 0  # measured zero gaze coordinate
    assert not valid[:,102].any() and not data[:,102].any()  # real EDF sample gap
    assert valid[:3,101].all() and not valid[3:,101].any()  # only right eye lost
    assert not valid[:,1000:].any()  # unreviewed second EDF block is not filled
    assert np.array_equal(data[:,100], [10,0,20,30,40,50])
    with pytest.raises(RuntimeError, match="overlap"):
        aux.map_samples(samples, [segment,segment], first_samp=1000,n_times=1200,sfreq=1000.)


# -----------------------------------------------------------------------------
# In-memory replacement: all retained native samples and acquisition info agree.
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("run", [1,5])
def test_apply_retains_meg_triggers_annotations_and_run_highpass(tmp_path, run):
    """Replace only authorized auxiliaries, with consistent data across runs."""
    native_names = ["MEG0111", "STI101"] + (aux.CHANNELS+["MISC008"] if run==5 else [])
    info = mne.create_info(native_names,1000.,["mag","stim"]+(["misc"]*7 if run==5 else []))
    with info._unlock(): info["highpass"] = .1 if run==1 else 0.
    values = np.arange(len(native_names)*200,dtype=float).reshape(len(native_names),200)
    raw = mne.io.RawArray(values,info,first_samp=30000,verbose="ERROR")
    raw.set_annotations(mne.Annotations([.05],[.02],["BAD_blink"]))
    before = raw.copy()
    data=np.full((6,200),23,dtype=np.float32);valid=np.ones(data.shape,dtype=bool)
    data[:,20:30]=0;valid[:,20:30]=False
    path=tmp_path/"arrays.npz";np.savez_compressed(path,data=data,valid=valid)
    source=f"sub24_run{run:02d}_100926.fif"
    plan={"final_channel_order":["MEG0111",*aux.CHANNELS,"STI101"],
        "reference_acquired_auxiliary_source":"sub24_run05_100926.fif",
        "reconstructed_runs":[{"source_filename":source,"array":path.name,
            "array_sha256":file_hash(path),"zero_fill_samples_by_channel":dict.fromkeys(aux.CHANNELS,10)}]}
    manifest=tmp_path/"manifest.json";manifest.write_text(json.dumps(plan))
    cfg={"epoching":{"auxiliary_reconstruction_manifest":"manifest.json"},
        "resolved_paths":{"derivatives_root":str(tmp_path)}}
    evidence={"definition":plan,"manifest_sha256":file_hash(manifest)}
    qc=aux.apply_plan(raw,cfg,source,evidence)
    assert raw.ch_names == plan["final_channel_order"]
    assert qc["retained_samples_verified"]
    np.testing.assert_array_equal(raw.get_data(picks=["MEG0111","STI101"]),before.get_data(picks=["MEG0111","STI101"]))
    np.testing.assert_array_equal(raw.get_data(picks=aux.CHANNELS),data)
    assert annotation_rows(raw)==annotation_rows(before)
    assert raw.info["highpass"]==before.info["highpass"]
    assert "MISC008" not in raw.ch_names
    assert before.ch_names==native_names  # the source fixture remains intact
    if run==5:np.testing.assert_array_equal(before.get_data(picks=aux.CHANNELS),values[2:8])
    path.write_bytes(b"tampered")
    with pytest.raises(RuntimeError,match="array changed"):
        aux.apply_plan(before,cfg,source,evidence)


def test_old_variants_are_unaffected_and_other_subjects_fail(tmp_path):
    """Optional adapter has no effect unless an explicit P24 plan is selected."""
    cfg={"epoching":{},"subject":{"id":"21"},"resolved_paths":{"derivatives_root":str(tmp_path)}}
    assert aux.validate_plan(cfg,[]) is None
    cfg["epoching"]={"auxiliary_reconstruction_manifest":"manifest.json","preparation_mode":"inline_no_ica","processing_run":"new"}
    with pytest.raises(RuntimeError,match="P24"):
        aux.validate_plan(cfg,[])
    cfg["epoching"]["auxiliary_reconstruction_manifest"]="../raw/manifest.json"
    with pytest.raises(RuntimeError,match="relative"):
        aux.plan_path(cfg)


def test_five_run_325_channel_epoch_concatenation_preserves_acquisition_metadata(tmp_path):
    """Regress the actual fifth-run nchan failure using all 306 MEG sensors."""
    import pandas as pd
    channels_file=Path(__file__).resolve().parents[1]/"code/preprocessing_basic/drsa_canonical_meg_channels.txt"
    meg=[n.strip() for n in channels_file.read_text().splitlines() if n.strip() and not n.startswith("#")]
    common=[*meg,*[f"IAS{i:03d}" for i in range(11)],"STI101","SYS201"]
    final=[*meg,*[f"IAS{i:03d}" for i in range(11)],*aux.CHANNELS,"STI101","SYS201"]
    assert len(common)==319 and len(final)==325
    rows=[]
    for run in range(1,6):
        target=tmp_path/f"run{run}.npz";np.savez_compressed(target,data=np.full((6,400),run,dtype=np.float32),valid=np.ones((6,400),bool))
        rows.append({"source_filename":f"sub24_run{run:02d}_100926.fif","array":target.name,
            "array_sha256":file_hash(target),"zero_fill_samples_by_channel":dict.fromkeys(aux.CHANNELS,0)})
    plan={"final_channel_order":final,"reference_acquired_auxiliary_source":"sub24_run05_100926.fif","reconstructed_runs":rows}
    manifest=tmp_path/"manifest.json";manifest.write_text(json.dumps(plan))
    cfg={"epoching":{"auxiliary_reconstruction_manifest":"manifest.json"},"resolved_paths":{"derivatives_root":str(tmp_path)}}
    evidence={"definition":plan,"manifest_sha256":file_hash(manifest)}
    epochs=[]
    for run in range(1,6):
        names=common+(aux.CHANNELS+["MISC008"] if run==5 else [])
        kinds=[("mag" if n.endswith("1") else "grad") if n.startswith("MEG") else ("stim" if n=="STI101" else "misc") for n in names]
        info=mne.create_info(names,1000.,kinds)
        with info._unlock():info["highpass"] = .1 if run<5 else 0.
        data=np.random.default_rng(run).normal(size=(len(names),400))*1e-15
        raw=mne.io.RawArray(data,info,first_samp=run*10000,verbose="ERROR")
        raw.set_annotations(mne.Annotations([.05],[.02],["BAD_blink"]))
        aux.apply_plan(raw,cfg,rows[run-1]["source_filename"],evidence)
        epochs.append(mne.Epochs(raw,np.array([[raw.first_samp+200,0,91]]),event_id={"deviant":91},
            tmin=0.,tmax=.05,baseline=None,detrend=1,reject_by_annotation=True,
            reject={"mag":4e-12,"grad":5e-10},picks="all",preload=True,
            metadata=pd.DataFrame({"run":[run],"online_highpass_hz":[.1 if run<5 else 0.]}),verbose="ERROR"))
    joined=mne.concatenate_epochs(epochs,verbose="ERROR")
    assert joined.info["nchan"]==325 and len(joined)==5
    assert joined.metadata["online_highpass_hz"].tolist()==[.1,.1,.1,.1,0.]
    assert "MISC008" not in joined.ch_names
    for run in range(1,6):np.testing.assert_array_equal(joined[run-1].get_data(picks=aux.CHANNELS),epochs[run-1].get_data(picks=aux.CHANNELS))
