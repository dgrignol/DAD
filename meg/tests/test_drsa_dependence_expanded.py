"""Check expanded diagnostic provenance, PE panels and full-cohort publication.

Usage: python -m pytest meg/tests/test_drsa_dependence_expanded.py -q
Synthetic immutable bundles stand in for Stage-06/07 data; no MATLAB, cluster
submission, or production writes occur. Real MATLAB checks are separate in
validate_drsa_dependence_models.m.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.io import savemat

CODE = Path(__file__).resolve().parents[1] / "code/preprocessing_basic"
sys.path.insert(0, str(CODE))
from continuous_cleaning import file_hash
from drsa_dependence import DEVIANT, LEGACY_SOURCES, model_axes, read_condition
from drsa_dependence_inputs import brain_result_path, validate_brain_inputs
from plot_drsa_dependence import run_group


@pytest.fixture
def expanded(tmp_path):
    """Two aligned participants with legacy and scalar/multifeature PE sources."""
    sources = LEGACY_SOURCES + ["pe_latent_curvature", "pe_prior_unconstrained_50", "pe_belief_update"]
    groups = ["frontal"]
    source_run = tmp_path / "processing_runs/source"
    brain = source_run / "expansion"
    root = brain / "diagnostics"
    axis = np.arange(70) / 25
    timing = dict(first_complete_occlusion_sec=0.5, reappearance_sec=1.5)
    for subject in ("sub-04", "sub-05"):
        handoff = dict(subject_bids_id=subject, session_bids_id="ses-01",
                       identity_values=[1, 2, 3, 4], time_axis_sec=axis.tolist(),
                       channel_labels=["MEG0111", "MEG0112", "MEG0113"], sampling_rate_hz=25,
                       occlusion_timing=timing, trajectory_models={"observed": "/mac/data/trajectory.mat"},
                       source_manifest=str(source_run / f"reports/{subject}/ses-01/source_manifest.json"),
                       export_report_json=str(source_run / f"reports/{subject}/ses-01/06_export.json"),
                       condition_outputs=[dict(name=DEVIANT, neural_mat=str(source_run / f"analysis/{subject}/neural.mat"),
                                               sidecar_json=str(source_run / f"analysis/{subject}/neural.json"))])
        directory = root / subject / "ses-01"
        directory.mkdir(parents=True)
        (directory / "handoff.json").write_text(json.dumps(handoff))
        axes = model_axes(DEVIANT, sources)
        matrices = {key: np.eye(len(axis)) for key in axes}
        metadata = {}
        for src in sources:
            feature = {"pe_latent_curvature": "latent_curvature", "pe_belief_update": "belief_update"}.get(src, "position")
            count = {"pe_latent_curvature": 1, "pe_prior_unconstrained_50": 3}.get(src, 2)
            metadata[src] = dict(feature=feature, feature_count=count, rdm_distance="euclidean")
            if src not in LEGACY_SOURCES:
                for key in (f"{src}_self", f"observed_vs_{src}", f"predicted_vs_{src}"):
                    matrices[key] = np.full((len(axis), len(axis)), np.nan)
        axes["neural_frontal"] = ("MEG (frontal)", "MEG (frontal)")
        matrices["neural_frontal"] = np.eye(len(axis))
        savemat(directory / f"{DEVIANT}.mat", dict(condition=DEVIANT, matrices=matrices,
                time_axis_sec=axis, identities=[1,2,3,4], shared_prefix_frames=3,
                selected_channel_labels={"frontal": handoff["channel_labels"]}, model_metadata=metadata),long_field_names=True)
        report = dict(status="completed", conditions=[DEVIANT], model_sources=sources, sensor_groups=groups,
                      smooth_neural=True, smooth_window_samples=15, contract_sha256=subject,
                      inputs=dict(processing_run="source",source_export_contract_sha256="source"),
                      matrix_axes={DEVIANT:{key:dict(row=labels[0],column=labels[1]) for key,labels in axes.items()}},
                      sha256={p.name:file_hash(p) for p in directory.iterdir()})
        (directory / "report.json").write_text(json.dumps(report))
        reports = brain / "reports" / subject / "ses-01"
        reports.mkdir(parents=True)
        (reports / "07_matlab_drsa_handoff.json").write_text(json.dumps(handoff))
        records, inventory = [], []
        for src in sources:
            path = brain_result_path(brain / "analysis",subject,"ses-01",DEVIANT,src,"frontal")
            path.parent.mkdir(parents=True)
            savemat(path, dict(dRSA_neural_position=np.full((len(axis),len(axis)),np.nan),
                    handoff=handoff, time_axis_sec=axis, condition_str=DEVIANT, model_source=src,
                    model_feature=metadata[src]["feature"], sensor_group_tag="frontal",
                    selected_channel_labels=handoff["channel_labels"],
                    params=dict(dRSAtype="corr",neuralDistMeasure="correlation",modelDistMeasure="euclidean")),long_field_names=True)
            inventory.append(str(path))
            records.append(dict(analysis=f"alt:{src}",status="completed",return_code=0,
                    command=["/matlab","-batch",f"doSmoothNeural=true;smoothWindowSamples=15;path_type='{src}';dRSAtypeToRun='corr';"]))
        (reports / "07_matlab_drsa_report.json").write_text(json.dumps(dict(
            status="completed",subject_bids_id=subject,session_bids_id="ses-01",selected_analyses=["alt"],
            model_sources=sources,matlab_records=records,outputs=dict(mat_results=inventory))))
    args = argparse.Namespace(overwrite=False, resume=True, validate_only=False, lag_window_sec=2.0,
                              model_sources=sources, brain_results_root=str(brain))
    cfg = dict(layout=dict(session_bids_id="ses-01"),epoching=dict(processing_run="source"))
    paths = dict(analysis_dir=source_run / "analysis/sub-04/ses-01")
    return root, brain, cfg, args, paths, sources


def test_pe_group_keeps_axes_nan_support_and_legacy_overview(expanded):
    root, brain, cfg, args, paths, sources = expanded
    result = run_group(cfg,args,paths,root,["sub-04","sub-05"],[DEVIANT],["frontal"])
    output = Path(result["directory"])
    assert result["status"] == "completed"
    assert len(result["panel_axes"]) == 4
    pe = np.load(output / f"{DEVIANT}_frontal_pe_latent_curvature.npz")
    assert not pe["pe_latent_curvature_self_matrix_n"].any()
    assert np.isnan(pe["pe_latent_curvature_self_curves_sem"]).all()
    assert result["panel_axes"][f"{DEVIANT}_frontal_pe_latent_curvature"]["observed_vs_pe_latent_curvature"] == dict(row="Observed position",column="Latent curvature PE")
    assert (output / f"{DEVIANT}_frontal.npz").exists()
    assert len(result["brain_provenance"]) == 2
    assert run_group(cfg,args,paths,root,["sub-04","sub-05"],[DEVIANT],["frontal"]) == result


def test_missing_cohort_or_wrong_axes_cannot_publish(expanded):
    root, brain, cfg, args, paths, sources = expanded
    with pytest.raises(FileNotFoundError):
        run_group(cfg,args,paths,root,["sub-04","sub-05","sub-06"],[DEVIANT],["frontal"])
    assert not (root / "group").exists()
    path = brain_result_path(brain / "analysis","sub-05","ses-01",DEVIANT,"pe_latent_curvature","frontal")
    path.unlink()
    args.validate_only = True
    with pytest.raises(FileNotFoundError):
        run_group(cfg,args,paths,root,["sub-04","sub-05"],[DEVIANT],["frontal"])
    assert not (root / "group").exists()


@pytest.mark.parametrize("field,value,match", [
    ("identity_values",[4,3,2,1],"axes"),
    ("source_manifest","/other/processing_runs/wrong/source_manifest.json","source run"),
])
def test_wrong_source_handoff_is_rejected(expanded,field,value,match):
    root, brain, cfg, args, paths, sources = expanded
    h = json.loads((root / "sub-04/ses-01/handoff.json").read_text())
    d = json.loads((root / "sub-04/ses-01/report.json").read_text())
    path = brain / "reports/sub-04/ses-01/07_matlab_drsa_handoff.json"
    saved = json.loads(path.read_text());saved[field]=value;path.write_text(json.dumps(saved))
    with pytest.raises(ValueError,match=match):
        validate_brain_inputs(brain,h,d,[DEVIANT],sources,["frontal"])


def test_worker_threads_cannot_exceed_allocation(expanded, monkeypatch):
    """Reject oversubscription before participant work or output creation."""
    import run_drsa_dependence_cohort as cohort
    root, brain, *_ = expanded
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setattr(sys, "argv", ["cohort", "--brain-results-root", str(brain),
        "--diagnostics-root", str(brain / "new-output"), "--workers", "2", "--threads", "2"])
    with pytest.raises(SystemExit) as error:
        cohort.main()
    assert error.value.code == 2
    assert not (brain / "new-output").exists()


def test_selected_channels_and_smoothing_are_checked_before_cached_use(expanded):
    """A cache cannot bypass matching readable channel labels or smoothing."""
    root, brain, cfg, args, paths, sources = expanded
    handoff = json.loads((root / "sub-04/ses-01/handoff.json").read_text())
    report = json.loads((root / "sub-04/ses-01/report.json").read_text())
    with pytest.raises(ValueError, match="axes"):
        validate_brain_inputs(brain, handoff, report, [DEVIANT], sources, ["frontal"], {},
                              {DEVIANT: {"frontal": ["wrong-channel"]}})
    report["smooth_window_samples"] = 7
    with pytest.raises(ValueError, match="smoothing"):
        validate_brain_inputs(brain, handoff, report, [DEVIANT], sources, ["frontal"])


def test_parallel_summaries_match_serial_with_undefined_support(monkeypatch):
    """Spawned summary workers preserve every value, support count and SEM."""
    from plot_drsa_dependence import summarize_panels
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        monkeypatch.setenv(key, "1")
    axis = np.arange(24) / 20
    rng = np.random.default_rng(7)
    matrices = [rng.standard_normal((24, 24)) for _ in range(3)]
    matrices[0][:5, :] = np.nan
    stacks = {"forward": matrices, "reverse": [matrix.T for matrix in matrices]}
    intervals = {"early": axis < .5, "late": axis >= .5}
    serial = summarize_panels(stacks, axis, intervals, 1.0, workers=1)
    parallel = summarize_panels(stacks, axis, intervals, 1.0, workers=2)
    for panel in serial:
        for key in serial[panel]:
            np.testing.assert_array_equal(serial[panel][key], parallel[panel][key])


def test_group_requires_the_same_identity_set_across_participants(expanded):
    """Locally aligned participant files cannot mix different cohort identities."""
    from scipy.io import loadmat
    root, brain, cfg, args, paths, sources = expanded
    directory = root / "sub-05/ses-01"
    handoff = json.loads((directory / "handoff.json").read_text())
    handoff["identity_values"] = [2, 3, 4, 5]
    (directory / "handoff.json").write_text(json.dumps(handoff))
    (brain / "reports/sub-05/ses-01/07_matlab_drsa_handoff.json").write_text(json.dumps(handoff))
    path = directory / f"{DEVIANT}.mat"
    data = loadmat(path, simplify_cells=True)
    data = {key:value for key,value in data.items() if not key.startswith("__")}
    data["identities"] = handoff["identity_values"]
    savemat(path, data, long_field_names=True)
    for source in sources:
        path = brain_result_path(brain / "analysis", "sub-05", "ses-01", DEVIANT, source, "frontal")
        data = loadmat(path, simplify_cells=True)
        data = {key:value for key,value in data.items() if not key.startswith("__")}
        data["handoff"] = handoff
        savemat(path, data, long_field_names=True)
    report = json.loads((directory / "report.json").read_text())
    report["sha256"] = {name:file_hash(directory/name) for name in report["sha256"]}
    (directory / "report.json").write_text(json.dumps(report))
    with pytest.raises(ValueError, match="requested group"):
        run_group(cfg, args, paths, root, ["sub-04", "sub-05"], [DEVIANT], ["frontal"])
    assert not (root / "group").exists()


def test_mnn_and_raw_handoffs_cannot_be_mixed(expanded):
    """Identical axes alone do not establish identical neural RDM metrics."""
    from drsa_dependence_inputs import validate_handoff
    root,brain,cfg,args,paths,sources=expanded
    saved=json.loads((root/'sub-04/ses-01/handoff.json').read_text())
    current=dict(saved,neural_normalization={'mode':'mnn','contract_sha256':'fixed-noise',
                                           'transforms_mat':'/processing_runs/mnn/noise.mat'})
    with pytest.raises(ValueError,match='normalization contracts differ'):
        validate_handoff(saved,current)
    with pytest.raises(ValueError,match='MNN transform does not match production'):
        read_condition(root/'sub-04/ses-01/occluded_deviant.mat',DEVIANT,['frontal'],current,sources)


def test_mnn_brain_matrix_requires_matching_saved_normalization(expanded):
    from drsa_dependence_inputs import read_brain_matrix
    from scipy.io import loadmat
    root,brain,cfg,args,paths,sources=expanded
    handoff=json.loads((root/'sub-04/ses-01/handoff.json').read_text())
    handoff['neural_normalization']={'mode':'mnn','contract_sha256':'fixed-noise',
                                   'transforms_mat':'/processing_runs/mnn/noise.mat'}
    path=brain_result_path(brain/'analysis','sub-04','ses-01',DEVIANT,'observed','frontal')
    data={k:v for k,v in loadmat(path,simplify_cells=True).items() if not k.startswith('__')}
    data['handoff']=handoff;data['normalization_info']={'mode':'none'}
    savemat(path,data,long_field_names=True)
    with pytest.raises(ValueError,match='normalization differs'):
        read_brain_matrix(path,'sub-04',DEVIANT,'observed','frontal',
                          np.asarray(handoff['time_axis_sec']),handoff['identity_values'],handoff=handoff)
