"""
test_preprocessing_stage_smoke.py

Purpose
Run lightweight synthetic artifact-annotation and ICA-application checks. These
tests validate channel/metadata preservation without Maxwell filtering or full
participant computation.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_preprocessing_stage_smoke.py -q
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import mne
import numpy as np
from mne.preprocessing import ICA


CODE_DIR = Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))

from preprocessing_stage_utils import AcquisitionRunRecord  # noqa: E402


def _load_script(name: str, filename: str):
    """Import a numbered stage under a stable test module name."""
    spec = importlib.util.spec_from_file_location(name, CODE_DIR / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MAXWELL = _load_script("dad_maxwell_smoke", "03_A_MaxFilt.py")
ARTIFACT = _load_script("dad_artifact_smoke", "03_B_ArtefactAnnotation.py")
ICA_STAGE = _load_script("dad_ica_smoke", "03_C_ICA.py")


def _synthetic_raw(seed: int = 7) -> mne.io.RawArray:
    """Build short MEG+MISC data suitable for annotation/ICA smoke tests."""
    rng = np.random.default_rng(seed)
    sfreq = 400.0
    ch_names = [
        "MEG0011",
        "MEG0021",
        "MEG0031",
        "MEG0041",
        "MISC001",
        "MISC002",
    ]
    ch_types = ["mag", "mag", "mag", "mag", "misc", "misc"]
    info = mne.create_info(ch_names, sfreq=sfreq, ch_types=ch_types)
    data = rng.normal(scale=1e-13, size=(len(ch_names), int(sfreq * 8)))
    data[4:] = rng.normal(scale=1e-3, size=data[4:].shape)
    raw = mne.io.RawArray(data, info, verbose="ERROR")
    raw.info["bads"] = ["MEG0041"]
    return raw


def test_artifact_annotation_preserves_misc_and_bads(tmp_path):
    """Physiological annotation must not discard eye-timing MISC channels."""
    raw = _synthetic_raw()
    input_path = tmp_path / "sub02_run01_130526_raw_sss.fif"
    output_fif = tmp_path / "sub02_run01_130526_ann.fif"
    output_csv = tmp_path / "sub02_run01_130526_ann.csv"
    raw.save(input_path, overwrite=True, verbose="ERROR")
    record = AcquisitionRunRecord(
        source_path="/raw/sub02_run01_130526.fif",
        source_filename="sub02_run01_130526.fif",
        source_stem="sub02_run01_130526",
        acquisition_run_label=1,
        raw_relative_path="sub02/meg/sub02_run01_130526.fif",
    )
    report, threshold, _images = ARTIFACT._annotate_one_run(
        record=record,
        input_path=input_path,
        output_fif=output_fif,
        output_csv=output_csv,
        artifact_cfg={
            "muscle_threshold": 10.0,
            "muscle_channel_type": "mag",
            "muscle_min_length_good_sec": 0.2,
            "muscle_filter_freq_hz": [110.0, 140.0],
            "eog_mode": "auto",
            "eog_channel_names": [],
            "blink_channel": None,
            "blink_tmin_sec": -0.25,
            "blink_tmax_sec": 0.25,
        },
        initial_threshold=10.0,
        interactive=False,
        show=False,
        muscle_figure_path=tmp_path / "muscle.png",
        eog_figure_path=tmp_path / "eog.png",
        overwrite=True,
    )
    saved = mne.io.read_raw_fif(output_fif, preload=False, verbose="ERROR")
    assert threshold == 10.0
    assert saved.info["bads"] == ["MEG0041"]
    assert report["misc_channels_preserved"] == ["MISC001", "MISC002"]
    assert "MISC001" in saved.ch_names

    resumed, _images = ARTIFACT._validate_resumable_annotation(
        record=record,
        input_path=input_path,
        output_fif=output_fif,
        output_csv=output_csv,
        threshold=10.0,
        artifact_cfg={
            "eog_mode": "auto",
            "eog_channel_names": [],
        },
        raw_root=tmp_path / "raw",
    )
    assert resumed["execution_status"] == "resumed_existing"
    assert resumed["misc_channels_preserved"] == ["MISC001", "MISC002"]


def test_manual_ica_apply_preserves_misc_and_run_identity(tmp_path):
    """ICA projection keeps non-MEG channels and lineage-derived filenames."""
    raw = _synthetic_raw(seed=11)
    input_path = tmp_path / "sub02_run01_130526_ann.fif"
    raw.save(input_path, overwrite=True, verbose="ERROR")
    fit_raw = raw.copy().pick("meg")
    fit_raw.filter(1.0, 40.0, verbose="ERROR")
    model = ICA(
        n_components=2,
        method="fastica",
        random_state=97,
        max_iter=100,
    )
    model.fit(fit_raw, verbose="ERROR")
    record = AcquisitionRunRecord(
        source_path="/raw/sub02_run01_130526.fif",
        source_filename="sub02_run01_130526.fif",
        source_stem="sub02_run01_130526",
        acquisition_run_label=1,
        raw_relative_path="sub02/meg/sub02_run01_130526.fif",
    )
    reports, lineage = ICA_STAGE.apply_ica_to_runs(
        model=model,
        exclusions=[],
        records=[record],
        input_paths=[input_path],
        sss_dir=tmp_path,
        raw_root=tmp_path / "raw",
        overwrite=True,
    )
    output_path = tmp_path / "sub02_run01_130526_ica.fif"
    saved = mne.io.read_raw_fif(output_path, preload=False, verbose="ERROR")
    assert reports[0]["acquisition_run_label"] == 1
    assert lineage[0]["source_filename"] == "sub02_run01_130526.fif"
    assert saved.info["bads"] == ["MEG0041"]
    assert [name for name in saved.ch_names if name.startswith("MISC")] == [
        "MISC001",
        "MISC002",
    ]

    resumed_reports, resumed_lineage = ICA_STAGE.apply_ica_to_runs(
        model=model,
        exclusions=[],
        records=[record],
        input_paths=[input_path],
        sss_dir=tmp_path,
        raw_root=tmp_path / "raw",
        overwrite=False,
        resume=True,
    )
    assert resumed_reports[0]["execution_status"] == "resumed_existing"
    assert resumed_lineage[0]["status"] == "resumed_existing"


def test_noninteractive_ica_diagnostics_use_headless_browser(
    tmp_path,
    monkeypatch,
):
    """Source diagnostics must not select Qt when no display is requested."""
    raw = _synthetic_raw(seed=19).pick("meg")
    raw.filter(1.0, 40.0, verbose="ERROR")
    model = ICA(
        n_components=2,
        method="fastica",
        random_state=97,
        max_iter=100,
    )
    model.fit(raw, verbose="ERROR")

    # Synthetic MEG channels have no sensor locations, so isolate the browser
    # regression by replacing only the unrelated topography renderer.
    monkeypatch.setattr(
        model,
        "plot_components",
        lambda show=False: plt.subplots(figsize=(4, 3))[0],
    )
    images = ICA_STAGE.generate_component_diagnostics(
        model=model,
        fit_data=raw,
        figure_dir=tmp_path,
        subject_label="sub-99",
        processing_pass="review",
        raw_root=tmp_path / "raw",
        show=False,
        include_properties=False,
    )
    assert len(images) == 2
    assert all(path.is_file() for path in images)


def test_maxwell_before_after_qc_uses_existing_fif_outputs(tmp_path):
    """Report comparison plotting must not require rerunning Maxwell filtering."""
    raw = _synthetic_raw(seed=23)
    filtered = raw.copy()
    filtered._data[:4] *= 0.5
    raw_path = tmp_path / "sub02_run01_130526_raw.fif"
    filtered_path = tmp_path / "sub02_run01_130526_raw_sss.fif"
    figure_path = tmp_path / "qc_before_after_psd.png"
    raw.save(raw_path, overwrite=True, verbose="ERROR")
    filtered.save(filtered_path, overwrite=True, verbose="ERROR")
    # Plotting moved to maxwell_qc before SQUID integration; exercise the
    # current public plot helper using saved FIFs, rather than the removed
    # private 03_A wrapper. No detector or SSS operation is required.
    from maxwell_qc import save_psds
    plots, window = save_psds(raw=mne.io.read_raw_fif(raw_path, verbose="ERROR"),
        sss=mne.io.read_raw_fif(filtered_path, verbose="ERROR"),
        source_filename="sub02_run01_130526.fif", detected_bads=["MEG0041"],
        applied_bads=["MEG0041"], prefix=tmp_path / "qc", raw_root=tmp_path / "protected_raw")
    assert figure_path.name in [row["filename"] for row in plots]
    assert "MEG0041" not in window["comparison_channels"]
    assert figure_path.is_file()
    assert figure_path.stat().st_size > 0
