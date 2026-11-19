"""Verify final ICA display scope, exact commands and non-persistent review tools.

Usage
-----
    python -B -m pytest -p no:cacheprovider meg/tests/test_final_ica_review.py -q

Synthetic Raw pieces reproduce the MNE 1.10 concatenation-length failure and
exercise real source projections/property axes. Task tests intercept source
construction only after checking accepted native anchors/window semantics.
Report tests preserve model/metadata/tuning bytes and refuse changed inputs;
no DAD participant, model decision, scheduler or production epoch is touched.
"""

from copy import deepcopy
from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import shlex
import sys
import subprocess
from types import SimpleNamespace

import matplotlib.pyplot as plt
import mne
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
import ica_review_plots as plots
import review_final_ica as review
import final_ica_native_diagnostics as native
import render_stage03_mne_reports as renderer
from continuous_cleaning import coordinates, file_hash, stage_configuration
from final_ica_review_book import write_component_book


# -----------------------------------------------------------------------------
# Real concatenated-source projections and restored shared-axis tick locators
# -----------------------------------------------------------------------------

def example():
    """Return two differently anchored pieces, their concatenation and a real ICA."""
    rng = np.random.RandomState(97)
    info = mne.create_info(["MEG0111", "MEG0121"], 100., "mag")
    pieces = [mne.io.RawArray(rng.normal(size=(2, size)) * 1e-12, info,
        first_samp=first, verbose="ERROR") for size, first in ((600, 1000), (900, 5000))]
    raw = mne.concatenate_raws([piece.copy() for piece in pieces], verbose="ERROR")
    raw.filter(1., 40., verbose="ERROR")
    model = mne.preprocessing.ICA(n_components=2, method="fastica", random_state=97, max_iter="auto")
    model.fit(raw, verbose="ERROR")
    return raw, model


def test_display_copy_preserves_all_source_samples_and_joins_without_refit():
    raw, model = example()
    before = raw.get_data().copy()
    display = plots.continuous_display_copy(raw)
    source = model.get_sources(display)
    assert source.n_times == raw.n_times == model.n_samples_ == 1500
    assert source.get_data().shape == (2, 1500)
    assert list(display.annotations.description) == list(raw.annotations.description)
    np.testing.assert_array_equal(raw.get_data(), before)
    # Independent per-piece projections have valid single-piece metadata and
    # must match every sample of the flattened display projection.
    first = raw.copy().crop(0., 5.99)
    second = raw.copy().crop(6., 14.99)
    expected = np.concatenate([model.get_sources(piece).get_data() for piece in (first, second)], axis=1)
    np.testing.assert_allclose(source.get_data(), expected, atol=1e-12)


def test_property_panel_uses_all_segments_and_restores_numeric_ticks(monkeypatch):
    raw, model = example()
    import mne.viz.ica as viz
    # Synthetic channels carry no physical sensor locations; only replace the
    # map artist, leaving real ICA projection, segmentation, PSD and means.
    monkeypatch.setattr(viz, "_plot_ica_topomap", lambda *args, **kwargs: kwargs["axes"].plot([0, 1], [0, 1]))
    display = plots.continuous_display_copy(raw)
    figure = model.plot_properties(display, picks=0, psd_args={"fmax": 40.}, show=False, verbose="ERROR")[0]
    metadata = plots.describe_properties(figure, display.info["sfreq"])
    axes = {ax.get_label(): ax for ax in figure.axes}
    assert len(axes["erp"].get_xticks()) >= 4
    assert metadata["segments_plotted"] > 3  # Three was the legacy first-piece limit.
    assert metadata["stimulus_locked"] is False
    assert metadata["last_sample_sec"] == 1.99
    assert "not stimulus-locked" in axes["image"].get_title()
    assert "ERP" not in axes["image"].get_title()
    plt.close(figure)


@pytest.mark.parametrize("date", [False, True])
def test_native_window_names_every_BAD_and_preserves_original_time_and_signals(tmp_path, date, monkeypatch):
    from test_stage4_continuous import raw_example
    raw = raw_example(date=date, duration=20.)
    raw.set_annotations(mne.Annotations([1., 7., 8.5, 12.], [.2, .5, .25, .3],
        ["BAD_muscle", "BAD_blink_auto", "BAD_muscle", "BAD_eye_movement_auto"]))
    model = mne.preprocessing.ICA(2, method="fastica", random_state=97, max_iter="auto")
    model.fit(raw, verbose="ERROR")
    before = raw.get_data().copy()
    record = SimpleNamespace(source_filename="sub08_run01_test.fif", source_stem="sub08_run01_test")
    metadata, preview = plots.context_selection(raw, record)
    assert metadata["anchor_category"] == "BAD_blink_auto"
    assert metadata["window_start_sample_raw"] == raw.first_samp + 5 * 400
    assert {row["description"] for row in metadata["bad_annotations"]} == {"BAD_blink_auto", "BAD_muscle"}
    path, plotted = plots.native_context(model, raw, record, tmp_path)
    assert path.is_file() and {key: plotted[key] for key in metadata} == metadata
    assert plotted["frontal_magnetometer_display"].startswith("per_channel_window_mean_subtracted")
    np.testing.assert_array_equal(raw.get_data(), before)
    muscle, _ = plots.context_selection(raw, record, "BAD_muscle")
    assert muscle["window_start_sample_raw"] == raw.first_samp
    # Each report command must retain the displayed centering operation. Old
    # raw-level images must not silently become centered when reproduced.
    monkeypatch.setattr(review, "discover_raw_run_records", lambda cfg: [record])
    cfg, payload = {"subject": {"id": "08"}}, {"prefit_model_sha256": "a" * 64}
    path.with_suffix(".json").write_text(json.dumps(plotted))
    recipe = review.make_plot_recipes(cfg, payload, [path])[path.name]
    assert "--center-frontal" in shlex.split(recipe["commands"]["local"])
    path.with_suffix(".json").write_text(json.dumps(metadata))
    recipe = review.make_plot_recipes(cfg, payload, [path])[path.name]
    assert "--center-frontal" not in shlex.split(recipe["commands"]["local"])


# -----------------------------------------------------------------------------
# Per-image commands, offline instructions and one human review interface
# -----------------------------------------------------------------------------

def test_every_book_image_has_pinned_commands_and_explicit_legacy_scope(tmp_path, monkeypatch):
    record = SimpleNamespace(source_filename="sub08_run01_test.fif", source_stem="sub08_run01_test")
    monkeypatch.setattr(review, "discover_raw_run_records", lambda cfg: [record])
    cfg = {"subject": {"id": "08"}}
    payload = {"subject_bids_id": "sub-08", "n_components_fitted": 1,
        "component_properties_requested": True, "prefit_model_sha256": "a" * 64,
        "fit_execution": {"n_samples_fitted": 3000},
        "fit_input_summaries": [{"intervals": [{"status": "included", "fit_samples": 1000}]}]}
    paths = []
    for name in ("sources", "topographies_01", "ICA000_properties_01"):
        path = tmp_path / f"03_C_ICA_review_{name}.png"
        figure, axis = plt.subplots(); axis.plot([0, 1], [1, 0]); figure.savefig(path); plt.close(figure)
        paths.append(path)
    payload["plot_recipes"] = review.make_plot_recipes(cfg, payload, paths)
    for recipe in payload["plot_recipes"].values():
        tokens = shlex.split(recipe["commands"]["local"])
        assert tokens[tokens.index("--expected-model") + 1] == "a" * 64
    prop = payload["plot_recipes"][paths[-1].name]
    assert prop["metadata"]["scope"] == "legacy_first_fit_stretch"
    assert "--scope full" in prop["full_commands"]["local"]
    assert "--plot task --component 0" in prop["task_commands"]["local"]
    output = tmp_path / "book.html"
    write_component_book(payload, paths, output, tmp_path / "raw")
    text = output.read_text()
    assert text.count("Reproduce locally / open interactively") == len(paths)
    assert "How to review this participant" in text and "not a task ERF" in text
    payload["review_interface_version"] = 2
    payload["scientific_acceptance"] = "separate_explicit_stage4_review"
    sections, images = renderer._ica_content(payload, tmp_path, "review")
    assert not images and "03_C_ICA_review_components.html" in sections[0][1]
    # Processing metadata must remain inspectable when a legacy property plot
    # is missing, while clearly retaining the incomplete-review requirement.
    paths[-1].unlink()
    sections, images = renderer._ica_content(payload, tmp_path, "review")
    assert not images and "diagnostics are incomplete" in "".join(text for _, text in sections)
    payload["review_interface_version"] = 1
    with pytest.raises(FileNotFoundError, match="Expected one ICA properties"):
        renderer._ica_content(payload, tmp_path, "review")


# -----------------------------------------------------------------------------
# Accepted task anchors are native; trial diagnostics never write production data
# -----------------------------------------------------------------------------

def test_task_diagnostic_uses_accepted_onsets_and_configured_window(tmp_path, monkeypatch):
    from test_stage4_continuous import raw_example
    import trial_preparation
    import visual_timing_utils
    raw = raw_example(date=True, duration=12.)
    raw.set_annotations(mne.Annotations([8.], [.1], ["BAD_sensor_context"], ch_names=[("MEG0111",)]))
    source = tmp_path / "native_raw.fif"; raw.save(source, verbose="ERROR")
    model = mne.preprocessing.ICA(2, method="fastica", random_state=97, max_iter="auto")
    model.fit(raw, verbose="ERROR")
    record = SimpleNamespace(source_filename="one.fif")
    event = raw.first_samp + 4 * 400 + 8  # accepted +20 ms, not the original pulse
    rows = [{"source_fif": "one.fif", "trial_key": "one", "condition_code": 71,
        "raw_sample_frame": event - 8, "epoch_event_sample_raw": event,
        "included_by_trial_policy": False, "final_analysis_eligible": True}]
    monkeypatch.setattr(trial_preparation, "load_prepared_decisions", lambda cfg: [{"original": True}])
    monkeypatch.setattr(visual_timing_utils, "load_reviewed_timing", lambda *args: (deepcopy(rows), {"accepted": "synthetic-gate"}))
    real_epochs = mne.Epochs
    observed = []

    def epochs(raw, events, **kwargs):
        observed.append((events.copy(), kwargs.copy()))
        return real_epochs(raw, events, **kwargs)

    monkeypatch.setattr(mne, "Epochs", epochs)
    cfg = {"layout": {"subject_bids_id": "sub-08"}, "resolved_paths": {"derivatives_root": str(tmp_path)},
        "epoching": {"event_id_map": {"always_visible": 71, "occluded_nondeviant": 81, "occluded_deviant": 91},
            "tmin_sec": -.5, "tmax_sec": 2.667, "baseline_sec": None, "detrend": 1, "reject_by_annotation": True}}
    before = file_hash(source)
    figure, metadata = review.task_figure(cfg, model, [record], [source], 0)
    assert observed[0][0].tolist() == [[event, 0, 71]]
    assert observed[0][1]["tmin"] == -.5 and observed[0][1]["tmax"] == 2.667
    assert observed[0][1]["detrend"] is None and observed[0][1]["baseline"] is None
    assert metadata["detrend"] == 1 and metadata["IC_channel_type"].startswith("misc")
    # Verify the displayed curve applies the configured linear detrend to the
    # correctly anchored component samples while preserving honest IC units.
    from scipy.signal import detrend
    stored = mne.io.read_raw_fif(source, verbose="ERROR")
    projected = model.get_sources(plots.continuous_display_copy(stored)).get_data()[0]
    left = event - stored.first_samp - 200
    stop = event - stored.first_samp + round(2.667 * 400) + 1
    expected = detrend(projected[left:stop], type="linear")
    np.testing.assert_allclose(figure.axes[3].lines[0].get_ydata(), expected, atol=1e-10)
    assert metadata["timing_gate"] == {"accepted": "synthetic-gate"}
    assert metadata["counts"]["always_visible"] == 1
    assert metadata["model_components_removed"] == [] and file_hash(source) == before
    assert list(tmp_path.glob("*-epo.fif")) == []
    assert stored.annotations.ch_names[0] == ("MEG0111",)
    plt.close(figure)
    # A batch projects the native acquisition once, independent of IC count,
    # and uses exactly the same events/operators as the single-component view.
    calls = []
    original_sources = model.get_sources
    def sources(raw):
        calls.append(raw.n_times)
        return original_sources(raw)
    monkeypatch.setattr(model, 'get_sources', sources)
    data = native.extract_native_diagnostics(cfg, model, [record], [source],
        components=[0, 1], chunk_sec=3.)
    assert len(calls) == 4 and sum(calls) == stored.n_times
    np.testing.assert_allclose(data['values'][0, 0], expected, atol=1e-10)
    assert data['metadata']['retained_trial_metadata'][0]['epoch_event_sample_raw'] == event
    assert data['windows'][0]['metadata']['window_start_sample_raw'] >= stored.first_samp
    assert data['windows'][0]['data'].shape[0] == 2
    # Default report generation uses that one extraction for both components
    # and visibly renders continuous AND task views, not optional commands only.
    cfg['resolved_paths']['data_root'] = str(tmp_path / 'read-only-raw')
    cfg['subject'] = {'id': '08'}
    record.source_stem = 'one'
    extracted = []
    def extract(*args, **kwargs):
        extracted.append(True)
        return data
    monkeypatch.setattr(native, 'extract_native_diagnostics', extract)
    paths, manifest = native.generate_native_diagnostics(cfg, model, [record], [source], tmp_path / 'native-figures', 'b'*64)
    assert len(extracted) == 1 and len(paths) == 4
    assert manifest['continuous_and_task_default'] is True
    assert all(path.is_file() and path.with_suffix('.json').is_file() for path in paths)
    payload = {'prefit_model_sha256': 'b'*64, 'review_interface_version': 3,
        'subject_bids_id': 'sub-08', 'n_components_fitted': 2, 'component_properties_requested': True}
    monkeypatch.setattr(review, 'discover_raw_run_records', lambda cfg: [record])
    payload['plot_recipes'] = review.make_plot_recipes(cfg, payload, paths)
    book = tmp_path / 'native-book.html'
    coverage = write_component_book(payload, paths, book, tmp_path / 'read-only-raw')
    text = book.read_text()
    assert '--plot continuous --component 0' in text and '--plot task --component 0' in text
    assert 'Scroll the full native acquisitions' in text
    assert coverage['missing_continuous'] == [] and coverage['missing_task'] == []
    assert coverage['complete'] is False  # Property evidence is independently required.
    assert file_hash(source) == before


def saved_review(tmp_path, monkeypatch):
    """Create a real fitted model/report in temporary storage for refresh guards."""
    from test_stage4b_adoption import fixture
    from preprocessing_stage_utils import stage_output_dirs, write_json
    args, cfg, current, plan, participant, record, oldstamp, tuning = fixture(
        tmp_path, monkeypatch, manual_empty=False, date=True)
    monkeypatch.setattr(review, "discover_raw_run_records", lambda cfg: [record])
    # Upstream policy gates have their own real integration tests; keep the
    # current test focused on exact model/source identity and non-persistence.
    monkeypatch.setattr(review, "require_review", lambda *args: {})
    dirs = stage_output_dirs(cfg)
    source = dirs["sss"] / (record.source_stem + "_ann.fif")
    stage = review.load_stage()
    fit, summaries = stage.prepare_ica_fit_data([source], cfg["ica"])
    model = mne.preprocessing.ICA(method=cfg["ica"]["method"], n_components=cfg["ica"]["n_components"],
        random_state=cfg["ica"]["random_state"], max_iter=cfg["ica"]["max_iter"])
    model.fit(fit, verbose="ERROR")
    path = dirs["ica_models"] / "sub-04_ses-01_ica_prefit.fif"
    model.save(path, verbose="ERROR")
    fingerprint = stage.build_ica_fingerprint(input_paths=[source], ica_cfg=cfg["ica"])
    metadata_path = path.with_name(path.stem + "_metadata.json")
    write_json(metadata_path, {"model_sha256": file_hash(path), "fingerprint": fingerprint}, raw_root=tmp_path / "raw")
    payload = {"subject_bids_id": "sub-04", "processing_pass": "review", "status": "review_required",
        "scientific_acceptance": "separate_explicit_stage4_review", "manual_exclude": None,
        "prefit_model": str(path), "prefit_model_sha256": file_hash(path), "fingerprint": fingerprint,
        "configuration_sha256": stage_configuration(cfg, "03_C_ICA", include_edge=False),
        "fit_input_summaries": summaries, "n_components_fitted": int(model.n_components_),
        "fit_execution": {"n_samples_fitted": int(model.n_samples_), "mapped_fit_samples": int(fit.n_times),
            "n_iterations": int(model.n_iter_), "max_iter": model.max_iter},
        "component_properties_requested": True}
    write_json(dirs["reports"] / "03_C_ICA_review_report.json", payload, raw_root=tmp_path / "raw")
    for name in ("report.html", "components.html"):
        (dirs["reports"] / f"03_C_ICA_review_{name}").write_text("old immutable human report")
    for name in ["sources", "topographies_01"] + [f"ICA{i:03d}_properties_01" for i in range(model.n_components_)]:
        figure, ax = plt.subplots(figsize=(1, 1)); ax.plot([0, 1], [0, 1])
        figure.savefig(dirs["figures"] / f"03_C_ICA_review_{name}.png"); plt.close(figure)
    args.label = "synthetic-v2"; args.regenerate = False; args.expected_model = None
    return args, cfg, source, path, metadata_path, tuning, dirs


def test_html_refresh_preserves_model_metadata_tuning_and_creates_no_receipt(tmp_path, monkeypatch):
    args, cfg, source, model, metadata, tuning, dirs = saved_review(tmp_path, monkeypatch)
    before = {path: file_hash(path) for path in (source, model, metadata, tuning)}
    monkeypatch.setattr(mne.preprocessing.ICA, "fit", lambda *args, **kwargs: pytest.fail("Report refresh must never fit"))
    monkeypatch.setattr(renderer, "render_stage_report", lambda **kwargs: None)
    result = review.refresh_participant((args, 4, "render"))
    assert result["status"] == "reports_updated_not_accepted"
    assert before == {path: file_hash(path) for path in before}
    archive = dirs["reports"] / "history/ica_review/synthetic-v2"
    assert (archive / "03_C_ICA_review_components.html").read_text() == "old immutable human report"
    book = (dirs["reports"] / "03_C_ICA_review_components.html").read_text()
    assert "--expected-model" in book and "How to review this participant" in book
    report = json.loads((dirs["reports"] / "03_C_ICA_review_report.json").read_text())
    assert report["diagnostic_refresh"]["model_refitted"] is False
    assert report["manual_exclude"] is None and report["prefit_model_sha256"] == before[model]


def test_independent_png_cannot_overwrite_a_report_through_its_json_sidecar(tmp_path, monkeypatch):
    """A new PNG basename must not silently replace an existing review receipt/report."""
    args, cfg, source, model, metadata, tuning, dirs = saved_review(tmp_path, monkeypatch)
    report = dirs["reports"] / "03_C_ICA_review_report.json"
    before = file_hash(report)
    records = review.discover_raw_run_records(cfg)
    monkeypatch.setattr(review, "validate_review", lambda *args: (cfg, json.loads(report.read_text()), records,
        [source], model, metadata))
    monkeypatch.setattr(review, "fit_display", lambda *args: pytest.fail("An occupied sidecar must fail before plotting"))
    with pytest.raises(SystemExit):
        review.main(["--subject", "4", "--plot", "properties", "--component", "0",
                     "--save", str(report.with_suffix(".png"))])
    assert file_hash(report) == before and not report.with_suffix(".png").exists()


def test_new_shared_html_waits_for_positive_stable_evidence_stamp(tmp_path, monkeypatch):
    """Reproduce SMB's transient zero mtime; never publish that stale binding."""
    path = tmp_path / "book.html"; path.write_text("generated interface")
    actual = review.file_stamp(path)
    reads, pauses = [], []

    def stamp(path):
        reads.append(path)
        return {**actual, "mtime_ns": 0} if len(reads) == 1 else actual

    monkeypatch.setattr(review, "file_stamp", stamp)
    monkeypatch.setattr(review.time, "sleep", pauses.append)
    result = review.settle_review_outputs([path])
    assert len(reads) == 3 and len(pauses) == 2
    assert result[str(path)] == actual and actual["mtime_ns"] > 0


def test_native_only_refuses_modified_cached_property_evidence(tmp_path, monkeypatch):
    """Adding native views must not adopt stale or modified full-fit pictures."""
    args, cfg, source, model, metadata, tuning, dirs = saved_review(tmp_path, monkeypatch)
    report = json.loads((dirs['reports'] / '03_C_ICA_review_report.json').read_text())
    report.update(diagnostic_scope='full_fit', component_review_coverage={'complete': True}, evidence_files=[])
    for path in dirs['figures'].glob('03_C_ICA_review_ICA*_properties_*.png'):
        path.with_suffix('.json').write_text(json.dumps({'scope':'full_fit',
            'display_samples':report['fit_execution']['n_samples_fitted']}))
        report['evidence_files'].append({'relative_path':str(path.relative_to(Path(cfg['resolved_paths']['derivatives_root']))),
            'stamp':review.file_stamp(path)})
    review.validate_full_property_evidence(cfg, report)
    path = next(dirs['figures'].glob('03_C_ICA_review_ICA*_properties_*.png'))
    path.write_bytes(path.read_bytes() + b'changed')
    with pytest.raises(RuntimeError, match='property evidence changed'):
        review.validate_full_property_evidence(cfg, report)


@pytest.mark.parametrize("changed", ["model", "source"])
def test_changed_inputs_stop_report_refresh_before_archiving(tmp_path, monkeypatch, changed):
    args, cfg, source, model, metadata, tuning, dirs = saved_review(tmp_path, monkeypatch)
    if changed == "model":
        with model.open("ab") as handle:
            handle.write(b"changed")
    else:
        source.touch()
    with pytest.raises(RuntimeError, match="changed|differ"):
        review.refresh_participant((args, 4, "render"))
    assert not (dirs["reports"] / "history/ica_review/synthetic-v2").exists()


def test_regeneration_requires_all_component_properties_and_preserves_existing_model(tmp_path, monkeypatch):
    """A report lacking old properties can become complete only after regeneration."""
    args, cfg, source, model, metadata, tuning, dirs = saved_review(tmp_path, monkeypatch)
    report_path = dirs["reports"] / "03_C_ICA_review_report.json"
    payload = json.loads(report_path.read_text())
    payload["component_properties_requested"] = False
    report_path.write_text(json.dumps(payload))
    args.regenerate = True
    before = {path: file_hash(path) for path in (source, model, metadata, tuning)}
    stage = review.load_stage()

    def generate(**kwargs):
        assert kwargs["include_properties"] is True
        assert kwargs["fit_data"].n_times == payload["fit_execution"]["n_samples_fitted"]
        return list(dirs["figures"].glob("03_C_ICA_review_*.png"))

    monkeypatch.setattr(stage, "generate_component_diagnostics", generate)
    monkeypatch.setattr(review, "load_stage", lambda: stage)
    monkeypatch.setattr(mne.preprocessing.ICA, "fit", lambda *args, **kwargs: pytest.fail("Report regeneration must never fit"))
    monkeypatch.setattr(renderer, "render_stage_report", lambda **kwargs: None)
    def native_figures(cfg, model, records, inputs, figure_dir, model_sha):
        from shutil import copy2
        paths = []
        for component in range(model.n_components_):
            original = figure_dir / f'03_C_ICA_review_ICA{component:03}_properties_01.png'
            for kind in ('continuous', 'task'):
                path = figure_dir / f'03_C_ICA_review_ICA{component:03}_{kind}_01.png'
                copy2(original, path)
                path.with_suffix('.json').write_text(json.dumps({'model_sha256': model_sha}))
                paths.append(path)
        manifest = {'components_complete': list(range(model.n_components_))}
        (figure_dir / '03_C_ICA_review_native_diagnostics.json').write_text(json.dumps(manifest))
        return paths, manifest
    monkeypatch.setattr(native, 'generate_native_diagnostics', native_figures)
    review.refresh_participant((args, 4, "render"))
    updated = json.loads(report_path.read_text())
    assert updated["component_properties_requested"] is True
    assert updated["component_review_coverage"]["complete"] is True
    assert updated["diagnostic_scope"] == "full_fit"
    assert updated["diagnostic_refresh"]["figures_regenerated"] is True
    assert before == {path: file_hash(path) for path in before}
    assert list((dirs["reports"] / "history/ica_review/synthetic-v2/figures").glob("*.png"))


@pytest.mark.parametrize("participant,success,mode", [(8, True, 'native'), (8, True, 'full'), (15, False, 'native')])
def test_cluster_refresh_dispatches_report_only_and_preserves_cohort_exclusion(tmp_path, participant, success, mode):
    """Use harmless runtime stubs to verify the array never calls fitting/apply."""
    binary = tmp_path / "bin"; binary.mkdir()
    project = tmp_path / "project"; project.mkdir()
    mamba = binary / "micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n"); mamba.chmod(0o755)
    python = binary / "python"
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$*" > "$DAD_TEST_OUTPUT"\n'); python.chmod(0o755)
    calls = tmp_path / "calls.txt"
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
        DAD_PROJECT_ROOT=str(project), DAD_ICA_REPORT_LABEL="test-v2", MICROMAMBA_BIN=str(mamba),
        SLURM_ARRAY_TASK_ID=str(participant), SLURM_CPUS_PER_TASK="4", DAD_TEST_OUTPUT=str(calls),
        TMPDIR=str(tmp_path / "tmp"))
    env['DAD_ICA_REPORT_MODE'] = mode
    wrapper = Path(__file__).resolve().parents[1] / "cluster/submit_final_ica_reports.sh"
    result = subprocess.run(["bash", str(wrapper)], env=env, capture_output=True, text=True)
    assert (result.returncode == 0) is success, result.stdout + result.stderr
    if not success:
        assert not calls.exists()
        return
    command = shlex.split(calls.read_text())
    assert command[1].endswith("review_final_ica.py")
    assert "--refresh" in command and ('--native-only' if mode == 'native' else '--regenerate') in command
    assert command[command.index("--subject") + 1] == "8"
    assert command[command.index("--label") + 1] == "test-v2_sub8"
    assert not any("run_pipeline" in token or "03_C_ICA.py" in token or "--exclude" == token for token in command)
