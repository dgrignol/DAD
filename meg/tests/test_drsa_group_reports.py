"""Validate group-only Stage-07 HTML and report rebuilding from saved figures.

Usage
    python -m pytest meg/tests/test_drsa_group_reports.py -q

Synthetic derivative trees contain tiny group PNGs and saved inventories with
cluster paths. Tests cover the normal group-average CLI, the HTML-only CLI,
portable links with escaped filenames, and failure before report publication.
No MATLAB or external participant data is required.
"""

from __future__ import annotations

import base64
import html
from html.parser import HTMLParser
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
from urllib.parse import unquote

import numpy as np
import pytest
from scipy.io import savemat


REPO_ROOT = Path(__file__).resolve().parents[2]
CODE_DIR = REPO_ROOT / "meg/code/preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))
STAGE07 = importlib.import_module("07_run_matlab_drsa")
PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="
)


class ReportLinks(HTMLParser):
    """Collect image sources and links without interpreting HTML as text."""

    def __init__(self):
        super().__init__()
        self.images = []
        self.links = []

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag == "img":
            self.images.append(attributes["src"])
        if tag == "a":
            self.links.append(attributes["href"])


@pytest.fixture
def report_tree(tmp_path):
    """Create local figures with inventories naming a different cluster mount."""
    root = tmp_path / "derivatives with spaces"
    cfg = {
        "resolved_paths": {
            "project_root": str(REPO_ROOT / "meg"),
            "derivatives_analysis_dir": str(root / "analysis/sub-04/ses-01"),
            "derivatives_reports_dir": str(root / "reports/sub-04/ses-01"),
            "data_root": str(tmp_path / "raw"),
        },
        "layout": {"subject_bids_id": "sub-04", "session_bids_id": "ses-01"},
    }
    paths = STAGE07._analysis_paths(cfg)
    saved_dir = Path("/mnt/storage/tier2/morwur/derivatives/analysis/group/ses-01/drsa")
    for variant, kind in (("raw", "raw_average"), ("smoothed", "average_of_smoothed")):
        name = "group & mean" + ("_smoothed" if variant == "smoothed" else "") + ".png"
        figure = paths["group_drsa_dir"] / "group_average" / name
        figure.parent.mkdir(parents=True, exist_ok=True)
        figure.write_bytes(PNG)
        metadata = {
            "stage": "matlab_drsa_group_average",
            "status": "completed",
            "variant": variant,
            "session_bids_id": "ses-01",
            "subjects_requested": ["sub-04", "sub-05"],
            "group_drsa_dir": str(saved_dir),
            "plot_smoothing_window": 9,
            "outputs": [{
                "kind": kind,
                "condition": "always_visible",
                "figure": str(saved_dir / "group_average" / name),
                "subjects": ["sub-04", "sub-05"],
                "n_subjects": 2,
            }],
            # A missing participant figure must never be needed by this mode.
            "participant_figures": [{"figure": "/missing/participant.png"}],
            "skipped_inputs": [],
        }
        key = "group_report_json" if variant == "raw" else "group_report_smoothed_json"
        paths[key].parent.mkdir(parents=True, exist_ok=True)
        paths[key].write_text(json.dumps(metadata), encoding="utf-8")
    return cfg, paths


def _run_cli(tmp_path, derivative_root, *arguments):
    """Run the real Stage-07 entry point against a temporary local site."""
    raw = tmp_path / "raw"
    (raw / "sub04_000000/input_files").mkdir(parents=True, exist_ok=True)
    site = tmp_path / "site.yaml"
    site.write_text(
        f"project:\n  root: \"{REPO_ROOT / 'meg'}\"\n"
        f"  data_root: \"{raw}\"\n  derivatives_root: \"{derivative_root}\"\n",
        encoding="utf-8",
    )
    env = dict(os.environ, DAD_EPOCH_PREPARATION="saved_ica")
    env.pop("DAD_PROCESSING_RUN", None)
    env.pop("DAD_DRSA_DIAGNOSTICS_ROOT", None)
    result = subprocess.run(
        [
            sys.executable, str(CODE_DIR / "07_run_matlab_drsa.py"),
            "--config", str(REPO_ROOT / "meg/params/pipeline_default.yaml"),
            "--neural-normalization", "none",
            "--site-config", str(site),
            "--subject-config", str(REPO_ROOT / "meg/params/sub-04.yaml"),
            *arguments,
        ],
        env=env, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_html_only_cli_rebases_paths_and_embeds_only_group_images(report_tree, tmp_path):
    """Saved cluster metadata produces portable HTML without reading any MATs."""
    _, paths = report_tree
    folder = paths["group_report_json"].parent
    inputs = [*folder.glob("*.json"), *paths["group_drsa_dir"].rglob("*.png")]
    original = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in inputs}

    _run_cli(tmp_path, tmp_path / "derivatives with spaces", "--group-html-only", "--print-summary")

    for key in ("group_report_html", "group_report_smoothed_html"):
        report = paths[key]
        text = report.read_text()
        parser = ReportLinks()
        parser.feed(text)
        assert len(parser.images) == len(parser.links) == 1
        assert parser.images == parser.links
        assert "%20" in parser.images[0] and "%26" in parser.images[0]
        figure = (report.parent / unquote(parser.images[0])).resolve()
        assert figure.read_bytes() == PNG
        assert html.escape(str(figure)) in text
        assert "Figure file " in text
        assert "participant.png" not in text
        assert "Participant figures" not in text
        embedded = ReportLinks()
        embedded.feed(report.with_name(f"{report.stem}_embedded.html").read_text())
        assert len(embedded.images) == 1
        assert base64.b64decode(embedded.images[0].split(",", 1)[1]) == PNG
        assert embedded.links == parser.links
    assert original == {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in inputs}
    assert not list(folder.glob(".07*"))


@pytest.mark.parametrize("failure", ["missing", "participant_kind", "outside_group"])
def test_rebuild_validates_all_variants_before_replacing_html(report_tree, failure):
    """A bad smoothed input leaves every published HTML file intact."""
    cfg, paths = report_tree
    report_files = []
    for key in ("group_report_html", "group_report_smoothed_html"):
        linked = paths[key]
        embedded = linked.with_name(f"{linked.stem}_embedded.html")
        for path in (linked, embedded):
            path.write_text("existing report")
            report_files.append(path)
    source = paths["group_report_smoothed_json"]
    metadata = json.loads(source.read_text())
    if failure == "missing":
        metadata["outputs"][0]["figure"] = str(Path(metadata["group_drsa_dir"]) / "missing.png")
    elif failure == "participant_kind":
        metadata["outputs"][0]["kind"] = "participant"
    else:
        metadata["outputs"][0]["figure"] = str(Path(metadata["group_drsa_dir"]) / "../participant.png")
    source.write_text(json.dumps(metadata))
    with pytest.raises((ValueError, FileNotFoundError)):
        STAGE07.rebuild_group_html(cfg)
    assert all(path.read_text() == "existing report" for path in report_files)
    assert not list(source.parent.glob(".07*"))


def test_normal_group_average_cli_renders_group_only(tmp_path):
    """The numeric group-average entry point uses the same group-only renderer."""
    root = tmp_path / "derivatives"
    for subject in ("sub-04", "sub-05"):
        directory = root / "analysis" / subject / "ses-01/drsa"
        result = directory / "matlab_real_neural/A_dRSA_alt/always_visible/corr/results/alt_occipital.mat"
        result.parent.mkdir(parents=True)
        savemat(result, {"dRSA_neural_position": np.ones((3, 3)), "time_axis_sec": np.arange(3) / 120})
        (directory / "participant.png").write_bytes(PNG)
        (directory / "participant_smoothed.png").write_bytes(PNG)
    _run_cli(
        tmp_path, root, "--group-average-only", "--group-subjects", "04", "05", "--smoothing-plots",
    )
    folder = root / "reports/group/ses-01"
    for source in folder.glob("*.json"):
        metadata = json.loads(source.read_text())
        assert len(metadata["outputs"]) == 1
        assert len(metadata["participant_figures"]) == 2
    for report in folder.glob("*.html"):
        parser = ReportLinks()
        text = report.read_text()
        parser.feed(text)
        assert len(parser.images) == 1
        assert "Figure file " in text
        assert "participant.png" not in text
        assert "participant_smoothed.png" not in text
