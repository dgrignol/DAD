"""
test_real_preprocessing_headers.py

Purpose
Perform read-only header checks against the mounted seven-subject raw dataset.
The test skips when MORWUR is unavailable and never writes derivatives.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_real_preprocessing_headers.py -q
"""

from __future__ import annotations

from pathlib import Path
import sys

import mne
import pytest


RAW_ROOT = Path("/Volumes/MORWUR/Projects/DAMIANO/Attention/data")
CODE_DIR = Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))

import importlib.util  # noqa: E402
from preprocessing_stage_utils import parse_acquisition_run_label  # noqa: E402


def _load_maxwell_stage():
    """Import the numbered Maxwell script for pure transform selection."""
    path = CODE_DIR / "03_A_MaxFilt.py"
    spec = importlib.util.spec_from_file_location("dad_real_header_maxwell", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MAXWELL = _load_maxwell_stage()


@pytest.mark.skipif(not RAW_ROOT.is_dir(), reason="MORWUR raw dataset not mounted")
def test_all_subject_task_headers_preserve_expected_ocular_modalities():
    """Inspect every task FIF lazily and verify EOG/MISC modality facts."""
    expected_counts = {"01": 7, **{f"{value:02d}": 5 for value in range(2, 8)}}
    expected_destinations = {
        "01": 2,
        "02": 5,
        "03": 3,
        "04": 4,
        "05": 2,
        "06": 3,
        "07": 3,
    }
    for subject_id, expected_count in expected_counts.items():
        subject_dirs = sorted(RAW_ROOT.glob(f"sub{subject_id}_*"))
        assert len(subject_dirs) == 1
        fif_paths = sorted(
            (subject_dirs[0] / "meg").glob(f"sub{subject_id}_run*.fif")
        )
        assert len(fif_paths) == expected_count
        transforms = {}
        centrality_points = None
        for path in fif_paths:
            raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
            run = parse_acquisition_run_label(path.name)
            assert raw.info.get("dev_head_t") is not None
            transforms[run] = raw.info["dev_head_t"]
            if centrality_points is None:
                centrality_points = MAXWELL._centrality_head_points(raw.info)
            # MNE 1.8/1.9 omitted these embedded Polhemus points, causing the
            # automatic Maxwell-origin fit to fail despite valid raw headers.
            extra_points = [
                point
                for point in raw.info.get("dig") or []
                if int(point["kind"]) == 4
            ]
            assert len(extra_points) > 10, (
                f"{path} exposed no usable head-shape points under "
                f"MNE {mne.__version__}; DAD preprocessing requires MNE >= 1.10"
            )
            types = raw.get_channel_types()
            eog = [
                channel
                for channel, kind in zip(raw.ch_names, types, strict=True)
                if kind == "eog"
            ]
            misc = [
                channel
                for channel, kind in zip(raw.ch_names, types, strict=True)
                if kind == "misc"
            ]
            if subject_id == "01":
                assert eog == ["EOG001", "EOG002"]
                assert misc == []
            else:
                assert eog == []
                assert misc == [
                    "MISC001",
                    "MISC002",
                    "MISC003",
                    "MISC004",
                    "MISC005",
                    "MISC006",
                ]
        destination_run, _details = MAXWELL.select_central_destination_run(
            transforms,
            head_points=centrality_points,
        )
        assert destination_run == expected_destinations[subject_id]
