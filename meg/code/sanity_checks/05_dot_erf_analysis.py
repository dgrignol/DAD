#!/usr/bin/env python3
"""
05_dot_erf_analysis.py

Purpose
Run a WORK_IN_PROGRESS event-related field (ERF) sanity-check workflow for the
DAD dot-visibility conditions. The workflow consumes completed Stage-04
epoched FIF files, computes per-subject condition ERFs on one explicit MEG
sensor type and one or more Stage-07 MEGIN sensor groups, averages those ERFs
at group level, and runs paired subject-level cluster permutation tests with
MNE. The primary test remains a spatio-temporal sensor x time test; a second
channel-average 1-D temporal cluster test is also saved for each contrast.

The analysis is intentionally separate from the dot-side decoding and
time-generalization workflows. It does not derive trajectory labels and it
does not modify any decoding outputs.

Usage examples
1) Validate inputs and resolved outputs without writing files:

   python /Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/code/sanity_checks/05_dot_erf_analysis.py \
     --config /Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/params/pipeline_default.yaml \
     --group-subjects 01 02 03 04 05 06 07 \
     --validate-only --print-summary

2) Build the default magnetometer ERF sanity-check report locally:

   python /Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/code/sanity_checks/05_dot_erf_analysis.py \
     --config /Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/params/pipeline_default.yaml \
     --group-subjects 01 02 03 04 05 06 07 \
     --sensor-type mag --sensor-groups wholebrain occipital \
     --baseline -0.1 0.0 --tmin -0.1 --tmax 2.667 \
     --analysis-sfreq 200 --overwrite

3) Run the same workflow on planar gradiometers:

   python /Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/code/sanity_checks/05_dot_erf_analysis.py \
     --config /Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/params/pipeline_default.yaml \
     --group-subjects 01 02 03 04 05 06 07 \
     --sensor-type grad --sensor-groups wholebrain occipital \
     --analysis-sfreq 200 --overwrite

4) Tune topomap latencies in seconds and include the verbose cluster table:

   python .../05_dot_erf_analysis.py --config .../pipeline_default.yaml \
     --group-subjects 01 02 03 04 05 06 07 \
     --sensor-type mag --sensor-groups occipital \
     --topomap-times 0.100 0.200 0.500 0.750 0.900 \
     --include-cluster-table --overwrite

5) Run against cluster paths from the deployed DAD directory:

   python meg/code/sanity_checks/05_dot_erf_analysis.py \
     --config meg/params/pipeline_default.yaml \
     --site-config meg/params/pipeline_cluster.yaml \
     --group-subjects 01 02 03 04 05 06 07 --overwrite

Inputs
- Completed Stage-04 epochs FIF files:
  derivatives/analysis/sub-XX/ses-01/meg/sub-XX_ses-01_task-Attention_run-all_epo.fif
- Stage-04 source manifests:
  derivatives/reports/sub-XX/ses-01/source_manifest.json
- Pipeline YAML configuration for subject/event map/path resolution.

Outputs
- Per-subject ERF arrays, trial-count summaries, JSON summaries, and diagnostic
  figures under:
  derivatives/sanity_checks/sub-XX/ses-01/dot_erf/
- Group ERF arrays, contrast arrays, channel-average contrast arrays, MNE
  cluster-stat outputs, and figures under:
  derivatives/sanity_checks/group/ses-01/dot_erf/
- Linked and fully embedded group HTML reports under:
  derivatives/reports/group/ses-01/dot_erf/

Scientific policy
- Stage-04 epochs must report ``baseline=None`` in source_manifest.json.
  Baseline correction is applied in this script and recorded in every summary
  and report caption.
- One MEG channel type is analyzed per run: either magnetometers (T; displayed
  as fT) or planar gradiometers (T/m; displayed as fT/cm). The workflow never
  mixes magnetometer and gradiometer units.
- Sensor groups reuse the Stage-07 MEGIN/Neuromag lobe selections in
  meg/code/analysis/drsa_stage_utils.m. The selected Stage-07 group is filtered
  to the requested MNE channel type and then intersected with good channels
  across included subjects, preserving channel order from the first included
  subject.
- ``occluded_combined`` is the equal-weight average of each subject's
  ``occluded_nondeviant`` and ``occluded_deviant`` ERFs. It is not a pooled
  trial average, so nondeviant/deviant count imbalance cannot make one
  occluded subtype dominate the combined condition.
- Cluster tests use subject-level paired differences and run two one-sample
  permutation tests against zero:
  1. MNE spatio-temporal sensor x time cluster test.
  2. MNE 1-D temporal cluster test after signed arithmetic averaging over the
     selected channels within subject.
  The paired contrasts are:
  1. ``occluded_deviant - occluded_nondeviant``
  2. ``occluded_combined - always_visible``
- Group condition overlays show only the three source conditions by default,
  omit individual-subject traces, and use across-subject SEM shading unless
  ``--group-overlay-band std`` or ``--group-overlay-band none`` is requested.
- Trial-count equalization is triggered only when a configured within-subject
  source-count imbalance threshold is exceeded. The default threshold is a
  pragmatic 1.25 majority/minority ratio: this corresponds to about a 12%
  difference in the expected averaging-noise standard deviation under the
  usual sqrt(n) approximation, while avoiding unnecessary resampling for the
  currently modest imbalances.

WORK_IN_PROGRESS interpretation warning
The generated report labels the analysis as WORK_IN_PROGRESS because
stimulus/path locations may differ across conditions, location imbalance may
drive ERF differences, path/location distributions need future plots and
comparisons, and inferential interpretation should remain cautious until those
checks exist.
"""

from __future__ import annotations

import argparse
import base64
from collections import Counter
import html
import json
import mimetypes
import os
from pathlib import Path
import re
import sys
from typing import Any, Iterable, Mapping

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
import mne  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402


PREPROCESSING_DIR = Path(__file__).resolve().parents[1] / "preprocessing_basic"
sys.path.insert(0, str(PREPROCESSING_DIR))

from drsa_result_utils import occlusion_timing_from_config  # noqa: E402
from pipeline_config import (  # noqa: E402
    assert_output_path_outside_raw,
    load_and_resolve_config,
)


CONDITION_ORDER = (
    "always_visible",
    "occluded_nondeviant",
    "occluded_deviant",
    "occluded_combined",
)
SOURCE_CONDITIONS = CONDITION_ORDER[:3]
DEFAULT_SUBJECTS = ("01", "02", "03", "04", "05", "06", "07")
FALLBACK_TOPOMAP_TIMES = (0.200, 0.975, 1.175, 1.783)
DEFAULT_SENSOR_GROUPS = ("wholebrain", "occipital")
KNOWN_SENSOR_GROUPS = ("wholebrain", "occipital", "parietal", "frontal", "temporal")
GROUP_OVERLAY_CONDITIONS = SOURCE_CONDITIONS
WORK_IN_PROGRESS_TEXT = (
    "WORK_IN_PROGRESS: stimulus locations may differ across conditions; "
    "location imbalance may drive ERF differences; path/location distributions "
    "per condition need future plots and comparisons; inferential "
    "interpretation should remain cautious until those checks exist."
)

MEGIN_SENSOR_GROUP_RE = re.compile(r"sensor_groups\.(\w+)\s*=\s*\{([^}]*)\};")
MEG_CHANNEL_RE = re.compile(r"'(MEG\d{4})'")

SENSOR_DISPLAY = {
    "mag": {
        "label": "magnetometers",
        "si_unit": "T",
        "plot_unit": "fT",
        "scale": 1e15,
    },
    "grad": {
        "label": "planar gradiometers",
        "si_unit": "T/m",
        "plot_unit": "fT/cm",
        "scale": 1e13,
    },
}

CONDITION_LABELS = {
    "always_visible": "Always visible",
    "occluded_nondeviant": "Occluded nondeviant",
    "occluded_deviant": "Occluded deviant",
    "occluded_combined": "Occluded combined",
}
CONDITION_COLORS = {
    "always_visible": "#1f77b4",
    "occluded_nondeviant": "#2ca02c",
    "occluded_deviant": "#d62728",
    "occluded_combined": "#9467bd",
}
CONTRAST_SPECS = (
    {
        "key": "occluded_deviant-minus-occluded_nondeviant",
        "label": "Occluded deviant - occluded nondeviant",
        "positive_condition": "occluded_deviant",
        "negative_condition": "occluded_nondeviant",
        "relevant_occlusion": True,
    },
    {
        "key": "occluded_combined-minus-always_visible",
        "label": "Occluded combined - always visible",
        "positive_condition": "occluded_combined",
        "negative_condition": "always_visible",
        "relevant_occlusion": True,
    },
)


# -----------------------------------------------------------------------------
# CLI and configuration helpers
# -----------------------------------------------------------------------------


def parse_n_permutations(value: str) -> int | str:
    """
    Parse MNE cluster permutation count.

    Accepted values are ``all`` for exact sign flips when MNE can enumerate
    them, or a positive integer for random permutations.
    """
    text = str(value).strip().lower()
    if text == "all":
        return "all"
    try:
        parsed = int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "--n-permutations must be 'all' or a positive integer."
        ) from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError(
            "--n-permutations must be 'all' or a positive integer."
        )
    return parsed


def normalize_subject_label(value: str) -> str:
    """Normalize subject input such as ``5`` or ``sub-05`` to ``sub-05``."""
    text = str(value).strip()
    if text.lower().startswith("sub-"):
        text = text[4:]
    elif text.lower().startswith("sub"):
        text = text[3:]
    if not text.isdigit():
        raise ValueError(f"Subject must be numeric or sub-XX; got {value!r}.")
    return f"sub-{int(text):02d}"


def subject_id_from_label(subject_label: str) -> str:
    """Return the two-digit subject id from a normalized ``sub-XX`` label."""
    return normalize_subject_label(subject_label)[4:]


def normalize_sensor_group(value: str) -> str:
    """
    Normalize Stage-07 sensor-group selectors to lowercase output tags.

    Stage 07 accepts ``wholeBrain``/``whole_brain`` style spellings. The ERF
    workflow mirrors that convention so users can reuse the same group names
    across dRSA and ERF sanity checks.
    """
    text = str(value).strip().lower().replace("_", "")
    if text == "wholebrain":
        return "wholebrain"
    if text == "pariatel":
        return "parietal"
    if text in KNOWN_SENSOR_GROUPS:
        return text
    raise ValueError(
        f"Unknown sensor group {value!r}; expected one of {', '.join(KNOWN_SENSOR_GROUPS)}."
    )


def normalize_sensor_groups(values: Iterable[str]) -> list[str]:
    """Normalize sensor-group CLI values while preserving first occurrence order."""
    normalized: list[str] = []
    for value in values:
        group = normalize_sensor_group(str(value))
        if group not in normalized:
            normalized.append(group)
    if not normalized:
        raise ValueError("At least one sensor group is required.")
    return normalized


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the command-line parser for the ERF sanity-check workflow."""
    meg_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="WORK_IN_PROGRESS dot-condition ERF sanity-check analysis."
    )
    parser.add_argument(
        "--config",
        default=str(meg_root / "params" / "pipeline_default.yaml"),
        help="Shared pipeline default YAML.",
    )
    parser.add_argument("--site-config", default=None)
    parser.add_argument(
        "--subject-config",
        default=None,
        help=(
            "Optional subject selector used for shared path resolution. Per-"
            "subject processing uses meg/params/sub-XX.yaml for each requested "
            "group subject."
        ),
    )
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument(
        "--group-subjects",
        nargs="+",
        default=list(DEFAULT_SUBJECTS),
        help="Subjects to include, e.g. 01 02 sub-03.",
    )
    parser.add_argument(
        "--sensor-type",
        choices=("mag", "grad"),
        default="mag",
        help="MEG channel type to analyze. Magnetometers and gradiometers are never mixed.",
    )
    parser.add_argument(
        "--sensor-groups",
        nargs="+",
        default=list(DEFAULT_SENSOR_GROUPS),
        help=(
            "Stage-07 MEGIN sensor groups to run. Defaults to wholebrain and "
            "occipital; also supports parietal, frontal, and temporal."
        ),
    )
    parser.add_argument(
        "--baseline",
        nargs=2,
        type=float,
        default=[-0.1, 0.0],
        metavar=("TMIN", "TMAX"),
        help="Baseline interval in seconds, applied inside this script.",
    )
    parser.add_argument("--tmin", type=float, default=-0.1)
    parser.add_argument("--tmax", type=float, default=2.667)
    parser.add_argument(
        "--analysis-sfreq",
        type=float,
        default=200.0,
        help="Sampling frequency for ERF arrays and cluster statistics.",
    )
    parser.add_argument(
        "--topomap-times",
        nargs="+",
        type=float,
        default=None,
        help=(
            "Latency values in seconds used for scalp topomaps. Defaults to "
            "0.2 s, 0.1 s before/after first complete occlusion, and 0.2 s "
            "after first reappearance."
        ),
    )
    parser.add_argument(
        "--group-overlay-band",
        choices=("sem", "std", "none"),
        default="sem",
        help=(
            "Across-subject variability band for group condition overlays. "
            "The default sem follows the Stage-08 across-subject SEM convention."
        ),
    )
    parser.add_argument("--cluster-forming-p", type=float, default=0.05)
    parser.add_argument("--cluster-alpha", type=float, default=0.05)
    parser.add_argument(
        "--tail",
        choices=("two-sided", "greater", "less"),
        default="two-sided",
        help="MNE cluster-test tail for signed ERF contrasts.",
    )
    parser.add_argument(
        "--n-permutations",
        type=parse_n_permutations,
        default="all",
        help="MNE permutations: 'all' for exact sign flips or a positive integer.",
    )
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument(
        "--imbalance-ratio-threshold",
        type=float,
        default=1.25,
        help=(
            "Within-subject source trial-count ratio above which equalized "
            "sensitivity analyses are run."
        ),
    )
    parser.add_argument(
        "--equalization-repeats",
        type=int,
        default=100,
        help="Repeated fixed-seed random subsampling draws for sensitivity analyses.",
    )
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument(
        "--include-cluster-table",
        action="store_true",
        help="Include the verbose cluster-results table in the HTML report.",
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--print-summary", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser


def validate_cli_args(args: argparse.Namespace) -> None:
    """Fail early for invalid numeric and overwrite/resume options."""
    try:
        args.sensor_groups = normalize_sensor_groups(args.sensor_groups)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    if args.resume and args.overwrite:
        raise SystemExit("--resume and --overwrite are mutually exclusive.")
    if args.tmax <= args.tmin:
        raise SystemExit("--tmax must be greater than --tmin.")
    baseline_start, baseline_end = [float(value) for value in args.baseline]
    if baseline_end <= baseline_start:
        raise SystemExit("--baseline requires TMAX greater than TMIN.")
    if args.analysis_sfreq <= 0:
        raise SystemExit("--analysis-sfreq must be positive.")
    if not 0 < args.cluster_forming_p < 1:
        raise SystemExit("--cluster-forming-p must be between 0 and 1.")
    if not 0 < args.cluster_alpha < 1:
        raise SystemExit("--cluster-alpha must be between 0 and 1.")
    if args.n_jobs < 1:
        raise SystemExit("--n-jobs must be at least 1.")
    if args.imbalance_ratio_threshold < 1:
        raise SystemExit("--imbalance-ratio-threshold must be at least 1.")
    if args.equalization_repeats < 1:
        raise SystemExit("--equalization-repeats must be at least 1.")


def _subject_config_path(
    args: argparse.Namespace,
    *,
    subject_label: str,
    n_requested_subjects: int,
) -> Path:
    """
    Resolve the subject YAML for one participant.

    A user-supplied ``--subject-config`` is honored only in single-subject runs
    to avoid accidentally applying one participant selector to every requested
    group subject.
    """
    if args.subject_config and n_requested_subjects == 1:
        return Path(args.subject_config).expanduser().resolve()
    meg_root = Path(__file__).resolve().parents[2]
    return meg_root / "params" / f"{normalize_subject_label(subject_label)}.yaml"


def load_subject_config(
    args: argparse.Namespace,
    *,
    subject_label: str,
    n_requested_subjects: int,
) -> tuple[dict[str, Any], list[str]]:
    """Load resolved pipeline configuration for one requested subject."""
    subject_config = _subject_config_path(
        args,
        subject_label=subject_label,
        n_requested_subjects=n_requested_subjects,
    )
    cfg, warnings, _config_path, _subject_path = load_and_resolve_config(
        args.config,
        subject_config,
        args.tuning_config,
        site_config_path=args.site_config,
        required_sections=("epoching",),
    )
    return cfg, warnings


def load_shared_config(args: argparse.Namespace) -> tuple[dict[str, Any], list[str]]:
    """Load a representative config used only for shared group/report paths."""
    meg_root = Path(__file__).resolve().parents[2]
    subject_config = (
        Path(args.subject_config).expanduser().resolve()
        if args.subject_config
        else meg_root / "params" / "sub-01.yaml"
    )
    cfg, warnings, _config_path, _subject_path = load_and_resolve_config(
        args.config,
        subject_config,
        args.tuning_config,
        site_config_path=args.site_config,
        required_sections=("epoching",),
    )
    return cfg, warnings


def resolve_topomap_times(args: argparse.Namespace, cfg: Mapping[str, Any]) -> None:
    """
    Resolve default topomap latencies from configured occlusion timing.

    Users can still pass explicit seconds through ``--topomap-times``. When the
    option is omitted, the defaults target 200 ms after trial start, 100 ms
    before/after first complete occlusion, and 200 ms after first reappearance.
    """
    if args.topomap_times is not None:
        args.topomap_times = [float(value) for value in args.topomap_times]
        return
    try:
        timing = occlusion_timing_from_config(cfg)
        first_complete = float(timing["first_complete_occlusion_sec"])
        reappearance = float(timing["reappearance_sec"])
        args.topomap_times = [
            0.200,
            first_complete - 0.100,
            first_complete + 0.100,
            reappearance + 0.200,
        ]
    except Exception:
        args.topomap_times = list(FALLBACK_TOPOMAP_TIMES)


# -----------------------------------------------------------------------------
# Stage-07 sensor-group definitions
# -----------------------------------------------------------------------------


def stage07_sensor_group_definition_path() -> Path:
    """Return the MATLAB utility that defines the canonical MEGIN lobe groups."""
    return Path(__file__).resolve().parents[1] / "analysis" / "drsa_stage_utils.m"


def load_stage07_sensor_groups(path: Path | None = None) -> dict[str, list[str]]:
    """
    Parse the Stage-07 MEGIN/Neuromag sensor groups from ``drsa_stage_utils.m``.

    Data flow
    - Read the MATLAB utility used by Stage 07.
    - Extract the hard-coded ``sensor_groups.<name> = {...};`` channel lists.
    - Reconstruct ``wholebrain`` as the stable union of frontal, temporal,
      parietal, and occipital, matching the MATLAB implementation.
    """
    definition_path = path or stage07_sensor_group_definition_path()
    text = definition_path.read_text(encoding="utf-8")
    groups: dict[str, list[str]] = {}
    for raw_group, body in MEGIN_SENSOR_GROUP_RE.findall(text):
        group = normalize_sensor_group(raw_group)
        channels = MEG_CHANNEL_RE.findall(body)
        if channels:
            groups[group] = channels

    for group in ("frontal", "temporal", "parietal", "occipital"):
        if group not in groups:
            raise RuntimeError(
                f"Stage-07 sensor group {group!r} was not found in {definition_path}."
            )

    wholebrain: list[str] = []
    for group in ("frontal", "temporal", "parietal", "occipital"):
        for channel in groups[group]:
            if channel not in wholebrain:
                wholebrain.append(channel)
    groups["wholebrain"] = wholebrain
    return groups


def sensor_group_output_tag(sensor_group: str) -> str:
    """Return the normalized sensor-group tag used in filenames and reports."""
    return normalize_sensor_group(sensor_group)


def select_sensor_group_channels(
    subject_headers: list[Mapping[str, Any]],
    *,
    sensor_group: str,
    sensor_group_definitions: Mapping[str, list[str]],
) -> tuple[list[str], dict[str, Any]]:
    """
    Select a Stage-07 sensor group on the current MNE sensor-type layout.

    The Stage-07 lists include magnetometers and planar gradiometers together.
    This function first filters those lists to the requested sensor type already
    scanned in ``subject_headers``, then intersects the candidates with good
    channels common to every included subject.
    """
    group = normalize_sensor_group(sensor_group)
    if group not in sensor_group_definitions:
        raise RuntimeError(f"Sensor group {group!r} is not defined by Stage 07.")

    first_all = list(subject_headers[0]["sensor_names_all"])
    first_all_set = set(first_all)
    common_good = set(subject_headers[0]["sensor_names_good"])
    for header in subject_headers[1:]:
        common_good &= set(header["sensor_names_good"])

    requested = list(sensor_group_definitions[group])
    requested_set = set(requested)
    candidates = [channel for channel in first_all if channel in requested_set]
    selected = [channel for channel in candidates if channel in common_good]
    dropped_not_sensor_type = [channel for channel in requested if channel not in first_all_set]
    dropped_not_common_good = [channel for channel in candidates if channel not in common_good]

    if not selected:
        raise RuntimeError(
            f"No common-good {subject_headers[0]['sensor_type']} channels remain "
            f"for Stage-07 sensor group {group!r}."
        )

    selection = {
        "sensor_group": group,
        "sensor_group_definition_file": str(stage07_sensor_group_definition_path()),
        "sensor_group_policy": (
            "Stage-07 MEGIN/Neuromag lobe list filtered to the requested MNE "
            "sensor type, then intersected with good channels across included subjects."
        ),
        "available_sensor_groups": list(KNOWN_SENSOR_GROUPS),
        "n_stage07_group_channels": int(len(requested)),
        "stage07_group_channels": requested,
        "n_candidate_channels_for_sensor_type": int(len(candidates)),
        "candidate_channels_for_sensor_type": candidates,
        "n_selected_channels": int(len(selected)),
        "selected_channels": selected,
        "n_dropped_stage07_channels_not_this_sensor_type_or_layout": int(
            len(dropped_not_sensor_type)
        ),
        "dropped_stage07_channels_not_this_sensor_type_or_layout": dropped_not_sensor_type,
        "n_dropped_candidate_channels_not_common_good": int(len(dropped_not_common_good)),
        "dropped_candidate_channels_not_common_good": dropped_not_common_good,
    }
    return selected, selection


# -----------------------------------------------------------------------------
# Output path contract
# -----------------------------------------------------------------------------


def subject_input_paths(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """Return canonical Stage-04 input paths for one subject."""
    subject = str(cfg["layout"]["subject_bids_id"])
    session = str(cfg["layout"]["session_bids_id"])
    stem = f"{subject}_{session}_task-Attention_run-all"
    analysis_meg = Path(cfg["resolved_paths"]["derivatives_analysis_dir"]) / "meg"
    reports = Path(cfg["resolved_paths"]["derivatives_reports_dir"])
    return {
        "epochs_fif": analysis_meg / f"{stem}_epo.fif",
        "source_manifest": reports / "source_manifest.json",
    }


def subject_output_paths(
    cfg: Mapping[str, Any],
    sensor_type: str,
    sensor_group: str,
) -> dict[str, Path]:
    """Return per-subject ERF sanity-check output paths."""
    derivatives_root = Path(cfg["resolved_paths"]["derivatives_root"])
    subject = str(cfg["layout"]["subject_bids_id"])
    session = str(cfg["layout"]["session_bids_id"])
    out_dir = derivatives_root / "sanity_checks" / subject / session / "dot_erf"
    stem = f"{subject}_{session}_dot_erf_{sensor_type}_{sensor_group_output_tag(sensor_group)}"
    return {
        "out_dir": out_dir,
        "erfs_npz": out_dir / f"{stem}_erfs.npz",
        "summary_json": out_dir / f"{stem}_summary.json",
        "trial_counts_csv": out_dir / f"{stem}_trial_counts.csv",
        "condition_overlay_figure": out_dir / f"{stem}_condition_overlay.png",
        "gfp_rms_figure": out_dir / f"{stem}_gfp_rms.png",
    }


def group_output_paths(
    cfg: Mapping[str, Any],
    sensor_type: str,
    sensor_group: str,
) -> dict[str, Path]:
    """Return group-level ERF sanity-check output paths."""
    derivatives_root = Path(cfg["resolved_paths"]["derivatives_root"])
    session = str(cfg["layout"]["session_bids_id"])
    out_dir = derivatives_root / "sanity_checks" / "group" / session / "dot_erf"
    report_dir = derivatives_root / "reports" / "group" / session / "dot_erf"
    stem = f"dot_erf_{sensor_type}_{sensor_group_output_tag(sensor_group)}"
    paths: dict[str, Path] = {
        "out_dir": out_dir,
        "report_dir": report_dir,
        "group_erfs_npz": out_dir / f"{stem}_group_erfs.npz",
        "group_contrasts_npz": out_dir / f"{stem}_group_contrasts.npz",
        "group_channel_average_contrasts_npz": (
            out_dir / f"{stem}_group_channel_average_contrasts.npz"
        ),
        "group_summary_json": out_dir / f"{stem}_group_summary.json",
        "trial_counts_csv": out_dir / f"{stem}_trial_counts.csv",
        "imbalance_csv": out_dir / f"{stem}_imbalance_diagnostics.csv",
        "trial_counts_figure": out_dir / f"{stem}_trial_counts.png",
        "sensor_layout_figure": out_dir / f"{stem}_selected_sensor_layout.png",
        "condition_overlay_figure": out_dir / f"{stem}_condition_overlay_gfp_rms.png",
        "report_html": report_dir / f"{stem}_group_report.html",
        "report_embedded_html": report_dir / f"{stem}_group_report_embedded.html",
    }
    for condition in CONDITION_ORDER:
        paths[f"condition_{condition}_butterfly_topomap"] = (
            out_dir / f"{stem}_{condition}_butterfly_topomaps.png"
        )
    for contrast in CONTRAST_SPECS:
        key = str(contrast["key"])
        paths[f"contrast_{key}_butterfly"] = out_dir / f"{stem}_{key}_contrast_butterfly.png"
        paths[f"cluster_{key}_figure"] = out_dir / f"{stem}_{key}_cluster_highlight.png"
        paths[f"cluster_{key}_npz"] = out_dir / f"{stem}_{key}_cluster_stats.npz"
        paths[f"cluster_{key}_json"] = out_dir / f"{stem}_{key}_cluster_stats.json"
        paths[f"cluster_{key}_equalized_npz"] = (
            out_dir / f"{stem}_{key}_equalized_sensitivity_cluster_stats.npz"
        )
        paths[f"cluster_{key}_equalized_json"] = (
            out_dir / f"{stem}_{key}_equalized_sensitivity_cluster_stats.json"
        )
        paths[f"cluster1d_{key}_figure"] = (
            out_dir / f"{stem}_{key}_channel_average_cluster_highlight.png"
        )
        paths[f"cluster1d_{key}_npz"] = (
            out_dir / f"{stem}_{key}_channel_average_cluster_stats.npz"
        )
        paths[f"cluster1d_{key}_json"] = (
            out_dir / f"{stem}_{key}_channel_average_cluster_stats.json"
        )
        paths[f"cluster1d_{key}_equalized_npz"] = (
            out_dir / f"{stem}_{key}_channel_average_equalized_sensitivity_cluster_stats.npz"
        )
        paths[f"cluster1d_{key}_equalized_json"] = (
            out_dir / f"{stem}_{key}_channel_average_equalized_sensitivity_cluster_stats.json"
        )
    return paths


def expected_output_files(
    *,
    subject_headers: list[dict[str, Any]],
    group_paths: Mapping[str, Path],
    sensor_type: str,
    sensor_group: str,
) -> list[Path]:
    """Return files used by overwrite/resume guards."""
    outputs: list[Path] = []
    for header in subject_headers:
        outputs.extend(
            path
            for key, path in subject_output_paths(
                header["cfg"], sensor_type, sensor_group
            ).items()
            if key != "out_dir"
        )
    outputs.extend(
        path
        for key, path in group_paths.items()
        if key not in {"out_dir", "report_dir"}
        and not key.endswith("_equalized_npz")
        and not key.endswith("_equalized_json")
    )
    return outputs


# -----------------------------------------------------------------------------
# Stage-04 and header validation
# -----------------------------------------------------------------------------


def validate_stage04_epoching_contract(
    *,
    source_manifest: Path,
    epochs_fif: Path,
    cfg: Mapping[str, Any],
) -> dict[str, Any]:
    """
    Validate that consumed epochs are recorded as completed Stage-04 outputs.

    Cluster manifests may store cluster-absolute paths while local runs access
    the same derivative tree through ``/Volumes``. The check therefore requires
    completion status and basename agreement, not identical path prefixes.
    """
    manifest = json.loads(source_manifest.read_text(encoding="utf-8"))
    epoching_stage = (
        manifest.get("derivative_lineage", {})
        .get("stages", {})
        .get("epoching")
    )
    if not epoching_stage:
        raise RuntimeError(
            f"Source manifest lacks derivative_lineage.stages.epoching: {source_manifest}"
        )
    if epoching_stage.get("status") != "completed":
        raise RuntimeError(
            f"Stage-04 epoching is not completed in {source_manifest}: "
            f"{epoching_stage.get('status')!r}"
        )

    qc = epoching_stage.get("qc", {})
    manifest_epochs = Path(str(qc.get("epochs_fif", "")))
    if manifest_epochs.name != epochs_fif.name:
        raise RuntimeError(
            "Stage-04 manifest epochs basename does not match consumed epochs: "
            f"manifest={manifest_epochs.name!r}, consumed={epochs_fif.name!r}"
        )
    if qc.get("baseline_sec", None) is not None:
        raise RuntimeError(
            "Stage-04 manifest reports baseline correction, but ERF workflow "
            "expects unbaselined epochs and applies an explicit baseline itself."
        )
    configured_event_map = dict(cfg["epoching"]["event_id_map"])
    manifest_event_map = dict(qc.get("event_id_map", {}))
    if manifest_event_map != configured_event_map:
        raise RuntimeError(
            "Stage-04 manifest event_id_map differs from resolved config: "
            f"manifest={manifest_event_map}, config={configured_event_map}"
        )
    return {
        "source_manifest": str(source_manifest),
        "stage04_status": str(epoching_stage.get("status")),
        "manifest_epochs_fif": str(manifest_epochs),
        "consumed_epochs_fif": str(epochs_fif),
        "event_id_map": manifest_event_map,
        "baseline_sec": qc.get("baseline_sec", None),
        "n_final_epochs": qc.get("n_final_epochs"),
    }


def _condition_identity_counts(epochs_condition: mne.Epochs) -> dict[str, int]:
    """Return retained trial counts by sequence identity for one condition."""
    metadata = epochs_condition.metadata
    if metadata is None or "trial_sequence_identity" not in metadata:
        return {}
    counts = Counter(metadata["trial_sequence_identity"].to_numpy(dtype=int).tolist())
    return {str(int(key)): int(value) for key, value in sorted(counts.items())}


def scan_subject_header(
    args: argparse.Namespace,
    *,
    subject_label: str,
    n_requested_subjects: int,
) -> dict[str, Any]:
    """
    Inspect one subject's epochs without loading sample data.

    Data flow
    - Resolve the subject-specific config and canonical Stage-04 paths.
    - Validate manifest lineage and configured condition availability.
    - Record retained trial counts and good selected channels for the requested
      MEG sensor type.
    """
    cfg, warnings = load_subject_config(
        args,
        subject_label=subject_label,
        n_requested_subjects=n_requested_subjects,
    )
    paths = subject_input_paths(cfg)
    epochs_fif = paths["epochs_fif"]
    source_manifest = paths["source_manifest"]
    if not epochs_fif.exists():
        raise FileNotFoundError(f"Epochs FIF not found: {epochs_fif}")
    if not source_manifest.exists():
        raise FileNotFoundError(f"Source manifest not found: {source_manifest}")

    stage04 = validate_stage04_epoching_contract(
        source_manifest=source_manifest,
        epochs_fif=epochs_fif,
        cfg=cfg,
    )
    epochs = mne.read_epochs(epochs_fif, preload=False, verbose="ERROR")

    baseline_start, baseline_end = [float(value) for value in args.baseline]
    if baseline_start < float(epochs.tmin) or baseline_end > float(epochs.tmax):
        raise RuntimeError(
            f"Requested baseline {args.baseline} lies outside {epochs_fif} "
            f"epoch range [{epochs.tmin}, {epochs.tmax}]."
        )
    if min(args.tmin, baseline_start) < float(epochs.tmin) or args.tmax > float(epochs.tmax):
        raise RuntimeError(
            f"Requested analysis/baseline window [{min(args.tmin, baseline_start)}, "
            f"{args.tmax}] lies outside {epochs_fif} epoch range "
            f"[{epochs.tmin}, {epochs.tmax}]."
        )

    configured_conditions = set(str(key) for key in cfg["epoching"]["event_id_map"])
    missing_config_conditions = [
        condition for condition in SOURCE_CONDITIONS if condition not in configured_conditions
    ]
    missing_epoch_conditions = [
        condition for condition in SOURCE_CONDITIONS if condition not in epochs.event_id
    ]
    if missing_config_conditions or missing_epoch_conditions:
        raise RuntimeError(
            f"{cfg['layout']['subject_bids_id']} lacks required ERF conditions: "
            f"missing_config={missing_config_conditions}, "
            f"missing_epochs={missing_epoch_conditions}"
        )

    all_picks = mne.pick_types(
        epochs.info,
        meg=args.sensor_type,
        eeg=False,
        eog=False,
        misc=False,
        stim=False,
        exclude=[],
    )
    good_picks = mne.pick_types(
        epochs.info,
        meg=args.sensor_type,
        eeg=False,
        eog=False,
        misc=False,
        stim=False,
        exclude="bads",
    )
    all_sensor_names = [epochs.ch_names[int(idx)] for idx in all_picks]
    good_sensor_names = [epochs.ch_names[int(idx)] for idx in good_picks]
    if not good_sensor_names:
        raise RuntimeError(
            f"No good {args.sensor_type} channels found in {epochs_fif}."
        )

    trial_counts: dict[str, int] = {}
    identity_counts: dict[str, dict[str, int]] = {}
    for condition in SOURCE_CONDITIONS:
        epochs_condition = epochs[condition]
        trial_counts[condition] = int(len(epochs_condition))
        identity_counts[condition] = _condition_identity_counts(epochs_condition)

    return {
        "subject": str(cfg["layout"]["subject_bids_id"]),
        "subject_id": str(cfg["subject"]["id"]),
        "session": str(cfg["layout"]["session_bids_id"]),
        "cfg": cfg,
        "warnings": warnings,
        "epochs_fif": epochs_fif,
        "source_manifest": source_manifest,
        "stage04_validation": stage04,
        "event_id": dict(epochs.event_id),
        "epoch_tmin_sec": float(epochs.tmin),
        "epoch_tmax_sec": float(epochs.tmax),
        "epoch_sfreq_hz": float(epochs.info["sfreq"]),
        "n_channels_total": int(len(epochs.ch_names)),
        "bads": list(epochs.info.get("bads", [])),
        "sensor_type": args.sensor_type,
        "sensor_names_all": all_sensor_names,
        "sensor_names_good": good_sensor_names,
        "n_sensor_channels_all": int(len(all_sensor_names)),
        "n_sensor_channels_good": int(len(good_sensor_names)),
        "trial_counts": trial_counts,
        "identity_counts": identity_counts,
    }


def collect_subject_headers(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Scan requested subjects and return valid headers plus skipped records."""
    subjects = [normalize_subject_label(value) for value in args.group_subjects]
    valid: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    for subject in subjects:
        try:
            valid.append(
                scan_subject_header(
                    args,
                    subject_label=subject,
                    n_requested_subjects=len(subjects),
                )
            )
        except Exception as exc:
            record = {
                "subject": subject,
                "reason": f"{type(exc).__name__}: {exc}",
            }
            if not args.allow_missing:
                raise RuntimeError(
                    f"Failed to validate {subject}. Use --allow-missing to "
                    "build from available subjects."
                ) from exc
            skipped.append(record)
    if len(valid) < 2:
        raise RuntimeError("At least two valid subjects are required for group ERF statistics.")
    return valid, skipped


def common_good_channels(subject_headers: list[Mapping[str, Any]]) -> list[str]:
    """
    Return the channel intersection used for fixed-layout group arrays.

    Channel order follows the first valid subject's sensor order. Every returned
    channel is present and not marked bad in every included subject.
    """
    common = set(subject_headers[0]["sensor_names_good"])
    for header in subject_headers[1:]:
        common &= set(header["sensor_names_good"])
    ordered = [
        channel for channel in subject_headers[0]["sensor_names_all"] if channel in common
    ]
    if not ordered:
        raise RuntimeError("No common-good channels remain across included subjects.")
    return ordered


def make_reference_info(header: Mapping[str, Any], channel_names: list[str]) -> mne.Info:
    """Create picked reference Info with the ERF analysis sampling frequency."""
    epochs = mne.read_epochs(header["epochs_fif"], preload=False, verbose="ERROR")
    picks = mne.pick_channels(epochs.ch_names, include=channel_names, ordered=True)
    info = mne.pick_info(epochs.info, picks, copy=True)
    info["bads"] = []
    return info


# -----------------------------------------------------------------------------
# Trial-count imbalance diagnostics
# -----------------------------------------------------------------------------


def _safe_ratio(values: Iterable[int]) -> float:
    """Return max/min for positive counts, or infinity when any count is zero."""
    counts = [int(value) for value in values]
    if not counts or min(counts) <= 0:
        return float("inf")
    return float(max(counts) / min(counts))


def imbalance_diagnostics(
    subject_headers: list[Mapping[str, Any]],
    *,
    threshold: float,
) -> pd.DataFrame:
    """
    Build per-subject source-count imbalance diagnostics for each contrast.

    ``occluded_combined`` is not pooled in the primary analysis, but the
    sensitivity diagnostic considers always-visible, nondeviant, and deviant
    source counts because the equalized sensitivity draw samples all three
    source conditions to the same within-subject count.
    """
    rows: list[dict[str, Any]] = []
    for header in subject_headers:
        counts = dict(header["trial_counts"])
        dev_nondev_ratio = _safe_ratio(
            [counts["occluded_deviant"], counts["occluded_nondeviant"]]
        )
        rows.append(
            {
                "subject": header["subject"],
                "contrast": "occluded_deviant-minus-occluded_nondeviant",
                "counts_in_ratio": json.dumps(
                    {
                        "occluded_deviant": counts["occluded_deviant"],
                        "occluded_nondeviant": counts["occluded_nondeviant"],
                    },
                    sort_keys=True,
                ),
                "majority_to_minority_ratio": dev_nondev_ratio,
                "threshold": float(threshold),
                "trigger_equalized_sensitivity": bool(dev_nondev_ratio > threshold),
            }
        )
        combined_ratio = _safe_ratio(
            [
                counts["always_visible"],
                counts["occluded_nondeviant"],
                counts["occluded_deviant"],
            ]
        )
        rows.append(
            {
                "subject": header["subject"],
                "contrast": "occluded_combined-minus-always_visible",
                "counts_in_ratio": json.dumps(
                    {
                        "always_visible": counts["always_visible"],
                        "occluded_nondeviant": counts["occluded_nondeviant"],
                        "occluded_deviant": counts["occluded_deviant"],
                    },
                    sort_keys=True,
                ),
                "majority_to_minority_ratio": combined_ratio,
                "threshold": float(threshold),
                "trigger_equalized_sensitivity": bool(combined_ratio > threshold),
            }
        )
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Subject-level ERF computation
# -----------------------------------------------------------------------------


def _prepare_condition_epochs(
    epochs_all: mne.Epochs,
    *,
    condition: str,
    channel_names: list[str],
    baseline: tuple[float, float],
    load_tmin: float,
    analysis_tmin: float,
    analysis_tmax: float,
    analysis_sfreq: float,
) -> mne.Epochs:
    """
    Pick, crop, resample, and baseline-correct one condition.

    Baseline is applied after resampling and before the final analysis crop so
    the baseline interval can extend earlier than the plotted/statistical
    window.
    """
    epochs_condition_source = epochs_all[condition].copy()
    picks = mne.pick_channels(
        epochs_condition_source.ch_names,
        include=channel_names,
        ordered=True,
    )
    picked_info = mne.pick_info(epochs_condition_source.info, picks, copy=True)
    picked_info["bads"] = []

    # ``Epochs.pick`` and ``Epochs.crop`` require preloaded data in this MNE
    # version. ``get_data(picks=..., tmin=..., tmax=...)`` lets us load only the
    # common-good sensor subset and requested time window, then wrap it as an
    # EpochsArray for MNE's resampling/baseline machinery.
    data = epochs_condition_source.get_data(
        picks=picks,
        tmin=load_tmin,
        tmax=analysis_tmax,
        copy=True,
    )
    source_times = np.asarray(epochs_condition_source.times, dtype=float)
    time_mask = (source_times >= float(load_tmin)) & (source_times <= float(analysis_tmax))
    if not np.any(time_mask):
        raise RuntimeError(
            f"No samples remain for condition {condition!r} in requested window "
            f"{load_tmin}..{analysis_tmax} s."
        )
    tmin_loaded = float(source_times[time_mask][0])
    epochs_condition = mne.EpochsArray(
        data,
        picked_info,
        events=epochs_condition_source.events.copy(),
        tmin=tmin_loaded,
        event_id={condition: int(epochs_condition_source.event_id[condition])},
        baseline=None,
        metadata=(
            epochs_condition_source.metadata.copy()
            if epochs_condition_source.metadata is not None
            else None
        ),
        selection=epochs_condition_source.selection.copy(),
        verbose="ERROR",
    )
    if not np.isclose(float(epochs_condition.info["sfreq"]), float(analysis_sfreq)):
        epochs_condition.resample(float(analysis_sfreq), npad="auto")
    epochs_condition.apply_baseline(baseline)
    crop_tmin = max(float(analysis_tmin), float(epochs_condition.tmin))
    crop_tmax = min(float(analysis_tmax), float(epochs_condition.tmax))
    epochs_condition.crop(tmin=crop_tmin, tmax=crop_tmax)
    return epochs_condition


def trial_count_rows_for_subject(header: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Flatten condition and identity trial counts for CSV/report tables."""
    rows: list[dict[str, Any]] = []
    for condition in SOURCE_CONDITIONS:
        identity_counts = dict(header["identity_counts"].get(condition, {}))
        identity_values = [int(value) for value in identity_counts.values()]
        rows.append(
            {
                "subject": header["subject"],
                "condition": condition,
                "n_trials": int(header["trial_counts"][condition]),
                "n_identities": int(len(identity_counts)),
                "min_trials_per_identity": int(min(identity_values)) if identity_values else 0,
                "max_trials_per_identity": int(max(identity_values)) if identity_values else 0,
                "median_trials_per_identity": (
                    float(np.median(identity_values)) if identity_values else np.nan
                ),
                "counts_by_identity": json.dumps(identity_counts, sort_keys=True),
            }
        )
    return rows


def compute_subject_erfs(
    header: Mapping[str, Any],
    *,
    args: argparse.Namespace,
    channel_names: list[str],
    sensor_group: str,
    sensor_selection: Mapping[str, Any],
    keep_trial_data: bool,
) -> dict[str, Any]:
    """
    Compute all condition ERFs for one subject on the common-good channel set.

    Data flow
    - Read the Stage-04 epoch FIF.
    - Process each source condition independently to limit peak memory use.
    - Average retained trials per source condition.
    - Construct ``occluded_combined`` as an equal-weight mean of the two
      occluded source ERFs.
    - Optionally retain baseline-corrected trial arrays for equalized
      sensitivity analyses when imbalance diagnostics require them.
    """
    baseline = (float(args.baseline[0]), float(args.baseline[1]))
    load_tmin = min(float(args.tmin), baseline[0])
    epochs_all = mne.read_epochs(header["epochs_fif"], preload=False, verbose="ERROR")
    condition_erfs: dict[str, np.ndarray] = {}
    condition_trial_data: dict[str, np.ndarray] = {}
    times: np.ndarray | None = None

    for condition in SOURCE_CONDITIONS:
        epochs_condition = _prepare_condition_epochs(
            epochs_all,
            condition=condition,
            channel_names=channel_names,
            baseline=baseline,
            load_tmin=load_tmin,
            analysis_tmin=float(args.tmin),
            analysis_tmax=float(args.tmax),
            analysis_sfreq=float(args.analysis_sfreq),
        )
        current_times = epochs_condition.times.copy()
        if times is None:
            times = current_times
        elif len(current_times) != len(times) or not np.allclose(current_times, times):
            raise RuntimeError(
                f"{header['subject']} {condition} time axis does not match prior conditions."
            )
        data = epochs_condition.get_data(copy=True)
        condition_trial_data[condition] = data if keep_trial_data else np.empty((0,))
        condition_erfs[condition] = data.mean(axis=0)

    assert times is not None
    condition_erfs["occluded_combined"] = 0.5 * (
        condition_erfs["occluded_nondeviant"] + condition_erfs["occluded_deviant"]
    )

    erfs = np.stack([condition_erfs[condition] for condition in CONDITION_ORDER], axis=0)
    trial_count_rows = trial_count_rows_for_subject(header)
    paths = subject_output_paths(header["cfg"], args.sensor_type, sensor_group)
    paths["out_dir"].mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        paths["erfs_npz"],
        times_sec=times,
        channel_names=np.asarray(channel_names, dtype=object),
        condition_names=np.asarray(CONDITION_ORDER, dtype=object),
        erfs_si=erfs,
        sensor_type=np.asarray(args.sensor_type),
        sensor_group=np.asarray(sensor_group_output_tag(sensor_group)),
        si_unit=np.asarray(SENSOR_DISPLAY[args.sensor_type]["si_unit"]),
        baseline_sec=np.asarray(args.baseline, dtype=float),
    )
    pd.DataFrame(trial_count_rows).to_csv(paths["trial_counts_csv"], index=False)

    summary = {
        "stage": "dot_erf_subject",
        "analysis_status": "WORK_IN_PROGRESS",
        "work_in_progress_warning": WORK_IN_PROGRESS_TEXT,
        "subject": header["subject"],
        "session": header["session"],
        "epochs_fif": str(header["epochs_fif"]),
        "source_manifest": str(header["source_manifest"]),
        "stage04_validation": header["stage04_validation"],
        "conditions": list(CONDITION_ORDER),
        "source_conditions": list(SOURCE_CONDITIONS),
        "trial_counts": dict(header["trial_counts"]),
        "identity_counts": dict(header["identity_counts"]),
        "occluded_combined_construction": (
            "0.5 * ERF(occluded_nondeviant) + 0.5 * ERF(occluded_deviant); "
            "source trials are not pooled."
        ),
        "sensor_type": args.sensor_type,
        "sensor_type_label": SENSOR_DISPLAY[args.sensor_type]["label"],
        "sensor_group": sensor_group_output_tag(sensor_group),
        "sensor_selection": dict(sensor_selection),
        "si_unit": SENSOR_DISPLAY[args.sensor_type]["si_unit"],
        "plot_unit": SENSOR_DISPLAY[args.sensor_type]["plot_unit"],
        "unit_scale_for_figures": float(SENSOR_DISPLAY[args.sensor_type]["scale"]),
        "n_common_good_channels": int(len(channel_names)),
        "common_good_channels": list(channel_names),
        "baseline_interval_sec": [float(args.baseline[0]), float(args.baseline[1])],
        "time_window_sec": [float(times[0]), float(times[-1])],
        "analysis_sfreq_hz": float(1.0 / np.median(np.diff(times))),
        "n_times": int(len(times)),
        "outputs": {key: str(value) for key, value in paths.items() if key != "out_dir"},
    }
    paths["summary_json"].write_text(json.dumps(_jsonify(summary), indent=2), encoding="utf-8")

    result = {
        "header": header,
        "paths": paths,
        "times": times,
        "erfs": erfs,
        "trial_data": condition_trial_data if keep_trial_data else {},
        "summary": summary,
        "trial_count_rows": trial_count_rows,
    }
    return result


def load_subject_erfs_from_outputs(
    header: Mapping[str, Any],
    *,
    args: argparse.Namespace,
    channel_names: list[str],
    sensor_group: str,
) -> dict[str, Any]:
    """
    Load completed subject ERF outputs for partial resume runs.

    This is used when subject-level ERFs were written successfully but a later
    group-level step failed. It avoids repeating the expensive FIF reads while
    still validating that the saved subject arrays match the current sensor,
    baseline, condition, and channel contract.
    """
    paths = subject_output_paths(header["cfg"], args.sensor_type, sensor_group)
    required = (
        paths["erfs_npz"],
        paths["summary_json"],
        paths["trial_counts_csv"],
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Cannot resume {header['subject']}; missing {missing}")

    with np.load(paths["erfs_npz"], allow_pickle=True) as payload:
        times = np.asarray(payload["times_sec"], dtype=float)
        saved_channels = [str(value) for value in payload["channel_names"].tolist()]
        saved_conditions = [str(value) for value in payload["condition_names"].tolist()]
        erfs = np.asarray(payload["erfs_si"], dtype=float)
        saved_sensor_type = str(np.asarray(payload["sensor_type"]).item())
        saved_sensor_group = str(np.asarray(payload["sensor_group"]).item())
        saved_baseline = np.asarray(payload["baseline_sec"], dtype=float).tolist()

    if saved_channels != list(channel_names):
        raise RuntimeError(
            f"Cannot resume {header['subject']}; saved channel list differs from current common-good layout."
        )
    if saved_conditions != list(CONDITION_ORDER):
        raise RuntimeError(
            f"Cannot resume {header['subject']}; saved condition order is {saved_conditions}."
        )
    if saved_sensor_type != str(args.sensor_type):
        raise RuntimeError(
            f"Cannot resume {header['subject']}; saved sensor type {saved_sensor_type!r} "
            f"differs from requested {args.sensor_type!r}."
        )
    if saved_sensor_group != sensor_group_output_tag(sensor_group):
        raise RuntimeError(
            f"Cannot resume {header['subject']}; saved sensor group {saved_sensor_group!r} "
            f"differs from requested {sensor_group!r}."
        )
    if not np.allclose(saved_baseline, np.asarray(args.baseline, dtype=float)):
        raise RuntimeError(
            f"Cannot resume {header['subject']}; saved baseline {saved_baseline} "
            f"differs from requested {args.baseline}."
        )

    summary = json.loads(paths["summary_json"].read_text(encoding="utf-8"))
    trial_count_rows = pd.read_csv(paths["trial_counts_csv"]).to_dict("records")
    return {
        "header": header,
        "paths": paths,
        "times": times,
        "erfs": erfs,
        "trial_data": {},
        "summary": summary,
        "trial_count_rows": trial_count_rows,
    }


# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------


def mne_tail_value(tail: str) -> int:
    """Map readable tail names to MNE tail integers."""
    if tail == "two-sided":
        return 0
    if tail == "greater":
        return 1
    if tail == "less":
        return -1
    raise ValueError(f"Unsupported tail: {tail}")


def cluster_t_threshold(
    *,
    n_subjects: int,
    cluster_forming_p: float,
    tail: str,
) -> float:
    """Return the t-statistic threshold passed to MNE cluster permutation."""
    df = int(n_subjects) - 1
    if df < 1:
        raise ValueError("At least two subjects are required.")
    if tail == "two-sided":
        return float(stats.t.ppf(1.0 - (cluster_forming_p / 2.0), df=df))
    if tail == "greater":
        return float(stats.t.ppf(1.0 - cluster_forming_p, df=df))
    if tail == "less":
        return float(-stats.t.ppf(1.0 - cluster_forming_p, df=df))
    raise ValueError(f"Unsupported tail: {tail}")


def build_contrast_values(erfs_by_subject: np.ndarray) -> dict[str, np.ndarray]:
    """Build subject-level paired contrast arrays with shape subject x channel x time."""
    condition_to_index = {condition: idx for idx, condition in enumerate(CONDITION_ORDER)}
    contrasts: dict[str, np.ndarray] = {}
    for spec in CONTRAST_SPECS:
        pos = condition_to_index[str(spec["positive_condition"])]
        neg = condition_to_index[str(spec["negative_condition"])]
        contrasts[str(spec["key"])] = erfs_by_subject[:, pos, :, :] - erfs_by_subject[:, neg, :, :]
    return contrasts


def run_mne_cluster_1samp(
    values_chan_time: np.ndarray,
    *,
    reference_info: mne.Info,
    times: np.ndarray,
    args: argparse.Namespace,
) -> dict[str, Any]:
    """
    Run an MNE one-sample spatio-temporal cluster permutation test.

    Input ``values_chan_time`` is subject x channel x time and represents a
    paired subject-level difference ERF. MNE expects observations x time x
    space, so the array is transposed before testing.
    """
    values_chan_time = np.asarray(values_chan_time, dtype=float)
    if values_chan_time.ndim != 3:
        raise ValueError("values_chan_time must be subject x channel x time.")
    if values_chan_time.shape[0] < 2:
        raise ValueError("At least two subjects are required.")

    adjacency, adjacency_names = mne.channels.find_ch_adjacency(
        reference_info,
        ch_type=args.sensor_type,
    )
    reference_names = list(reference_info["ch_names"])
    if list(adjacency_names) != reference_names:
        name_to_adjacency_idx = {
            str(name): int(idx) for idx, name in enumerate(adjacency_names)
        }
        missing = [name for name in reference_names if name not in name_to_adjacency_idx]
        if missing:
            raise RuntimeError(
                "MNE sensor adjacency does not contain all selected channels: "
                f"{missing}"
            )
        adjacency_indices = [name_to_adjacency_idx[name] for name in reference_names]
        adjacency = adjacency[adjacency_indices, :][:, adjacency_indices]
        adjacency_names = reference_names
    threshold = cluster_t_threshold(
        n_subjects=values_chan_time.shape[0],
        cluster_forming_p=float(args.cluster_forming_p),
        tail=str(args.tail),
    )
    x_mne = np.transpose(values_chan_time, (0, 2, 1))
    t_obs, clusters, cluster_p_values, h0 = mne.stats.spatio_temporal_cluster_1samp_test(
        x_mne,
        threshold=threshold,
        n_permutations=args.n_permutations,
        tail=mne_tail_value(args.tail),
        adjacency=adjacency,
        seed=int(args.random_seed),
        n_jobs=int(args.n_jobs),
        max_step=1,
        out_type="indices",
        check_disjoint=True,
        verbose="ERROR",
    )
    significant_mask = np.zeros(t_obs.shape, dtype=bool)
    cluster_rows: list[dict[str, Any]] = []
    for cluster_idx, cluster in enumerate(clusters):
        time_indices = np.asarray(cluster[0], dtype=int)
        sensor_indices = np.asarray(cluster[1], dtype=int)
        if time_indices.size == 0:
            continue
        p_value = float(cluster_p_values[cluster_idx])
        is_significant = bool(p_value <= float(args.cluster_alpha))
        if is_significant:
            significant_mask[time_indices, sensor_indices] = True

        cluster_t = t_obs[time_indices, sensor_indices]
        peak_local = int(np.nanargmax(np.abs(cluster_t)))
        peak_time_idx = int(time_indices[peak_local])
        peak_sensor_idx = int(sensor_indices[peak_local])
        mass_source = np.abs(cluster_t) if args.tail == "two-sided" else cluster_t
        cluster_rows.append(
            {
                "cluster_index": int(cluster_idx),
                "p_value": p_value,
                "significant": is_significant,
                "n_points": int(time_indices.size),
                "n_unique_times": int(np.unique(time_indices).size),
                "n_unique_sensors": int(np.unique(sensor_indices).size),
                "time_start_sec": float(times[int(time_indices.min())]),
                "time_end_sec": float(times[int(time_indices.max())]),
                "peak_time_sec": float(times[peak_time_idx]),
                "peak_sensor": str(reference_info["ch_names"][peak_sensor_idx]),
                "peak_t": float(t_obs[peak_time_idx, peak_sensor_idx]),
                "mass": float(np.nansum(mass_source)),
            }
        )

    return {
        "t_obs": t_obs,
        "clusters": cluster_rows,
        "cluster_p_values": cluster_p_values,
        "h0": h0,
        "significant_mask": significant_mask,
        "adjacency": adjacency,
        "adjacency_names": list(adjacency_names),
        "threshold": float(threshold),
        "tail_mne": int(mne_tail_value(args.tail)),
        "n_permutations_effective": (
            str(args.n_permutations)
            if isinstance(args.n_permutations, str)
            else int(args.n_permutations)
        ),
    }


def channel_average_values(values_chan_time: np.ndarray) -> np.ndarray:
    """
    Average signed subject-level contrast values over selected channels.

    The returned subject x time matrix intentionally uses an arithmetic channel
    mean, not RMS/GFP, so the 1-D permutation test answers the requested
    selected-channel average question with the same contrast sign convention as
    the sensor x time analysis.
    """
    values_chan_time = np.asarray(values_chan_time, dtype=float)
    if values_chan_time.ndim != 3:
        raise ValueError("values_chan_time must be subject x channel x time.")
    return np.mean(values_chan_time, axis=1)


def run_mne_cluster_1d(
    values_time: np.ndarray,
    *,
    times: np.ndarray,
    args: argparse.Namespace,
) -> dict[str, Any]:
    """
    Run a one-sample temporal cluster permutation test on channel-average ERFs.

    Input ``values_time`` is subject x time and represents a paired subject-level
    contrast averaged over the selected channels before inference.
    """
    values_time = np.asarray(values_time, dtype=float)
    if values_time.ndim != 2:
        raise ValueError("values_time must be subject x time.")
    if values_time.shape[0] < 2:
        raise ValueError("At least two subjects are required.")

    threshold = cluster_t_threshold(
        n_subjects=values_time.shape[0],
        cluster_forming_p=float(args.cluster_forming_p),
        tail=str(args.tail),
    )
    t_obs, clusters, cluster_p_values, h0 = mne.stats.permutation_cluster_1samp_test(
        values_time,
        threshold=threshold,
        n_permutations=args.n_permutations,
        tail=mne_tail_value(args.tail),
        seed=int(args.random_seed),
        n_jobs=int(args.n_jobs),
        out_type="indices",
        verbose="ERROR",
    )

    significant_mask = np.zeros(t_obs.shape, dtype=bool)
    cluster_rows: list[dict[str, Any]] = []
    for cluster_idx, cluster in enumerate(clusters):
        time_indices = np.asarray(
            cluster[0] if isinstance(cluster, tuple) else cluster,
            dtype=int,
        )
        if time_indices.size == 0:
            continue
        p_value = float(cluster_p_values[cluster_idx])
        is_significant = bool(p_value <= float(args.cluster_alpha))
        if is_significant:
            significant_mask[time_indices] = True

        cluster_t = np.asarray(t_obs[time_indices], dtype=float)
        peak_local = int(np.nanargmax(np.abs(cluster_t)))
        peak_time_idx = int(time_indices[peak_local])
        mass_source = np.abs(cluster_t) if args.tail == "two-sided" else cluster_t
        cluster_rows.append(
            {
                "cluster_index": int(cluster_idx),
                "p_value": p_value,
                "significant": is_significant,
                "n_points": int(time_indices.size),
                "n_unique_times": int(np.unique(time_indices).size),
                "time_start_sec": float(times[int(time_indices.min())]),
                "time_end_sec": float(times[int(time_indices.max())]),
                "peak_time_sec": float(times[peak_time_idx]),
                "peak_t": float(t_obs[peak_time_idx]),
                "mass": float(np.nansum(mass_source)),
            }
        )

    return {
        "t_obs": t_obs,
        "clusters": cluster_rows,
        "cluster_p_values": cluster_p_values,
        "h0": h0,
        "significant_mask": significant_mask,
        "threshold": float(threshold),
        "tail_mne": int(mne_tail_value(args.tail)),
        "n_permutations_effective": (
            str(args.n_permutations)
            if isinstance(args.n_permutations, str)
            else int(args.n_permutations)
        ),
        "statistic_scope": "one-sample temporal cluster on signed channel-average subject contrasts",
    }


def equalized_contrast_values(
    subject_results: list[Mapping[str, Any]],
    *,
    contrast_key: str,
    repeats: int,
    random_seed: int,
) -> np.ndarray:
    """
    Build repeated random-subsampling sensitivity contrasts.

    For ``occluded_deviant - occluded_nondeviant`` each repeat samples both
    source conditions to the smaller count. For ``occluded_combined - always``
    each repeat samples always-visible, nondeviant, and deviant to their common
    smaller count, then reconstructs the equal-weight combined ERF.
    """
    per_subject: list[np.ndarray] = []
    for subject_index, result in enumerate(subject_results):
        trial_data = result.get("trial_data", {})
        if not trial_data:
            raise RuntimeError("Equalized sensitivity requested but trial arrays were not retained.")
        rng = np.random.default_rng(int(random_seed) + (subject_index * 100003))
        repeat_diffs: list[np.ndarray] = []
        if contrast_key == "occluded_deviant-minus-occluded_nondeviant":
            dev = np.asarray(trial_data["occluded_deviant"], dtype=float)
            nondev = np.asarray(trial_data["occluded_nondeviant"], dtype=float)
            n_keep = int(min(dev.shape[0], nondev.shape[0]))
            for _repeat in range(int(repeats)):
                dev_idx = rng.choice(dev.shape[0], size=n_keep, replace=False)
                nondev_idx = rng.choice(nondev.shape[0], size=n_keep, replace=False)
                repeat_diffs.append(dev[dev_idx].mean(axis=0) - nondev[nondev_idx].mean(axis=0))
        elif contrast_key == "occluded_combined-minus-always_visible":
            always = np.asarray(trial_data["always_visible"], dtype=float)
            nondev = np.asarray(trial_data["occluded_nondeviant"], dtype=float)
            dev = np.asarray(trial_data["occluded_deviant"], dtype=float)
            n_keep = int(min(always.shape[0], nondev.shape[0], dev.shape[0]))
            for _repeat in range(int(repeats)):
                always_idx = rng.choice(always.shape[0], size=n_keep, replace=False)
                nondev_idx = rng.choice(nondev.shape[0], size=n_keep, replace=False)
                dev_idx = rng.choice(dev.shape[0], size=n_keep, replace=False)
                combined = 0.5 * (nondev[nondev_idx].mean(axis=0) + dev[dev_idx].mean(axis=0))
                repeat_diffs.append(combined - always[always_idx].mean(axis=0))
        else:
            raise ValueError(f"Unsupported contrast for equalization: {contrast_key}")
        per_subject.append(np.stack(repeat_diffs, axis=0).mean(axis=0))
    return np.stack(per_subject, axis=0)


# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------


def _scale_for_sensor(sensor_type: str) -> float:
    """Return display scale from SI units for one sensor type."""
    return float(SENSOR_DISPLAY[sensor_type]["scale"])


def _unit_for_sensor(sensor_type: str) -> str:
    """Return display unit label for one sensor type."""
    return str(SENSOR_DISPLAY[sensor_type]["plot_unit"])


def _copy_info_with_sfreq(info: mne.Info, times: np.ndarray) -> mne.Info:
    """Copy MNE Info and set sampling frequency to match an ERF time axis."""
    copied = info.copy()
    if len(times) > 1:
        sfreq = float(1.0 / np.median(np.diff(times)))
        with copied._unlock():
            copied["sfreq"] = sfreq
    copied["bads"] = []
    return copied


def make_evoked(
    data_chan_time: np.ndarray,
    *,
    reference_info: mne.Info,
    times: np.ndarray,
    comment: str,
    nave: int,
) -> mne.EvokedArray:
    """Create an EvokedArray preserving sensor geometry from reference Info."""
    info = _copy_info_with_sfreq(reference_info, times)
    return mne.EvokedArray(
        np.asarray(data_chan_time, dtype=float),
        info,
        tmin=float(times[0]),
        nave=int(nave),
        comment=comment,
    )


def channel_colors(channel_names: list[str]) -> dict[str, Any]:
    """Return stable, high-contrast colors for selected channel traces/layout dots."""
    n_channels = max(len(channel_names), 1)
    cmap = plt.get_cmap("turbo", n_channels)
    return {
        channel: cmap(index / max(n_channels - 1, 1))
        for index, channel in enumerate(channel_names)
    }


def sensor_layout_positions(
    reference_info: mne.Info,
    *,
    sensor_type: str,
) -> dict[str, tuple[float, float]]:
    """
    Return 2-D MNE layout centers keyed by channel name.

    MNE's MEG layout gives the familiar top-view sensor projection used by
    sensor-layout plots. The returned positions are in layout coordinates, not
    physical head-coordinate units.
    """
    layout = mne.channels.find_layout(reference_info, ch_type=sensor_type)
    positions: dict[str, tuple[float, float]] = {}
    for name, pos in zip(layout.names, layout.pos):
        positions[str(name)] = (float(pos[0] + pos[2] / 2.0), float(pos[1] + pos[3] / 2.0))
    return positions


def draw_sensor_layout(
    ax: plt.Axes,
    *,
    reference_info: mne.Info,
    channel_names: list[str],
    sensor_type: str,
    colors: Mapping[str, Any],
    title: str,
    show_frame: bool = True,
) -> None:
    """
    Draw a compact selected-sensor layout using the same channel colors as traces.

    The map intentionally omits channel labels by default because whole-brain
    layouts can contain nearly one hundred magnetometers or twice as many
    gradiometers. Complete channel names are recorded in JSON/report tables.
    """
    positions = sensor_layout_positions(reference_info, sensor_type=sensor_type)
    xy = np.asarray([positions[channel] for channel in channel_names if channel in positions])
    if xy.size == 0:
        ax.text(0.5, 0.5, "layout unavailable", ha="center", va="center", fontsize=7)
        ax.set_axis_off()
        return

    plot_channels = [channel for channel in channel_names if channel in positions]
    ax.scatter(
        xy[:, 0],
        xy[:, 1],
        s=18,
        c=[colors.get(channel, "#444444") for channel in plot_channels],
        edgecolors="black",
        linewidths=0.25,
        alpha=0.95,
    )
    ax.set_title(title, fontsize=8)
    ax.set_aspect("equal", adjustable="box")
    pad_x = max(float(np.ptp(xy[:, 0])) * 0.08, 0.02)
    pad_y = max(float(np.ptp(xy[:, 1])) * 0.08, 0.02)
    ax.set_xlim(float(np.min(xy[:, 0]) - pad_x), float(np.max(xy[:, 0]) + pad_x))
    ax.set_ylim(float(np.min(xy[:, 1]) - pad_y), float(np.max(xy[:, 1]) + pad_y))
    ax.set_xticks([])
    ax.set_yticks([])
    if not show_frame:
        ax.set_axis_off()


def plot_sensor_selection_layout(
    *,
    reference_info: mne.Info,
    channel_names: list[str],
    sensor_type: str,
    sensor_group: str,
    colors: Mapping[str, Any],
) -> plt.Figure:
    """Create a standalone selected-sensor layout figure for report auditing."""
    fig, ax = plt.subplots(figsize=(4.6, 4.6), constrained_layout=True)
    draw_sensor_layout(
        ax,
        reference_info=reference_info,
        channel_names=channel_names,
        sensor_type=sensor_type,
        colors=colors,
        title=f"Selected sensors | {sensor_group_output_tag(sensor_group)}",
    )
    return fig


def across_subject_band(values_subject_time: np.ndarray, mode: str) -> np.ndarray | None:
    """
    Return a NaN-safe across-subject variability band for line plots.

    ``sem`` follows the Stage-08 convention: sample SD (ddof=1) divided by the
    square root of the finite subject count at each time point.
    """
    mode = str(mode)
    if mode == "none":
        return None
    values = np.asarray(values_subject_time, dtype=float)
    finite_counts = np.sum(np.isfinite(values), axis=0)
    band = np.full(values.shape[1], np.nan, dtype=float)
    estimable = finite_counts > 1
    if np.any(estimable):
        with np.errstate(invalid="ignore"):
            std = np.nanstd(values[:, estimable], axis=0, ddof=1)
        if mode == "sem":
            band[estimable] = std / np.sqrt(finite_counts[estimable])
        elif mode == "std":
            band[estimable] = std
        else:
            raise ValueError(f"Unsupported band mode: {mode!r}")
    return band


def band_label(mode: str) -> str:
    """Return a compact legend label for a group variability band."""
    if mode == "sem":
        return "+/- SEM"
    if mode == "std":
        return "+/- SD"
    return "no band"


def add_time_markers(
    ax: plt.Axes,
    *,
    occlusion_timing: Mapping[str, Any],
    show_occlusion: bool,
) -> None:
    """Add stimulus onset and optional occlusion/reappearance timing markers."""
    ax.axvline(0.0, color="black", lw=1.0, ls="--", alpha=0.75, label="Stimulus onset")
    if not show_occlusion:
        return
    marker_specs = (
        ("first_complete_occlusion_sec", "First complete occlusion", "#7f7f7f", "-."),
        ("reappearance_sec", "First reappearance", "#8c564b", "--"),
    )
    for key, label, color, linestyle in marker_specs:
        if key not in occlusion_timing:
            continue
        ax.axvline(
            float(occlusion_timing[key]),
            color=color,
            lw=1.3,
            ls=linestyle,
            alpha=0.9,
            label=label,
        )


def contiguous_true_regions(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return inclusive start/end index pairs for contiguous True regions."""
    mask = np.asarray(mask, dtype=bool)
    regions: list[tuple[int, int]] = []
    start: int | None = None
    for idx, value in enumerate(mask):
        if value and start is None:
            start = idx
        elif not value and start is not None:
            regions.append((start, idx - 1))
            start = None
    if start is not None:
        regions.append((start, mask.size - 1))
    return regions


def _valid_topomap_times(times: np.ndarray, requested_times: Iterable[float]) -> list[float]:
    """Return requested topomap times that lie inside the ERF time axis."""
    tmin = float(times[0])
    tmax = float(times[-1])
    valid: list[float] = []
    for value in requested_times:
        time_value = float(value)
        if tmin <= time_value <= tmax:
            valid.append(time_value)
    if not valid:
        valid.append(float(times[int(np.argmin(np.abs(times - 0.0)))]))
    return valid


def plot_butterfly_with_topomaps(
    *,
    data_chan_time: np.ndarray,
    times: np.ndarray,
    reference_info: mne.Info,
    channel_names: list[str],
    sensor_type: str,
    sensor_group: str,
    title: str,
    topomap_times: Iterable[float],
    occlusion_timing: Mapping[str, Any],
    show_occlusion: bool,
    color: str,
    colors: Mapping[str, Any],
    nave: int,
) -> plt.Figure:
    """
    Plot a butterfly/channel time-series with scalp topomaps above it.

    The topomap strip intentionally mirrors a classic ERF/ERP visual layout:
    selected scalp maps provide spatial snapshots while the bottom axis shows
    all selected sensor waveforms over time.
    """
    data = np.asarray(data_chan_time, dtype=float)
    scale = _scale_for_sensor(sensor_type)
    unit = _unit_for_sensor(sensor_type)
    scaled = data * scale
    topo_times = _valid_topomap_times(times, topomap_times)
    n_topos = len(topo_times)
    evoked = make_evoked(
        data,
        reference_info=reference_info,
        times=times,
        comment=title,
        nave=nave,
    )

    fig = plt.figure(figsize=(max(12.0, 2.1 * n_topos), 7.0), constrained_layout=True)
    grid = GridSpec(2, n_topos, height_ratios=[1.0, 2.6], figure=fig)
    vmax = float(np.nanpercentile(np.abs(scaled), 99.0))
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = 1.0

    for idx, topo_time in enumerate(topo_times):
        ax_topo = fig.add_subplot(grid[0, idx])
        evoked.plot_topomap(
            times=[topo_time],
            ch_type=sensor_type,
            scalings={sensor_type: scale},
            units={sensor_type: unit},
            axes=ax_topo,
            colorbar=False,
            vlim=(-vmax, vmax),
            time_format=f"{topo_time:.3f} s",
            show=False,
        )

    ax = fig.add_subplot(grid[1, :])
    if sensor_group_output_tag(sensor_group) == "wholebrain":
        for channel_idx, channel in enumerate(channel_names):
            ax.plot(
                times,
                scaled[channel_idx],
                color=colors.get(channel, "0.35"),
                alpha=0.72,
                lw=0.75,
            )
    else:
        ax.plot(times, scaled.T, color="0.35", alpha=0.22, lw=0.65)
        rms = np.sqrt(np.mean(scaled * scaled, axis=0))
        ax.plot(times, rms, color=color, lw=2.1, label="RMS/GFP across selected sensors")
    inset = ax.inset_axes([0.018, 0.56, 0.17, 0.36])
    draw_sensor_layout(
        inset,
        reference_info=reference_info,
        channel_names=channel_names,
        sensor_type=sensor_type,
        colors=colors,
        title="Sensors",
        show_frame=False,
    )
    add_time_markers(ax, occlusion_timing=occlusion_timing, show_occlusion=show_occlusion)
    ax.axhline(0.0, color="black", lw=0.8, ls=":", alpha=0.7)
    ax.set_title(title)
    ax.set_xlabel("Time from condition onset (s)")
    ax.set_ylabel(f"Amplitude ({unit})")
    if sensor_group_output_tag(sensor_group) != "wholebrain":
        ax.legend(loc="upper right", frameon=True)
    return fig


def plot_subject_condition_overlay(
    *,
    result: Mapping[str, Any],
    sensor_type: str,
    occlusion_timing: Mapping[str, Any],
) -> plt.Figure:
    """Plot one subject's condition RMS/GFP overlay as a diagnostic."""
    erfs = np.asarray(result["erfs"], dtype=float)
    times = np.asarray(result["times"], dtype=float)
    scale = _scale_for_sensor(sensor_type)
    unit = _unit_for_sensor(sensor_type)
    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    for cond_idx, condition in enumerate(CONDITION_ORDER):
        scaled = erfs[cond_idx] * scale
        rms = np.sqrt(np.mean(scaled * scaled, axis=0))
        ax.plot(
            times,
            rms,
            lw=2.0,
            color=CONDITION_COLORS[condition],
            label=CONDITION_LABELS[condition],
        )
    add_time_markers(ax, occlusion_timing=occlusion_timing, show_occlusion=True)
    ax.set_title(f"{result['header']['subject']} | Condition RMS/GFP overlay")
    ax.set_xlabel("Time from condition onset (s)")
    ax.set_ylabel(f"RMS/GFP ({unit})")
    ax.legend(loc="upper right", frameon=True)
    return fig


def plot_group_condition_overlay(
    *,
    erfs_by_subject: np.ndarray,
    times: np.ndarray,
    reference_info: mne.Info,
    channel_names: list[str],
    sensor_type: str,
    sensor_group: str,
    occlusion_timing: Mapping[str, Any],
    band_mode: str,
    colors: Mapping[str, Any],
) -> plt.Figure:
    """Plot group condition RMS/GFP overlays with optional across-subject bands."""
    scale = _scale_for_sensor(sensor_type)
    unit = _unit_for_sensor(sensor_type)
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    condition_to_index = {condition: idx for idx, condition in enumerate(CONDITION_ORDER)}
    for condition in GROUP_OVERLAY_CONDITIONS:
        cond_idx = condition_to_index[condition]
        scaled = erfs_by_subject[:, cond_idx, :, :] * scale
        subject_rms = np.sqrt(np.mean(scaled * scaled, axis=1))
        mean_line = subject_rms.mean(axis=0)
        band = across_subject_band(subject_rms, band_mode)
        if band is not None:
            finite_band = np.isfinite(mean_line) & np.isfinite(band)
            if np.any(finite_band):
                ax.fill_between(
                    times,
                    mean_line - band,
                    mean_line + band,
                    where=finite_band,
                    color=CONDITION_COLORS[condition],
                    alpha=0.18,
                    linewidth=0.0,
                    label="_nolegend_",
                )
        ax.plot(
            times,
            mean_line,
            color=CONDITION_COLORS[condition],
            lw=2.4,
            label=CONDITION_LABELS[condition],
        )
    inset_colors = (
        colors
        if sensor_group_output_tag(sensor_group) == "wholebrain"
        else {channel: "#444444" for channel in channel_names}
    )
    inset = ax.inset_axes([0.018, 0.56, 0.17, 0.36])
    draw_sensor_layout(
        inset,
        reference_info=reference_info,
        channel_names=channel_names,
        sensor_type=sensor_type,
        colors=inset_colors,
        title="Sensors",
        show_frame=False,
    )
    add_time_markers(ax, occlusion_timing=occlusion_timing, show_occlusion=True)
    ax.set_title("Group condition overlays | RMS/GFP across selected sensors")
    ax.set_xlabel("Time from condition onset (s)")
    ax.set_ylabel(f"RMS/GFP ({unit})")
    ax.legend(title=band_label(band_mode), loc="upper right", frameon=True)
    return fig


def plot_trial_counts(trial_counts: pd.DataFrame) -> plt.Figure:
    """Plot retained trial counts for every subject and source condition."""
    pivot = trial_counts.pivot(index="subject", columns="condition", values="n_trials")
    fig, ax = plt.subplots(figsize=(11.5, 5.5), constrained_layout=True)
    x = np.arange(len(pivot.index), dtype=float)
    width = 0.24
    offsets = {
        "always_visible": -width,
        "occluded_nondeviant": 0.0,
        "occluded_deviant": width,
    }
    for condition in SOURCE_CONDITIONS:
        ax.bar(
            x + offsets[condition],
            pivot[condition].to_numpy(dtype=float),
            width=width,
            color=CONDITION_COLORS[condition],
            label=CONDITION_LABELS[condition],
        )
    ax.set_xticks(x)
    ax.set_xticklabels(list(pivot.index))
    ax.set_ylabel("Retained trials")
    ax.set_title("Retained Stage-04 trials by subject and condition")
    ax.legend(loc="upper right", frameon=True)
    return fig


def plot_contrast_butterfly(
    *,
    contrast_values: np.ndarray,
    times: np.ndarray,
    sensor_type: str,
    title: str,
    occlusion_timing: Mapping[str, Any],
) -> plt.Figure:
    """Plot group mean contrast channel waveforms plus contrast RMS."""
    group_mean = np.mean(contrast_values, axis=0)
    scale = _scale_for_sensor(sensor_type)
    unit = _unit_for_sensor(sensor_type)
    scaled = group_mean * scale
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    ax.plot(times, scaled.T, color="0.35", alpha=0.22, lw=0.65)
    rms = np.sqrt(np.mean(scaled * scaled, axis=0))
    ax.plot(times, rms, color="#d62728", lw=2.2, label="Contrast RMS/GFP")
    add_time_markers(ax, occlusion_timing=occlusion_timing, show_occlusion=True)
    ax.axhline(0.0, color="black", lw=0.8, ls=":", alpha=0.7)
    ax.set_title(title)
    ax.set_xlabel("Time from condition onset (s)")
    ax.set_ylabel(f"Contrast amplitude ({unit})")
    ax.legend(loc="upper right", frameon=True)
    return fig


def plot_cluster_highlight(
    *,
    cluster_result: Mapping[str, Any],
    contrast_values: np.ndarray,
    times: np.ndarray,
    channel_names: list[str],
    sensor_type: str,
    title: str,
    occlusion_timing: Mapping[str, Any],
) -> plt.Figure:
    """
    Plot t-statistic heatmap and contrast butterfly with significant clusters.

    Significant clusters are projected to a time mask in the lower panel. The
    heatmap keeps the full sensor x time result visible for auditability.
    """
    t_obs = np.asarray(cluster_result["t_obs"], dtype=float)
    sig_mask = np.asarray(cluster_result["significant_mask"], dtype=bool)
    group_mean = np.mean(contrast_values, axis=0)
    scale = _scale_for_sensor(sensor_type)
    unit = _unit_for_sensor(sensor_type)
    scaled = group_mean * scale

    fig, (ax_heat, ax_wave) = plt.subplots(
        2,
        1,
        figsize=(13, 8.5),
        sharex=True,
        gridspec_kw={"height_ratios": [1.3, 1.0]},
        constrained_layout=True,
    )
    vmax = float(np.nanpercentile(np.abs(t_obs), 99.0))
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = 1.0
    image = ax_heat.imshow(
        t_obs.T,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        extent=[float(times[0]), float(times[-1]), -0.5, len(channel_names) - 0.5],
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    if np.any(sig_mask):
        ax_heat.contour(
            times,
            np.arange(len(channel_names), dtype=float),
            sig_mask.T.astype(float),
            levels=[0.5],
            colors="black",
            linewidths=0.8,
        )
    cbar = fig.colorbar(image, ax=ax_heat, pad=0.01)
    cbar.set_label("t value")
    ax_heat.set_ylabel("Sensor index")
    ax_heat.set_title(title + " | MNE spatio-temporal cluster t map")
    add_time_markers(ax_heat, occlusion_timing=occlusion_timing, show_occlusion=True)

    ax_wave.plot(times, scaled.T, color="0.35", alpha=0.2, lw=0.65)
    rms = np.sqrt(np.mean(scaled * scaled, axis=0))
    ax_wave.plot(times, rms, color="#d62728", lw=2.2, label="Group contrast RMS/GFP")
    for start_idx, end_idx in contiguous_true_regions(sig_mask.any(axis=1)):
        ax_wave.axvspan(
            float(times[start_idx]),
            float(times[end_idx]),
            color="#ffcc66",
            alpha=0.35,
            lw=0,
        )
    add_time_markers(ax_wave, occlusion_timing=occlusion_timing, show_occlusion=True)
    ax_wave.axhline(0.0, color="black", lw=0.8, ls=":", alpha=0.7)
    ax_wave.set_xlabel("Time from condition onset (s)")
    ax_wave.set_ylabel(f"Contrast amplitude ({unit})")
    ax_wave.legend(loc="upper right", frameon=True)
    return fig


def plot_channel_average_cluster_highlight(
    *,
    cluster_result: Mapping[str, Any],
    values_time: np.ndarray,
    times: np.ndarray,
    sensor_type: str,
    title: str,
    occlusion_timing: Mapping[str, Any],
) -> plt.Figure:
    """
    Plot the signed selected-channel average contrast and 1-D cluster spans.

    The line is the across-subject mean of each subject's arithmetic channel
    average. The shaded uncertainty band is SEM across subjects; yellow spans
    mark significant 1-D temporal clusters.
    """
    values_time = np.asarray(values_time, dtype=float)
    scale = _scale_for_sensor(sensor_type)
    unit = _unit_for_sensor(sensor_type)
    scaled = values_time * scale
    mean_line = scaled.mean(axis=0)
    sem = across_subject_band(scaled, "sem")
    sig_mask = np.asarray(cluster_result["significant_mask"], dtype=bool)

    fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)
    if sem is not None:
        finite_band = np.isfinite(mean_line) & np.isfinite(sem)
        if np.any(finite_band):
            ax.fill_between(
                times,
                mean_line - sem,
                mean_line + sem,
                where=finite_band,
                color="#d62728",
                alpha=0.18,
                linewidth=0.0,
            )
    ax.plot(
        times,
        mean_line,
        color="#d62728",
        lw=2.4,
        label="Selected-channel average contrast (+/- SEM)",
    )
    for start_idx, end_idx in contiguous_true_regions(sig_mask):
        ax.axvspan(
            float(times[start_idx]),
            float(times[end_idx]),
            color="#ffcc66",
            alpha=0.35,
            lw=0,
        )
    add_time_markers(ax, occlusion_timing=occlusion_timing, show_occlusion=True)
    ax.axhline(0.0, color="black", lw=0.8, ls=":", alpha=0.7)
    ax.set_title(title + " | 1-D channel-average cluster test")
    ax.set_xlabel("Time from condition onset (s)")
    ax.set_ylabel(f"Channel-average contrast ({unit})")
    ax.legend(loc="upper right", frameon=True)
    return fig


def save_figure(fig: plt.Figure, path: Path) -> None:
    """Save and close a Matplotlib figure."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# HTML report helpers
# -----------------------------------------------------------------------------


def _jsonify(value: Any) -> Any:
    """Convert NumPy-heavy objects to JSON-serializable structures."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _jsonify(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonify(item) for item in value]
    return value


def _figure_html(
    report_path: Path,
    *,
    image_path: Path,
    caption: str,
    embed: bool,
) -> str:
    """Return figure HTML, optionally embedding the image as base64."""
    if embed:
        mime_type = mimetypes.guess_type(image_path.name)[0] or "image/png"
        encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
        src = f"data:{mime_type};base64,{encoded}"
    else:
        src = os.path.relpath(image_path, start=report_path.parent)
    return (
        "<figure>"
        f"<img src=\"{html.escape(src)}\" alt=\"{html.escape(image_path.name)}\">"
        f"<figcaption>{html.escape(caption)}</figcaption>"
        "</figure>"
    )


def _records_table(records: Iterable[Mapping[str, Any]], columns: list[str]) -> str:
    """Render selected record columns as an HTML table."""
    rows = list(records)
    header = "".join(f"<th>{html.escape(col)}</th>" for col in columns)
    body_rows: list[str] = []
    for row in rows:
        cells = "".join(
            f"<td>{html.escape(str(row.get(col, '')))}</td>" for col in columns
        )
        body_rows.append(f"<tr>{cells}</tr>")
    return (
        "<table><thead><tr>"
        + header
        + "</tr></thead><tbody>"
        + "\n".join(body_rows)
        + "</tbody></table>"
    )


def _cluster_table(cluster_results: Mapping[str, Mapping[str, Any]]) -> str:
    """Render concise cluster rows for every contrast."""
    rows: list[dict[str, Any]] = []
    for contrast_key, result in cluster_results.items():
        clusters = list(result.get("clusters", []))
        if not clusters:
            rows.append(
                {
                    "contrast": contrast_key,
                    "cluster": "none",
                    "p_value": "",
                    "significant": "",
                    "time_start_sec": "",
                    "time_end_sec": "",
                    "peak_time_sec": "",
                    "peak_sensor": "",
                    "peak_t": "",
                }
            )
            continue
        for cluster in clusters:
            rows.append(
                {
                    "contrast": contrast_key,
                    "cluster": cluster.get("cluster_index"),
                    "p_value": f"{float(cluster.get('p_value', np.nan)):.4f}",
                    "significant": cluster.get("significant"),
                    "time_start_sec": f"{float(cluster.get('time_start_sec', np.nan)):.3f}",
                    "time_end_sec": f"{float(cluster.get('time_end_sec', np.nan)):.3f}",
                    "peak_time_sec": f"{float(cluster.get('peak_time_sec', np.nan)):.3f}",
                    "peak_sensor": cluster.get("peak_sensor"),
                    "peak_t": f"{float(cluster.get('peak_t', np.nan)):.3f}",
                }
            )
    return _records_table(
        rows,
        [
            "contrast",
            "cluster",
            "p_value",
            "significant",
            "time_start_sec",
            "time_end_sec",
            "peak_time_sec",
            "peak_sensor",
            "peak_t",
        ],
    )


def write_html_report(
    path: Path,
    *,
    title: str,
    summary: Mapping[str, Any],
    trial_counts: pd.DataFrame,
    imbalance: pd.DataFrame,
    cluster_results: Mapping[str, Mapping[str, Any]],
    figure_paths: list[tuple[Path, str]],
    skipped: list[Mapping[str, Any]],
    embed_images: bool,
    include_cluster_table: bool,
) -> Path:
    """Write the linked or fully embedded group HTML report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    method_rows = [
        ("Analysis status", "WORK_IN_PROGRESS"),
        ("Baseline interval", summary["baseline_interval_sec"]),
        ("Time window", summary["time_window_sec"]),
        ("Sensor selection", summary["sensor_selection"]),
        ("Units/scaling", summary["units_and_scaling"]),
        ("Occluded combined", summary["occluded_combined_construction"]),
        ("Contrast sign convention", summary["contrast_sign_convention"]),
        ("Spatio-temporal cluster parameters", summary["cluster_parameters"]),
        ("1-D channel-average cluster parameters", summary["cluster1d_parameters"]),
        ("Group overlay", summary["group_overlay"]),
        ("Topomap times sec", summary["topomap_times_sec"]),
        ("Imbalance policy", summary["imbalance_policy"]),
    ]
    method_table = _records_table(
        [{"setting": key, "value": json.dumps(_jsonify(value), sort_keys=True)} for key, value in method_rows],
        ["setting", "value"],
    )
    trial_table = _records_table(
        trial_counts[
            [
                "subject",
                "condition",
                "n_trials",
                "n_identities",
                "min_trials_per_identity",
                "max_trials_per_identity",
                "median_trials_per_identity",
            ]
        ].to_dict("records"),
        [
            "subject",
            "condition",
            "n_trials",
            "n_identities",
            "min_trials_per_identity",
            "max_trials_per_identity",
            "median_trials_per_identity",
        ],
    )
    imbalance_table = _records_table(
        imbalance.to_dict("records"),
        [
            "subject",
            "contrast",
            "majority_to_minority_ratio",
            "threshold",
            "trigger_equalized_sensitivity",
        ],
    )
    skipped_table = (
        _records_table(skipped, ["subject", "reason"])
        if skipped
        else "<p>No subjects or conditions were skipped.</p>"
    )
    figures_html = "\n".join(
        _figure_html(path, image_path=fig_path, caption=caption, embed=embed_images)
        for fig_path, caption in figure_paths
    )

    cluster_table_section = (
        ["<h2>Cluster Results</h2>", _cluster_table(cluster_results)]
        if include_cluster_table
        else [
            "<h2>Cluster Results</h2>",
            (
                "<p>Verbose cluster table omitted by default. Use "
                "<code>--include-cluster-table</code> to include it in the report; "
                "full cluster rows are saved in JSON outputs.</p>"
            ),
        ]
    )
    machine_summary = _jsonify(summary)
    if not include_cluster_table:
        for section in (
            "clusters",
            "channel_average_clusters",
            "equalized_sensitivity_clusters",
            "equalized_sensitivity_channel_average_clusters",
        ):
            section_value = machine_summary.get(section, {})
            if isinstance(section_value, dict):
                for cluster_summary in section_value.values():
                    if isinstance(cluster_summary, dict) and "clusters" in cluster_summary:
                        cluster_summary["clusters"] = (
                            "omitted from HTML by default; see JSON cluster outputs"
                        )
    body_parts = [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\">",
            f"<title>{html.escape(title)}</title>",
            "<style>",
            "body{font-family:system-ui,sans-serif;max-width:1320px;margin:2rem auto;padding:0 1rem;line-height:1.45;}",
            ".banner{background:#fff3cd;border:2px solid #d39e00;padding:1rem;border-radius:6px;font-weight:700;}",
            "figure{margin:1.6rem 0;padding-bottom:1rem;border-bottom:1px solid #ddd;}",
            "img{max-width:100%;height:auto;border:1px solid #ddd;}",
            "figcaption{font-size:.95rem;color:#333;margin-top:.45rem;}",
            "table{border-collapse:collapse;margin:1rem 0;width:100%;font-size:.92rem;}",
            "td,th{border:1px solid #ddd;padding:.4rem;text-align:left;vertical-align:top;}",
            "th{background:#f4f4f4;}",
            "pre{background:#f6f8fa;padding:1rem;border-radius:6px;overflow:auto;}",
            "h2{margin-top:2.4rem;border-top:3px solid #222;padding-top:1rem;}",
            "</style>",
            "</head><body>",
            f"<h1>{html.escape(title)}</h1>",
            f"<div class=\"banner\">{html.escape(WORK_IN_PROGRESS_TEXT)}</div>",
            "<h2>Method Choices</h2>",
            method_table,
            "<h2>Trial Counts</h2>",
            trial_table,
            "<h2>Imbalance Diagnostics</h2>",
            imbalance_table,
            *cluster_table_section,
            "<h2>Skipped Subjects / Conditions</h2>",
            skipped_table,
            "<h2>Figures</h2>",
            figures_html,
            "<h2>Machine Summary</h2>",
            "<pre>" + html.escape(json.dumps(machine_summary, indent=2, sort_keys=True)) + "</pre>",
            "</body></html>",
    ]
    body = "\n".join(body_parts)
    path.write_text(body, encoding="utf-8")
    return path


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------


def write_cluster_outputs(
    *,
    paths: Mapping[str, Path],
    contrast_key: str,
    result: Mapping[str, Any],
    times: np.ndarray,
    channel_names: list[str],
    args: argparse.Namespace,
    suffix: str = "",
) -> None:
    """Write NPZ and JSON files for one cluster result."""
    if suffix:
        npz_path = paths[f"cluster_{contrast_key}_{suffix}_npz"]
        json_path = paths[f"cluster_{contrast_key}_{suffix}_json"]
    else:
        npz_path = paths[f"cluster_{contrast_key}_npz"]
        json_path = paths[f"cluster_{contrast_key}_json"]
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        npz_path,
        times_sec=times,
        channel_names=np.asarray(channel_names, dtype=object),
        t_obs=np.asarray(result["t_obs"], dtype=float),
        significant_mask=np.asarray(result["significant_mask"], dtype=bool),
        cluster_p_values=np.asarray(result["cluster_p_values"], dtype=float),
        h0=np.asarray(result["h0"], dtype=float),
        threshold=np.asarray(result["threshold"], dtype=float),
    )
    payload = {
        "contrast": contrast_key,
        "clusters": result["clusters"],
        "cluster_parameters": cluster_parameter_summary(args, result),
        "outputs": {"npz": str(npz_path), "json": str(json_path)},
    }
    json_path.write_text(json.dumps(_jsonify(payload), indent=2), encoding="utf-8")


def write_cluster1d_outputs(
    *,
    paths: Mapping[str, Path],
    contrast_key: str,
    result: Mapping[str, Any],
    times: np.ndarray,
    args: argparse.Namespace,
    suffix: str = "",
) -> None:
    """Write NPZ and JSON files for one channel-average temporal cluster result."""
    if suffix:
        npz_path = paths[f"cluster1d_{contrast_key}_{suffix}_npz"]
        json_path = paths[f"cluster1d_{contrast_key}_{suffix}_json"]
    else:
        npz_path = paths[f"cluster1d_{contrast_key}_npz"]
        json_path = paths[f"cluster1d_{contrast_key}_json"]
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        npz_path,
        times_sec=times,
        t_obs=np.asarray(result["t_obs"], dtype=float),
        significant_mask=np.asarray(result["significant_mask"], dtype=bool),
        cluster_p_values=np.asarray(result["cluster_p_values"], dtype=float),
        h0=np.asarray(result["h0"], dtype=float),
        threshold=np.asarray(result["threshold"], dtype=float),
    )
    payload = {
        "contrast": contrast_key,
        "clusters": result["clusters"],
        "cluster_parameters": cluster1d_parameter_summary(args, result),
        "outputs": {"npz": str(npz_path), "json": str(json_path)},
    }
    json_path.write_text(json.dumps(_jsonify(payload), indent=2), encoding="utf-8")


def cluster_parameter_summary(
    args: argparse.Namespace,
    result: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return report-ready cluster permutation parameters."""
    summary = {
        "statistic": "one-sample t-test of paired subject-level difference ERFs against zero",
        "mne_function": "mne.stats.spatio_temporal_cluster_1samp_test",
        "tail": args.tail,
        "tail_mne": mne_tail_value(args.tail),
        "cluster_forming_p": float(args.cluster_forming_p),
        "cluster_alpha": float(args.cluster_alpha),
        "n_permutations": (
            str(args.n_permutations)
            if isinstance(args.n_permutations, str)
            else int(args.n_permutations)
        ),
        "seed": int(args.random_seed),
        "adjacency": "mne.channels.find_ch_adjacency on selected common-good MEG info",
        "max_step": 1,
        "out_type": "indices",
    }
    if result is not None:
        summary["t_threshold"] = float(result["threshold"])
        summary["n_adjacency_channels"] = int(len(result["adjacency_names"]))
    return summary


def cluster1d_parameter_summary(
    args: argparse.Namespace,
    result: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return report-ready parameters for the 1-D channel-average cluster test."""
    summary = {
        "statistic": (
            "one-sample t-test of paired subject-level difference ERFs after "
            "signed arithmetic averaging over selected channels"
        ),
        "mne_function": "mne.stats.permutation_cluster_1samp_test",
        "tail": args.tail,
        "tail_mne": mne_tail_value(args.tail),
        "cluster_forming_p": float(args.cluster_forming_p),
        "cluster_alpha": float(args.cluster_alpha),
        "n_permutations": (
            str(args.n_permutations)
            if isinstance(args.n_permutations, str)
            else int(args.n_permutations)
        ),
        "seed": int(args.random_seed),
        "adjacency": "1-D temporal adjacency implicit in MNE permutation_cluster_1samp_test",
        "out_type": "indices",
    }
    if result is not None:
        summary.update(
            {
                "threshold_t": float(result["threshold"]),
                "n_clusters": int(len(result["clusters"])),
                "n_significant_clusters": int(
                    sum(1 for cluster in result["clusters"] if cluster.get("significant"))
                ),
            }
        )
    return summary


def build_group_outputs(
    *,
    args: argparse.Namespace,
    shared_cfg: Mapping[str, Any],
    subject_headers: list[dict[str, Any]],
    skipped: list[dict[str, Any]],
    channel_names: list[str],
    sensor_group: str,
    sensor_selection: Mapping[str, Any],
    reference_info: mne.Info,
    need_equalized_sensitivity: bool,
) -> dict[str, Any]:
    """Compute subject/group ERFs, cluster tests, figures, and reports."""
    group_paths = group_output_paths(shared_cfg, args.sensor_type, sensor_group)
    raw_root = Path(shared_cfg["resolved_paths"]["data_root"])
    assert_output_path_outside_raw(group_paths["out_dir"], raw_root)
    assert_output_path_outside_raw(group_paths["report_dir"], raw_root)

    subject_results: list[dict[str, Any]] = []
    for header in subject_headers:
        if args.resume and not need_equalized_sensitivity:
            try:
                result = load_subject_erfs_from_outputs(
                    header,
                    args=args,
                    channel_names=channel_names,
                    sensor_group=sensor_group,
                )
            except Exception:
                result = compute_subject_erfs(
                    header,
                    args=args,
                    channel_names=channel_names,
                    sensor_group=sensor_group,
                    sensor_selection=sensor_selection,
                    keep_trial_data=need_equalized_sensitivity,
                )
        else:
            result = compute_subject_erfs(
                header,
                args=args,
                channel_names=channel_names,
                sensor_group=sensor_group,
                sensor_selection=sensor_selection,
                keep_trial_data=need_equalized_sensitivity,
            )
        subject_results.append(result)

    times = np.asarray(subject_results[0]["times"], dtype=float)
    for result in subject_results[1:]:
        current_times = np.asarray(result["times"], dtype=float)
        if len(current_times) != len(times) or not np.allclose(current_times, times):
            raise RuntimeError(f"{result['header']['subject']} ERF time axis differs from group reference.")

    erfs_by_subject = np.stack(
        [np.asarray(result["erfs"], dtype=float) for result in subject_results],
        axis=0,
    )
    contrast_values = build_contrast_values(erfs_by_subject)
    channel_average_contrasts = {
        key: channel_average_values(value) for key, value in contrast_values.items()
    }
    occlusion_timing = occlusion_timing_from_config(shared_cfg)
    line_colors = channel_colors(channel_names)
    trial_counts_df = pd.DataFrame(
        [row for result in subject_results for row in result["trial_count_rows"]]
    )
    imbalance_df = imbalance_diagnostics(
        subject_headers,
        threshold=float(args.imbalance_ratio_threshold),
    )

    group_paths["out_dir"].mkdir(parents=True, exist_ok=True)
    group_paths["report_dir"].mkdir(parents=True, exist_ok=True)
    trial_counts_df.to_csv(group_paths["trial_counts_csv"], index=False)
    imbalance_df.to_csv(group_paths["imbalance_csv"], index=False)

    np.savez_compressed(
        group_paths["group_erfs_npz"],
        times_sec=times,
        subjects=np.asarray([header["subject"] for header in subject_headers], dtype=object),
        channel_names=np.asarray(channel_names, dtype=object),
        condition_names=np.asarray(CONDITION_ORDER, dtype=object),
        erfs_si=erfs_by_subject,
        group_mean_erfs_si=np.mean(erfs_by_subject, axis=0),
        sensor_type=np.asarray(args.sensor_type),
        sensor_group=np.asarray(sensor_group_output_tag(sensor_group)),
        si_unit=np.asarray(SENSOR_DISPLAY[args.sensor_type]["si_unit"]),
        baseline_sec=np.asarray(args.baseline, dtype=float),
    )
    np.savez_compressed(
        group_paths["group_contrasts_npz"],
        times_sec=times,
        subjects=np.asarray([header["subject"] for header in subject_headers], dtype=object),
        channel_names=np.asarray(channel_names, dtype=object),
        contrast_names=np.asarray([spec["key"] for spec in CONTRAST_SPECS], dtype=object),
        contrasts_si=np.stack(
            [contrast_values[str(spec["key"])] for spec in CONTRAST_SPECS],
            axis=1,
        ),
        group_mean_contrasts_si=np.stack(
            [np.mean(contrast_values[str(spec["key"])], axis=0) for spec in CONTRAST_SPECS],
            axis=0,
        ),
        sensor_type=np.asarray(args.sensor_type),
        sensor_group=np.asarray(sensor_group_output_tag(sensor_group)),
        si_unit=np.asarray(SENSOR_DISPLAY[args.sensor_type]["si_unit"]),
    )
    np.savez_compressed(
        group_paths["group_channel_average_contrasts_npz"],
        times_sec=times,
        subjects=np.asarray([header["subject"] for header in subject_headers], dtype=object),
        contrast_names=np.asarray([spec["key"] for spec in CONTRAST_SPECS], dtype=object),
        channel_average_contrasts_si=np.stack(
            [channel_average_contrasts[str(spec["key"])] for spec in CONTRAST_SPECS],
            axis=1,
        ),
        group_mean_channel_average_contrasts_si=np.stack(
            [
                np.mean(channel_average_contrasts[str(spec["key"])], axis=0)
                for spec in CONTRAST_SPECS
            ],
            axis=0,
        ),
        sensor_type=np.asarray(args.sensor_type),
        sensor_group=np.asarray(sensor_group_output_tag(sensor_group)),
        channel_names=np.asarray(channel_names, dtype=object),
        channel_average_policy=np.asarray(
            "signed arithmetic mean over selected channels before 1-D inference"
        ),
        si_unit=np.asarray(SENSOR_DISPLAY[args.sensor_type]["si_unit"]),
    )

    cluster_results: dict[str, dict[str, Any]] = {}
    cluster1d_results: dict[str, dict[str, Any]] = {}
    equalized_cluster_results: dict[str, dict[str, Any]] = {}
    equalized_cluster1d_results: dict[str, dict[str, Any]] = {}
    for spec in CONTRAST_SPECS:
        key = str(spec["key"])
        cluster_result = run_mne_cluster_1samp(
            contrast_values[key],
            reference_info=reference_info,
            times=times,
            args=args,
        )
        cluster_results[key] = cluster_result
        write_cluster_outputs(
            paths=group_paths,
            contrast_key=key,
            result=cluster_result,
            times=times,
            channel_names=channel_names,
            args=args,
        )
        cluster1d_result = run_mne_cluster_1d(
            channel_average_contrasts[key],
            times=times,
            args=args,
        )
        cluster1d_results[key] = cluster1d_result
        write_cluster1d_outputs(
            paths=group_paths,
            contrast_key=key,
            result=cluster1d_result,
            times=times,
            args=args,
        )
        if need_equalized_sensitivity:
            equalized_values = equalized_contrast_values(
                subject_results,
                contrast_key=key,
                repeats=int(args.equalization_repeats),
                random_seed=int(args.random_seed),
            )
            equalized_result = run_mne_cluster_1samp(
                equalized_values,
                reference_info=reference_info,
                times=times,
                args=args,
            )
            equalized_cluster_results[key] = equalized_result
            write_cluster_outputs(
                paths=group_paths,
                contrast_key=key,
                result=equalized_result,
                times=times,
                channel_names=channel_names,
                args=args,
                suffix="equalized",
            )
            equalized_1d_values = channel_average_values(equalized_values)
            equalized_1d_result = run_mne_cluster_1d(
                equalized_1d_values,
                times=times,
                args=args,
            )
            equalized_cluster1d_results[key] = equalized_1d_result
            write_cluster1d_outputs(
                paths=group_paths,
                contrast_key=key,
                result=equalized_1d_result,
                times=times,
                args=args,
                suffix="equalized",
            )

    figure_paths: list[tuple[Path, str]] = []
    fig = plot_sensor_selection_layout(
        reference_info=reference_info,
        channel_names=channel_names,
        sensor_type=args.sensor_type,
        sensor_group=sensor_group,
        colors=line_colors,
    )
    save_figure(fig, group_paths["sensor_layout_figure"])
    figure_paths.append(
        (
            group_paths["sensor_layout_figure"],
            (
                f"Selected Stage-07 {sensor_group_output_tag(sensor_group)} sensors "
                f"after {args.sensor_type} filtering and common-good intersection."
            ),
        )
    )

    group_mean = erfs_by_subject.mean(axis=0)
    for cond_idx, condition in enumerate(CONDITION_ORDER):
        fig_path = group_paths[f"condition_{condition}_butterfly_topomap"]
        fig = plot_butterfly_with_topomaps(
            data_chan_time=group_mean[cond_idx],
            times=times,
            reference_info=reference_info,
            channel_names=channel_names,
            sensor_type=args.sensor_type,
            sensor_group=sensor_group,
            title=f"Group ERF | {CONDITION_LABELS[condition]}",
            topomap_times=args.topomap_times,
            occlusion_timing=occlusion_timing,
            show_occlusion="occluded" in condition,
            color=CONDITION_COLORS[condition],
            colors=(
                line_colors
                if sensor_group_output_tag(sensor_group) == "wholebrain"
                else {channel: CONDITION_COLORS[condition] for channel in channel_names}
            ),
            nave=len(subject_headers),
        )
        save_figure(fig, fig_path)
        figure_paths.append(
            (
                fig_path,
                (
                    f"Group mean {CONDITION_LABELS[condition]} ERF. Butterfly traces show "
                    f"selected Stage-07 {sensor_group_output_tag(sensor_group)} "
                    f"{SENSOR_DISPLAY[args.sensor_type]['label']} in "
                    f"{SENSOR_DISPLAY[args.sensor_type]['plot_unit']}; topomaps use "
                    f"selected latencies; baseline {args.baseline[0]}..{args.baseline[1]} s."
                ),
            )
        )

    fig = plot_group_condition_overlay(
        erfs_by_subject=erfs_by_subject,
        times=times,
        reference_info=reference_info,
        channel_names=channel_names,
        sensor_type=args.sensor_type,
        sensor_group=sensor_group,
        occlusion_timing=occlusion_timing,
        band_mode=args.group_overlay_band,
        colors=line_colors,
    )
    save_figure(fig, group_paths["condition_overlay_figure"])
    figure_paths.append(
        (
            group_paths["condition_overlay_figure"],
            (
                "Condition RMS/GFP overlay. Individual subject traces are omitted. "
                "Lines are group means for source "
                f"conditions only with {band_label(args.group_overlay_band)} shading. "
                "Occlusion markers show first complete occlusion and first reappearance."
            ),
        )
    )

    fig = plot_trial_counts(trial_counts_df)
    save_figure(fig, group_paths["trial_counts_figure"])
    figure_paths.append(
        (
            group_paths["trial_counts_figure"],
            "Retained Stage-04 trials by subject and source condition.",
        )
    )

    for spec in CONTRAST_SPECS:
        key = str(spec["key"])
        fig = plot_contrast_butterfly(
            contrast_values=contrast_values[key],
            times=times,
            sensor_type=args.sensor_type,
            title=f"Group contrast | {spec['label']}",
            occlusion_timing=occlusion_timing,
        )
        save_figure(fig, group_paths[f"contrast_{key}_butterfly"])
        figure_paths.append(
            (
                group_paths[f"contrast_{key}_butterfly"],
                (
                    f"Group paired contrast ERF for {spec['label']}. Positive values "
                    "follow the displayed sign convention."
                ),
            )
        )
        fig = plot_cluster_highlight(
            cluster_result=cluster_results[key],
            contrast_values=contrast_values[key],
            times=times,
            channel_names=channel_names,
            sensor_type=args.sensor_type,
            title=f"Cluster test | {spec['label']}",
            occlusion_timing=occlusion_timing,
        )
        save_figure(fig, group_paths[f"cluster_{key}_figure"])
        figure_paths.append(
            (
                group_paths[f"cluster_{key}_figure"],
                (
                    f"MNE spatio-temporal cluster result for {spec['label']}. Heatmap "
                    "shows t values over time x sensor; black contours and shaded "
                    "time spans mark clusters with p <= cluster alpha."
                ),
            )
        )
        fig = plot_channel_average_cluster_highlight(
            cluster_result=cluster1d_results[key],
            values_time=channel_average_contrasts[key],
            times=times,
            sensor_type=args.sensor_type,
            title=f"Channel-average cluster test | {spec['label']}",
            occlusion_timing=occlusion_timing,
        )
        save_figure(fig, group_paths[f"cluster1d_{key}_figure"])
        figure_paths.append(
            (
                group_paths[f"cluster1d_{key}_figure"],
                (
                    f"1-D temporal cluster result for {spec['label']} after signed "
                    "arithmetic averaging over selected channels within each subject. "
                    "Line is the group mean channel-average contrast; band is SEM."
                ),
            )
        )

    # Subject diagnostics are generated after group figures so the report reads
    # from group-level results down to participant-level checks.
    for result in subject_results:
        fig = plot_subject_condition_overlay(
            result=result,
            sensor_type=args.sensor_type,
            occlusion_timing=occlusion_timing,
        )
        save_figure(fig, result["paths"]["condition_overlay_figure"])
        figure_paths.append(
            (
                result["paths"]["condition_overlay_figure"],
                f"{result['header']['subject']} condition RMS/GFP diagnostic.",
            )
        )
        fig = plot_subject_condition_overlay(
            result=result,
            sensor_type=args.sensor_type,
            occlusion_timing=occlusion_timing,
        )
        save_figure(fig, result["paths"]["gfp_rms_figure"])

    summary = {
        "stage": "dot_erf_group",
        "analysis_status": "WORK_IN_PROGRESS",
        "work_in_progress_warning": WORK_IN_PROGRESS_TEXT,
        "sensor_group": sensor_group_output_tag(sensor_group),
        "subjects": [header["subject"] for header in subject_headers],
        "requested_sensor_groups": list(args.sensor_groups),
        "skipped_subjects_or_conditions": skipped,
        "conditions": list(CONDITION_ORDER),
        "source_conditions": list(SOURCE_CONDITIONS),
        "baseline_interval_sec": [float(args.baseline[0]), float(args.baseline[1])],
        "time_window_sec": [float(times[0]), float(times[-1])],
        "analysis_sfreq_hz": float(1.0 / np.median(np.diff(times))),
        "sensor_selection": {
            "sensor_type": args.sensor_type,
            "sensor_type_label": SENSOR_DISPLAY[args.sensor_type]["label"],
            "sensor_group": sensor_group_output_tag(sensor_group),
            "policy": sensor_selection["sensor_group_policy"],
            "n_common_good_channels": int(len(channel_names)),
            "common_good_channels": list(channel_names),
            "stage07_selection": dict(sensor_selection),
        },
        "units_and_scaling": {
            "si_unit": SENSOR_DISPLAY[args.sensor_type]["si_unit"],
            "plot_unit": SENSOR_DISPLAY[args.sensor_type]["plot_unit"],
            "plot_scale_from_si": float(SENSOR_DISPLAY[args.sensor_type]["scale"]),
        },
        "occluded_combined_construction": (
            "Per subject, occluded_combined = 0.5 * ERF(occluded_nondeviant) "
            "+ 0.5 * ERF(occluded_deviant); source trials are not pooled."
        ),
        "contrast_sign_convention": {
            str(spec["key"]): f"{spec['positive_condition']} - {spec['negative_condition']}"
            for spec in CONTRAST_SPECS
        },
        "cluster_parameters": cluster_parameter_summary(args, next(iter(cluster_results.values()))),
        "cluster1d_parameters": cluster1d_parameter_summary(
            args, next(iter(cluster1d_results.values()))
        ),
        "group_overlay": {
            "conditions_plotted": list(GROUP_OVERLAY_CONDITIONS),
            "excluded_conditions": ["occluded_combined"],
            "band_mode": args.group_overlay_band,
            "band_interpretation": band_label(args.group_overlay_band),
            "individual_subject_lines": False,
        },
        "topomap_times_sec": [float(value) for value in args.topomap_times],
        "channel_average_statistic": (
            "For each contrast, each subject's signed contrast ERF is averaged "
            "arithmetically over selected channels before the 1-D temporal "
            "cluster test against zero."
        ),
        "imbalance_policy": {
            "primary_policy": "Use all retained source trials; construct occluded_combined as equal-weight subtype average.",
            "equalization_trigger_threshold": float(args.imbalance_ratio_threshold),
            "threshold_justification": (
                "A 1.25 majority/minority source-count ratio corresponds to about a "
                "12% difference in expected averaging-noise standard deviation under "
                "sqrt(n) scaling; current primary analysis avoids unnecessary "
                "subsampling unless imbalances exceed this pragmatic WIP threshold."
            ),
            "equalized_sensitivity_triggered": bool(need_equalized_sensitivity),
            "equalization_repeats": int(args.equalization_repeats),
            "equalized_sensitivity_results": (
                "Written and included in JSON outputs."
                if need_equalized_sensitivity
                else "Not triggered by inspected trial-count ratios."
            ),
        },
        "occlusion_timing": occlusion_timing,
        "outputs": {key: str(value) for key, value in group_paths.items() if key not in {"out_dir", "report_dir"}},
        "subject_outputs": [
            {key: str(value) for key, value in result["paths"].items() if key != "out_dir"}
            for result in subject_results
        ],
        "clusters": {
            key: {
                "clusters": result["clusters"],
                "n_significant_clusters": int(
                    sum(1 for cluster in result["clusters"] if cluster.get("significant"))
                ),
            }
            for key, result in cluster_results.items()
        },
        "channel_average_clusters": {
            key: {
                "clusters": result["clusters"],
                "n_significant_clusters": int(
                    sum(1 for cluster in result["clusters"] if cluster.get("significant"))
                ),
            }
            for key, result in cluster1d_results.items()
        },
        "equalized_sensitivity_clusters": {
            key: {
                "clusters": result["clusters"],
                "n_significant_clusters": int(
                    sum(1 for cluster in result["clusters"] if cluster.get("significant"))
                ),
            }
            for key, result in equalized_cluster_results.items()
        },
        "equalized_sensitivity_channel_average_clusters": {
            key: {
                "clusters": result["clusters"],
                "n_significant_clusters": int(
                    sum(1 for cluster in result["clusters"] if cluster.get("significant"))
                ),
            }
            for key, result in equalized_cluster1d_results.items()
        },
    }
    group_paths["group_summary_json"].write_text(
        json.dumps(_jsonify(summary), indent=2),
        encoding="utf-8",
    )
    write_html_report(
        group_paths["report_html"],
        title=(
            "WORK_IN_PROGRESS Dot ERF Sanity Check | "
            f"{SENSOR_DISPLAY[args.sensor_type]['label']} | "
            f"{sensor_group_output_tag(sensor_group)}"
        ),
        summary=summary,
        trial_counts=trial_counts_df,
        imbalance=imbalance_df,
        cluster_results=cluster_results,
        figure_paths=figure_paths,
        skipped=skipped,
        embed_images=False,
        include_cluster_table=bool(args.include_cluster_table),
    )
    write_html_report(
        group_paths["report_embedded_html"],
        title=(
            "WORK_IN_PROGRESS Dot ERF Sanity Check | "
            f"{SENSOR_DISPLAY[args.sensor_type]['label']} | "
            f"{sensor_group_output_tag(sensor_group)}"
        ),
        summary=summary,
        trial_counts=trial_counts_df,
        imbalance=imbalance_df,
        cluster_results=cluster_results,
        figure_paths=figure_paths,
        skipped=skipped,
        embed_images=True,
        include_cluster_table=bool(args.include_cluster_table),
    )
    return summary


def validation_summary(
    *,
    args: argparse.Namespace,
    shared_cfg: Mapping[str, Any],
    subject_headers: list[Mapping[str, Any]],
    skipped: list[Mapping[str, Any]],
    channel_names: list[str],
    sensor_group: str,
    sensor_selection: Mapping[str, Any],
    imbalance: pd.DataFrame,
) -> dict[str, Any]:
    """Build read-only validation output for --validate-only mode."""
    paths = group_output_paths(shared_cfg, args.sensor_type, sensor_group)
    return {
        "stage": "dot_erf_validation",
        "status": "validated",
        "analysis_status": "WORK_IN_PROGRESS",
        "work_in_progress_warning": WORK_IN_PROGRESS_TEXT,
        "subjects": [header["subject"] for header in subject_headers],
        "skipped": skipped,
        "epochs": {header["subject"]: str(header["epochs_fif"]) for header in subject_headers},
        "source_manifests": {
            header["subject"]: str(header["source_manifest"]) for header in subject_headers
        },
        "trial_counts": {
            header["subject"]: dict(header["trial_counts"]) for header in subject_headers
        },
        "sensor_type": args.sensor_type,
        "sensor_group": sensor_group_output_tag(sensor_group),
        "sensor_selection": dict(sensor_selection),
        "n_common_good_channels": int(len(channel_names)),
        "common_good_channels": list(channel_names),
        "baseline_interval_sec": [float(args.baseline[0]), float(args.baseline[1])],
        "time_window_request_sec": [float(args.tmin), float(args.tmax)],
        "analysis_sfreq_hz": float(args.analysis_sfreq),
        "topomap_times_sec": [float(value) for value in args.topomap_times],
        "group_overlay_band": args.group_overlay_band,
        "imbalance_diagnostics": imbalance.to_dict("records"),
        "equalized_sensitivity_would_run": bool(
            imbalance["trigger_equalized_sensitivity"].any()
        ),
        "cluster_parameters": cluster_parameter_summary(args),
        "cluster1d_parameters": cluster1d_parameter_summary(args),
        "group_outputs": {key: str(value) for key, value in paths.items()},
        "subject_outputs": {
            header["subject"]: {
                key: str(value)
                for key, value in subject_output_paths(
                    header["cfg"], args.sensor_type, sensor_group
                ).items()
            }
            for header in subject_headers
        },
    }


def main() -> None:
    """Run the end-to-end WORK_IN_PROGRESS ERF sanity-check workflow."""
    args = build_arg_parser().parse_args()
    validate_cli_args(args)
    shared_cfg, shared_warnings = load_shared_config(args)
    resolve_topomap_times(args, shared_cfg)
    raw_root = Path(shared_cfg["resolved_paths"]["data_root"])

    subject_headers, skipped = collect_subject_headers(args)
    sensor_group_definitions = load_stage07_sensor_groups()
    imbalance = imbalance_diagnostics(
        subject_headers,
        threshold=float(args.imbalance_ratio_threshold),
    )
    need_equalized_sensitivity = bool(imbalance["trigger_equalized_sensitivity"].any())

    if args.validate_only:
        summaries: list[dict[str, Any]] = []
        for sensor_group in args.sensor_groups:
            channel_names, sensor_selection = select_sensor_group_channels(
                subject_headers,
                sensor_group=sensor_group,
                sensor_group_definitions=sensor_group_definitions,
            )
            summary = validation_summary(
                args=args,
                shared_cfg=shared_cfg,
                subject_headers=subject_headers,
                skipped=skipped,
                channel_names=channel_names,
                sensor_group=sensor_group,
                sensor_selection=sensor_selection,
                imbalance=imbalance,
            )
            summaries.append(summary)
        if args.print_summary:
            payload: Any = summaries[0] if len(summaries) == 1 else summaries
            print(json.dumps(_jsonify(payload), indent=2, sort_keys=True))
        else:
            print(
                "Dot ERF input validation passed for sensor groups: "
                + ", ".join(args.sensor_groups)
            )
        return

    summaries = []
    for sensor_group in args.sensor_groups:
        group_paths = group_output_paths(shared_cfg, args.sensor_type, sensor_group)
        assert_output_path_outside_raw(group_paths["out_dir"], raw_root)
        assert_output_path_outside_raw(group_paths["report_dir"], raw_root)
        channel_names, sensor_selection = select_sensor_group_channels(
            subject_headers,
            sensor_group=sensor_group,
            sensor_group_definitions=sensor_group_definitions,
        )
        reference_info = make_reference_info(subject_headers[0], channel_names)
        expected_outputs = expected_output_files(
            subject_headers=subject_headers,
            group_paths=group_paths,
            sensor_type=args.sensor_type,
            sensor_group=sensor_group,
        )
        if args.resume and all(path.is_file() for path in expected_outputs):
            print(
                f"Dot ERF outputs already exist for {sensor_group}; "
                "resume mode reuses them."
            )
            for path in expected_outputs:
                print(path)
            summaries.append(
                {
                    "sensor_group": sensor_group,
                    "status": "reused_existing_outputs",
                    "group_summary": str(group_paths["group_summary_json"]),
                    "embedded_report": str(group_paths["report_embedded_html"]),
                }
            )
            continue
        if not args.resume and not args.overwrite:
            existing = [str(path) for path in expected_outputs if path.exists()]
            if existing:
                raise RuntimeError(
                    "Refusing to overwrite existing dot ERF outputs without "
                    f"--overwrite or --resume for sensor group {sensor_group}: {existing}"
                )

        summary = build_group_outputs(
            args=args,
            shared_cfg=shared_cfg,
            subject_headers=subject_headers,
            skipped=skipped,
            channel_names=channel_names,
            sensor_group=sensor_group,
            sensor_selection=sensor_selection,
            reference_info=reference_info,
            need_equalized_sensitivity=need_equalized_sensitivity,
        )
        summaries.append(summary)

    if args.print_summary:
        payload = summaries[0] if len(summaries) == 1 else summaries
        print(json.dumps(_jsonify(payload), indent=2, sort_keys=True))
    print("Dot ERF WORK_IN_PROGRESS sanity-check complete.")
    for summary in summaries:
        outputs = summary.get("outputs", {})
        print(f"Sensor group: {summary.get('sensor_group', summary.get('sensor_selection', {}).get('sensor_group', 'unknown'))}")
        if "group_summary_json" in outputs:
            print(f"Group summary: {outputs['group_summary_json']}")
        elif "group_summary" in summary:
            print(f"Group summary: {summary['group_summary']}")
        if "report_embedded_html" in outputs:
            print(f"Embedded report: {outputs['report_embedded_html']}")
        elif "embedded_report" in summary:
            print(f"Embedded report: {summary['embedded_report']}")
    if shared_warnings:
        print("Configuration warnings:")
        for warning in shared_warnings:
            print(f"- {warning}")


if __name__ == "__main__":
    main()
