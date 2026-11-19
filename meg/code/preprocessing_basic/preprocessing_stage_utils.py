#!/usr/bin/env python3
"""
preprocessing_stage_utils.py

Purpose
Provide shared, importable infrastructure for the configuration-driven DAD MEG
preprocessing stages. This module is not a numbered scientific transformation.
It keeps acquisition-run lineage, guarded output handling, tuning decisions,
and source-manifest updates consistent across 03_A, 03_B, and 03_C.

Usage examples
1) Build lineage records for the configured raw acquisition runs:

   from preprocessing_stage_utils import discover_raw_run_records
   records = discover_raw_run_records(cfg)

2) Write or update an interactive tuning file outside the raw tree:

   update_tuning_file(
       path="/tmp/sub-01_tuning.yaml",
       raw_root=cfg["resolved_paths"]["data_root"],
       subject_id="01",
       section="artifact_annotation",
       values={"muscle_threshold": 8.0},
   )

3) Save only one recording's review, retaining other concurrently saved runs:

   update_tuning_file(
       path="/tmp/sub-04_tuning.yaml", raw_root="/read-only/raw",
       subject_id="04", section="artifact_annotation",
       mapping_key="manual_eye_reviews",
       values={"sub04_run01_220526.fif": receipt},
   )

Data flow
resolved configuration -> raw acquisition-run discovery -> stable run records
-> guarded derivative paths -> stage reports/tuning/manifest lineage.

Safety assumptions
- Raw files are read-only inputs.
- Every write path is validated with
  ``pipeline_config.assert_output_path_outside_raw`` before parent-directory
  creation or file writing.
- Behavioral block labels embedded in a raw filename remain part of the source
  filename only. Acquisition-run identity is parsed exclusively from ``runNN``.
"""

from __future__ import annotations

import fcntl
import html
import json
import os
import re
import stat
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import yaml

from dataset_inventory import discover_meg_fif_files
from pipeline_config import assert_output_path_outside_raw


DERIVATIVE_LINEAGE_SCHEMA_VERSION = "1.0"


# -----------------------------------------------------------------------------
# Acquisition-run discovery and stable naming
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class AcquisitionRunRecord:
    """
    Describe one raw task acquisition run without inferring behavioral mapping.

    Attributes
    ----------
    source_path
        Absolute raw FIF path.
    source_filename
        Original FIF filename, preserved verbatim for provenance.
    source_stem
        Original filename without the final ``.fif`` suffix.
    acquisition_run_label
        Integer parsed only from the source filename's ``runNN`` token.
    raw_relative_path
        Path relative to the configured read-only raw root.
    """

    source_path: str
    source_filename: str
    source_stem: str
    acquisition_run_label: int
    raw_relative_path: str


def utc_now() -> str:
    """Return the current UTC timestamp in stable ISO8601 form."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def parse_acquisition_run_label(filename: str) -> int:
    """
    Parse the acquisition-run number from a raw or derivative filename.

    The parser intentionally ignores any ``blockNN`` tokens. This prevents the
    historical behavioral-block naming convention from becoming a positional
    MEG-to-behavior association.
    """
    match = re.search(r"(?:^|_)run(\d+)(?:_|\.|$)", filename, re.IGNORECASE)
    if match is None:
        raise ValueError(
            "Could not parse an acquisition-run token from filename: "
            f"{filename!r}"
        )
    return int(match.group(1))


def discover_raw_run_records(cfg: Mapping[str, Any]) -> list[AcquisitionRunRecord]:
    """
    Discover configured raw task FIFs and return acquisition-run records.

    Returns
    -------
    list of AcquisitionRunRecord
        Records sorted by acquisition-run label and then source filename.
    """
    resolved = cfg["resolved_paths"]
    paths = discover_meg_fif_files(
        resolved["subject_data_dir"],
        cfg["data_discovery"]["meg_fif_glob"],
    )
    if not paths:
        raise FileNotFoundError(
            "No configured task acquisition FIFs were discovered under "
            f"{resolved['subject_data_dir']}."
        )

    raw_root = Path(resolved["data_root"]).resolve()
    records = [
        AcquisitionRunRecord(
            source_path=str(path),
            source_filename=path.name,
            source_stem=path.stem,
            acquisition_run_label=parse_acquisition_run_label(path.name),
            raw_relative_path=str(path.relative_to(raw_root)),
        )
        for path in paths
    ]
    records.sort(
        key=lambda item: (item.acquisition_run_label, item.source_filename)
    )
    labels = [item.acquisition_run_label for item in records]
    if len(labels) != len(set(labels)):
        raise ValueError(
            "Task FIF acquisition-run labels must be unique within one subject: "
            f"{labels}"
        )
    return records


def derivative_path_for_run(
    record: AcquisitionRunRecord,
    output_dir: str | Path,
    suffix: str,
) -> Path:
    """
    Build a derivative path while preserving the complete raw source stem.

    ``suffix`` should include the leading underscore and file extension, for
    example ``"_raw_sss.fif"``, ``"_ann.fif"``, or ``"_ica.fif"``.
    """
    if not suffix.startswith("_"):
        raise ValueError(f"Derivative suffix must start with '_': {suffix!r}")
    return Path(output_dir).expanduser().resolve() / f"{record.source_stem}{suffix}"


# -----------------------------------------------------------------------------
# Guarded output handling
# -----------------------------------------------------------------------------


def guarded_path(path: str | Path, raw_root: str | Path) -> Path:
    """Resolve and validate one output path without creating it."""
    return assert_output_path_outside_raw(path, raw_root)


def prepare_output_path(
    path: str | Path,
    *,
    raw_root: str | Path,
    overwrite: bool,
) -> Path:
    """
    Validate one output and require explicit overwrite permission if it exists.

    Parent directories are created only after the raw-tree guard succeeds.
    """
    resolved = guarded_path(path, raw_root)
    if resolved.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists and overwrite is disabled: {resolved}"
        )
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return resolved


def ensure_output_directory(path: str | Path, raw_root: str | Path) -> Path:
    """Validate and create one derivative output directory."""
    resolved = guarded_path(path, raw_root)
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def _json_default(value: Any) -> Any:
    """
    Convert common scientific scalar/container types at the JSON boundary.

    MNE and NumPy may return ``np.int64`` values even when a report field is
    conceptually an integer. Normalization is centralized here so stage code
    can retain the source library's types while JSON remains standards
    compliant. Unknown objects still raise rather than being stringified
    silently.
    """
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(
        f"Object of type {type(value).__name__} is not JSON serializable"
    )


def write_json(
    path: str | Path,
    payload: Mapping[str, Any] | list[Any],
    *,
    raw_root: str | Path,
    overwrite: bool = True,
) -> Path:
    """Write guarded JSON with deterministic formatting."""
    resolved = prepare_output_path(
        path,
        raw_root=raw_root,
        overwrite=overwrite,
    )
    resolved.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
            default=_json_default,
        ),
        encoding="utf-8",
    )
    return resolved


def write_text(
    path: str | Path,
    text: str,
    *,
    raw_root: str | Path,
    overwrite: bool = True,
) -> Path:
    """Write guarded UTF-8 text."""
    resolved = prepare_output_path(
        path,
        raw_root=raw_root,
        overwrite=overwrite,
    )
    resolved.write_text(text, encoding="utf-8")
    return resolved


def write_simple_html_report(
    path: str | Path,
    *,
    title: str,
    summary_lines: Iterable[str],
    image_paths: Iterable[str | Path],
    raw_root: str | Path,
) -> Path:
    """
    Write a compact standalone HTML report linking stage QC images.

    Images are referenced by relative path when possible so reports remain
    portable inside the derivative subject directory.
    """
    report_path = guarded_path(path, raw_root)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    image_blocks: list[str] = []
    for image_path_value in image_paths:
        image_path = Path(image_path_value).expanduser().resolve()
        try:
            href = image_path.relative_to(report_path.parent)
        except ValueError:
            href = image_path
        image_blocks.append(
            "<figure>"
            f"<img src=\"{html.escape(str(href))}\" "
            "style=\"max-width:100%;height:auto\"/>"
            f"<figcaption>{html.escape(image_path.name)}</figcaption>"
            "</figure>"
        )
    body = "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\">",
            f"<title>{html.escape(title)}</title>",
            "<style>body{font-family:system-ui,sans-serif;max-width:1200px;"
            "margin:2rem auto;padding:0 1rem}pre{white-space:pre-wrap;"
            "background:#f5f5f5;padding:1rem}figure{margin:2rem 0}</style>",
            "</head><body>",
            f"<h1>{html.escape(title)}</h1>",
            f"<pre>{html.escape(chr(10).join(summary_lines))}</pre>",
            *image_blocks,
            "</body></html>",
        ]
    )
    report_path.write_text(body, encoding="utf-8")
    return report_path


# -----------------------------------------------------------------------------
# Interactive tuning persistence
# -----------------------------------------------------------------------------


def _read_optional_yaml_mapping(path: Path) -> dict[str, Any]:
    """Read an existing YAML mapping or return an empty mapping."""
    if not path.exists():
        return {}
    parsed = yaml.safe_load(path.read_text(encoding="utf-8"))
    if parsed is None:
        return {}
    if not isinstance(parsed, dict):
        raise ValueError(f"Tuning file must contain a YAML mapping: {path}")
    return parsed


def _plain_tuning_values(value: Any) -> Any:
    """Normalize scientific values while preserving native YAML keys/dates."""
    if isinstance(value, Mapping):
        return {_plain_tuning_values(key): _plain_tuning_values(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain_tuning_values(item) for item in value]
    if isinstance(value, (np.generic, np.ndarray, Path)):
        return _plain_tuning_values(_json_default(value))
    return value


def update_tuning_file(
    *,
    path: str | Path,
    raw_root: str | Path,
    subject_id: str,
    section: str,
    values: Mapping[str, Any],
    mapping_key: str | None = None,
) -> Path:
    """
    Merge one accepted interactive decision into a reusable tuning YAML file.

    The generated file can be supplied to ``run_pipeline.py --tuning-config``
    for the noninteractive apply pass. Normalize NumPy scalars/arrays and Paths
    with the existing report serializer before YAML encoding; scientific library
    return types must not turn a completed interactive review into a save error.

    When mapping_key is supplied, values contains entries of that nested mapping
    (for example one exact-source receipt), not a whole map loaded at GUI startup.
    Each supplied entry replaces its complete old value; other entries survive.
    A short advisory lock covers reading, merging and atomic replacement on
    macOS/Linux. The persistent hidden .lock file is only a write lock. This
    protects different runs saved by overlapping interactive processes; it does
    not merge competing annotation edits of the same recording.
    """
    resolved = guarded_path(path, raw_root)
    resolved.parent.mkdir(parents=True, exist_ok=True)
    lock_path = guarded_path(resolved.with_name(f".{resolved.name}.lock"), raw_root)
    # Read the latest file inside the lock, after the potentially long GUI
    # session. Never delete the lock file: another process may be waiting on it.
    with lock_path.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        payload = _read_optional_yaml_mapping(resolved)
        payload.setdefault("subject", {})["id"] = str(subject_id)
        current = payload.get(section)
        if current is None:
            current = {}
        if not isinstance(current, dict):
            raise ValueError(f"Existing tuning section {section!r} must be a mapping in {resolved}.")
        if mapping_key is None:
            current.update(dict(values))
        else:
            entries = current.get(mapping_key)
            if entries is None:
                entries = {}
            if not isinstance(entries, dict):
                raise ValueError(f"Tuning {section}.{mapping_key} must be a mapping in {resolved}.")
            entries.update(dict(values))
            current[mapping_key] = entries
        payload[section] = current
        # Serialize before touching the old file. A same-directory replacement
        # exposes a complete YAML to readers and retains existing permissions.
        encoded = yaml.safe_dump(_plain_tuning_values(payload), sort_keys=False)
        mode = stat.S_IMODE(resolved.stat().st_mode) if resolved.exists() else 0o644
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=resolved.parent,
                    prefix=f".{resolved.name}.", suffix=".tmp", delete=False) as handle:
                temporary = Path(handle.name)
                handle.write(encoded)
                handle.flush()
                os.fsync(handle.fileno())
            temporary.chmod(mode)
            os.replace(temporary, resolved)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    return resolved


# -----------------------------------------------------------------------------
# Source-manifest derivative lineage
# -----------------------------------------------------------------------------


def lineage_record(
    *,
    source: AcquisitionRunRecord,
    input_path: str | Path,
    output_path: str | Path | None,
    status: str,
    metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build one stage lineage record for a source acquisition run."""
    record = {
        **asdict(source),
        "input_path": str(Path(input_path).expanduser().resolve()),
        "output_path": (
            str(Path(output_path).expanduser().resolve())
            if output_path is not None
            else None
        ),
        "status": status,
    }
    if metadata:
        record["metadata"] = dict(metadata)
    return record


def update_source_manifest_stage(
    *,
    manifest_path: str | Path | None,
    raw_root: str | Path,
    stage_name: str,
    processing_pass: str,
    records: list[Mapping[str, Any]],
    status: str,
    qc: Mapping[str, Any] | None = None,
) -> str:
    """
    Add or replace one stage's derivative lineage in the durable manifest.

    Missing manifests are allowed for standalone stage execution and return an
    explicit status rather than silently dropping the provenance update.
    """
    if manifest_path is None:
        return "not_updated_missing_path"
    resolved = guarded_path(manifest_path, raw_root)
    if not resolved.exists():
        return "not_updated_missing_file"
    payload = json.loads(resolved.read_text(encoding="utf-8"))
    derivative_lineage = payload.setdefault(
        "derivative_lineage",
        {
            "schema_version": DERIVATIVE_LINEAGE_SCHEMA_VERSION,
            "stages": {},
        },
    )
    derivative_lineage.setdefault("stages", {})[stage_name] = {
        "updated_utc": utc_now(),
        "processing_pass": processing_pass,
        "status": status,
        "records": records,
        "qc": dict(qc or {}),
    }
    resolved.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return "updated"


def stage_output_dirs(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """
    Return canonical stage output directories below the derivative root.

    Directory creation is intentionally left to callers so every requested
    path can be guarded before side effects.
    """
    preprocessing = Path(
        cfg["resolved_paths"]["derivatives_preprocessing_dir"]
    ).resolve()
    reports = Path(cfg["resolved_paths"]["derivatives_reports_dir"]).resolve()
    return {
        "sss": preprocessing / "sss",
        "figures": preprocessing / "figures",
        "ica_models": preprocessing / "ica_models",
        "reports": reports,
    }


__all__ = [
    "AcquisitionRunRecord",
    "derivative_path_for_run",
    "discover_raw_run_records",
    "ensure_output_directory",
    "guarded_path",
    "lineage_record",
    "parse_acquisition_run_label",
    "prepare_output_path",
    "stage_output_dirs",
    "update_source_manifest_stage",
    "update_tuning_file",
    "utc_now",
    "write_json",
    "write_simple_html_report",
    "write_text",
]
