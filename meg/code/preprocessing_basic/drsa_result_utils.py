#!/usr/bin/env python3
"""
drsa_result_utils.py

Purpose
Shared Python helpers for lightweight post-processing of Stage-07 dRSA result
MAT files. The helpers read MATLAB v7 and v7.3 result files, normalize dRSA
matrix arrays, infer metadata from the Stage-07 folder layout, smooth matrices
for display-only summaries while preserving undefined matrix support, convert
occlusion frames to seconds, and create self-contained embedded HTML copies of
static reports.

Usage examples
1) Read all dRSA matrices from one Stage-07 result file:

   from pathlib import Path
   from drsa_result_utils import read_drsa_result_matrices

   matrices, time_axis = read_drsa_result_matrices(
       Path("/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/"
            "sub-01/ses-01/drsa/matlab_real_neural/07_dRSA/"
            "always_visible/corr/results/dRSA_always_visible_corr_results.mat")
   )

2) Convert the configured 1-based occlusion frames to seconds:

   timing = occlusion_timing_from_config(cfg)
   assert timing["first_complete_occlusion_sec"] == 129 / 120

Data flow
Stage-07 result MAT files -> numeric matrix variables and time/model metadata
-> Stage-08 line-summary calculations and reports.

Assumptions
- Result matrices are treated exactly as Stage 07 plots them: x is neural time
  and y is model time. Therefore a 2-D matrix is indexed
  ``matrix[model_time_index, neural_time_index]``.
- MATLAB v7.3 files require h5py. MATLAB v7 files require scipy.
- These helpers never write into the raw data tree; callers pass report/output
  paths through the shared guard before writing.
"""

from __future__ import annotations

import base64
import html
import mimetypes
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


CONDITIONS = ("always_visible", "occluded_nondeviant", "occluded_deviant")
DRSA_MATRIX_VARIABLES = (
    "dRSA_neural_position",
    "dRSA_neural_correlation",
    "dRSA_neural_self",
    "dRSA_model_pos_pos",
    "dRSA_model_pos_dir",
)
KNOWN_SENSOR_GROUPS = (
    "wholebrain",
    "occipital",
    "parietal",
    "frontal",
    "temporal",
)


def utc_now() -> str:
    """Return an ISO8601 UTC timestamp with an explicit Z suffix."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sanitize_label(value: Any) -> str:
    """
    Convert an arbitrary label into a conservative filename/anchor token.

    The Stage-08 report paths are generated from analysis metadata that may
    contain spaces, underscores, MATLAB variable names, or model labels. This
    helper keeps output filenames stable and shell-friendly.
    """
    text = str(value).strip().lower()
    safe = "".join(char if char.isalnum() else "-" for char in text).strip("-")
    while "--" in safe:
        safe = safe.replace("--", "-")
    return safe or "unknown"


def html_id(value: Any) -> str:
    """Create a stable HTML id from a display label."""
    return sanitize_label(value)


def image_href(report_path: Path, image_path: str | Path) -> str:
    """Return a report-local image href when possible, else an absolute path."""
    resolved_image = Path(image_path).expanduser().resolve()
    try:
        return str(resolved_image.relative_to(report_path.parent))
    except ValueError:
        return str(resolved_image)


def write_embedded_html_copy(
    report_path: str | Path,
    *,
    raw_root: str | Path,
) -> Path:
    """
    Write a sibling HTML file with local image references embedded as data URIs.

    Data flow:
    static HTML report -> resolve ``img src`` paths -> base64-encoded image
    bytes -> sibling ``*_embedded.html`` report.
    """
    from preprocessing_stage_utils import guarded_path

    report = Path(report_path).expanduser().resolve()
    output = report.with_name(f"{report.stem}_embedded{report.suffix}")
    guarded_path(output, raw_root)
    html_text = report.read_text(encoding="utf-8")

    def replace_src(match: re.Match[str]) -> str:
        quote = match.group(1)
        src = match.group(2)
        if src.startswith("data:"):
            return match.group(0)
        image_path = Path(src).expanduser()
        if not image_path.is_absolute():
            image_path = (report.parent / image_path).resolve()
        if not image_path.is_file():
            return match.group(0)
        mime_type = mimetypes.guess_type(image_path.name)[0] or "image/png"
        encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
        return f"src={quote}data:{mime_type};base64,{encoded}{quote}"

    embedded = re.sub(r"src=(['\"])([^'\"]+)\1", replace_src, html_text)
    output.write_text(embedded, encoding="utf-8")
    return output


def frame_to_seconds(frame_index: int | float, frame_frequency_hz: int | float) -> float:
    """
    Convert one 1-based stimulus frame index to seconds.

    The DAD trajectory/stimulus code uses 1-based frame indices, while the MEG
    and dRSA time axes start at 0 s. The conversion is therefore
    ``(frame - 1) / frame_frequency_hz``.
    """
    return (float(frame_index) - 1.0) / float(frame_frequency_hz)


def occlusion_timing_from_config(cfg: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve occlusion timing and derived interval boundaries from config."""
    timing_cfg = cfg.get("matlab_drsa", {}).get("occlusion_timing", {})
    fps = float(timing_cfg.get("frame_frequency_hz", 120.0))
    first_frame = int(timing_cfg.get("first_complete_occlusion_frame", 130))
    reappearance_frame = int(timing_cfg.get("reappearance_frame", 191))
    return {
        "frame_frequency_hz": fps,
        "first_complete_occlusion_frame": first_frame,
        "reappearance_frame": reappearance_frame,
        "first_complete_occlusion_sec": frame_to_seconds(first_frame, fps),
        "reappearance_sec": frame_to_seconds(reappearance_frame, fps),
    }


def occlusion_intervals(
    time_axis_sec: Any,
    occlusion_timing: Mapping[str, Any],
) -> dict[str, Any]:
    """
    Build model-time interval masks for pre, during, and post occlusion.

    Intervals follow the requested convention:
    - pre_occlusion: [0, first_complete_occlusion_sec)
    - during_occlusion: [first_complete_occlusion_sec, reappearance_sec)
    - post_occlusion: [reappearance_sec, end]
    """
    import numpy as np

    axis = np.asarray(time_axis_sec, dtype=float).squeeze()
    first = float(occlusion_timing["first_complete_occlusion_sec"])
    reappearance = float(occlusion_timing["reappearance_sec"])
    return {
        "pre_occlusion": (axis >= 0.0) & (axis < first),
        "during_occlusion": (axis >= first) & (axis < reappearance),
        "post_occlusion": axis >= reappearance,
    }


def _is_matlab_v73_file(path: Path) -> bool:
    """Detect MATLAB v7.3/HDF5 MAT files from their cheap header bytes."""
    try:
        with path.open("rb") as file_obj:
            header = file_obj.read(128)
    except OSError:
        return False
    return b"MATLAB 7.3" in header or header.startswith(b"\x89HDF")


def _load_hdf5_mat_arrays(path: Path) -> dict[str, Any]:
    """Load selected numeric arrays from a MATLAB v7.3/HDF5 result file."""
    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "Cannot read MATLAB v7.3 result MAT files because h5py is not "
            f"installed. Install h5py or use a Python environment that provides "
            f"it. File: {path}"
        ) from exc

    payload: dict[str, Any] = {}
    with h5py.File(path, "r") as handle:
        for name in (*DRSA_MATRIX_VARIABLES, "time_axis_sec"):
            if name not in handle:
                continue
            array = handle[name][()]
            # MATLAB v7.3 stores arrays in HDF5 dimension order; reversing axes
            # restores the usual MATLAB/scipy [time x time x model] layout.
            if getattr(array, "ndim", 0) > 1:
                array = array.transpose(tuple(reversed(range(array.ndim))))
            payload[name] = array
        for name in ("model_names", "model_dist"):
            if name in handle:
                payload[name] = _read_hdf5_text_values(handle, name)
    return payload


def _load_mat_arrays(path: Path) -> dict[str, Any]:
    """Load selected arrays from one MATLAB v7 or v7.3 result file."""
    if _is_matlab_v73_file(path):
        return _load_hdf5_mat_arrays(path)

    try:
        from scipy.io import loadmat
    except ImportError as exc:
        raise RuntimeError(
            "Reading MATLAB v7 result files requires scipy in the active Python "
            f"environment. File: {path}"
        ) from exc

    try:
        return loadmat(
            path,
            variable_names=[
                *DRSA_MATRIX_VARIABLES,
                "time_axis_sec",
                "model_names",
                "model_dist",
            ],
            squeeze_me=True,
            struct_as_record=False,
        )
    except (NotImplementedError, ValueError) as exc:
        message = str(exc).lower()
        if "7.3" in message or "hdf" in message or "unknown mat file type" in message:
            return _load_hdf5_mat_arrays(path)
        raise


def as_numeric_matrix(value: Any, *, label: str, path: Path):
    """
    Normalize one loaded MAT variable to a 2-D or 3-D float array.

    Stage-07 single-model outputs can squeeze to 2-D. Multi-model outputs are
    represented as [model_time x neural_time x model_index]. Some HDF5 readers
    expose dimensions in reversed order, so this function repairs the common
    [model_index x time x time] layout when necessary.
    """
    import numpy as np

    array = np.asarray(value, dtype=float).squeeze()
    if array.ndim not in {2, 3}:
        raise RuntimeError(
            f"{label} in {path} must be a 2-D or 3-D numeric matrix, "
            f"found shape {array.shape}."
        )
    if array.ndim == 3 and array.shape[0] != array.shape[1]:
        if array.shape[1] == array.shape[2]:
            array = np.moveaxis(array, 0, -1)
        else:
            raise RuntimeError(
                f"{label} in {path} has unsupported 3-D shape {array.shape}."
            )
    return array


def _read_hdf5_text_values(handle: Any, name: str) -> list[str]:
    """
    Read a MATLAB v7.3 char/cellstr dataset into a list of Python strings.

    MATLAB cell arrays are stored as HDF5 object references, while char arrays
    are stored as integer character-code matrices. This decoder handles the
    simple string metadata used by Stage-07 result files.
    """
    import numpy as np

    dataset = handle[name]
    data = dataset[()]

    def decode_char_array(value: Any) -> str:
        array = np.asarray(value)
        if array.dtype.kind in {"u", "i"}:
            chars = [chr(int(item)) for item in array.ravel(order="F") if int(item) != 0]
            return "".join(chars)
        if array.dtype.kind in {"S", "U"}:
            return "".join(str(item.decode("utf-8") if isinstance(item, bytes) else item) for item in array.ravel())
        return str(value)

    if getattr(data, "dtype", None) is not None and data.dtype.kind == "O":
        values: list[str] = []
        for ref in data.ravel(order="F"):
            if not ref:
                continue
            values.append(decode_char_array(handle[ref][()]))
        return values
    if getattr(data, "ndim", 0) <= 1:
        decoded = decode_char_array(data)
        return [decoded] if decoded else []
    # For a char matrix, each MATLAB column can represent one string. If that
    # guess is wrong, the fallback still returns a useful joined label.
    try:
        return [
            decode_char_array(data[:, index])
            for index in range(data.shape[1])
            if decode_char_array(data[:, index])
        ]
    except Exception:  # pragma: no cover - defensive for rare MAT encodings.
        decoded = decode_char_array(data)
        return [decoded] if decoded else []


def _normalize_text_values(value: Any) -> list[str]:
    """Normalize scipy/h5py-loaded MATLAB string metadata to a Python list."""
    import numpy as np

    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    array = np.asarray(value, dtype=object).squeeze()
    if array.ndim == 0:
        return [str(array.item())]
    return [str(item) for item in array.ravel(order="C")]


def read_drsa_result_file(path: Path) -> tuple[dict[str, Any], Any | None, dict[str, list[str]]]:
    """
    Read all Stage-08-consumable dRSA matrix variables from one result MAT file.

    Returns
    -------
    matrices
        Mapping from MATLAB variable name to a 2-D or 3-D float array.
    time_axis
        Optional 1-D time axis in seconds if the MAT file provides it.
    """
    import numpy as np

    payload = _load_mat_arrays(path)
    matrices: dict[str, Any] = {}
    for variable in DRSA_MATRIX_VARIABLES:
        if variable in payload:
            matrices[variable] = as_numeric_matrix(
                payload[variable],
                label=variable,
                path=path,
            )
    time_axis = None
    if "time_axis_sec" in payload:
        candidate = np.asarray(payload["time_axis_sec"], dtype=float).squeeze()
        if candidate.ndim == 1 and candidate.size > 1:
            time_axis = candidate
    metadata = {
        "model_names": _normalize_text_values(payload.get("model_names")),
        "model_dist": _normalize_text_values(payload.get("model_dist")),
    }
    return matrices, time_axis, metadata


def read_drsa_result_matrices(path: Path) -> tuple[dict[str, Any], Any | None]:
    """
    Read all Stage-08-consumable dRSA matrix variables from one result MAT file.

    This compatibility wrapper returns only numeric matrices and time metadata.
    Use ``read_drsa_result_file`` when model-name metadata is needed.
    """
    matrices, time_axis, _metadata = read_drsa_result_file(path)
    return matrices, time_axis


def smooth_matrix_for_plot(matrix_2d: Any, window_samples: int):
    """
    Apply NaN-aware Gaussian smoothing to one 2-D display matrix.

    This is display-only smoothing for exploratory summaries. It does not
    require or write smoothed MAT files. Cells that are non-finite before
    smoothing remain NaN afterward, so masks such as undefined prediction-error
    rows cannot leak into interval summaries or display-only group averages.
    """
    import numpy as np
    from scipy.ndimage import gaussian_filter

    matrix = np.asarray(matrix_2d, dtype=float)
    if int(window_samples) <= 1:
        return matrix.copy()
    sigma = max(float(window_samples) / 6.0, 0.5)
    finite = np.isfinite(matrix)
    values = np.where(finite, matrix, 0.0)
    weights = finite.astype(float)
    smooth_values = gaussian_filter(values, sigma=sigma, mode="nearest")
    smooth_weights = gaussian_filter(weights, sigma=sigma, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        smoothed = smooth_values / smooth_weights
    smoothed[smooth_weights <= 1e-12] = np.nan
    smoothed[~finite] = np.nan
    return smoothed


def smooth_array_for_plot(matrix: Any, window_samples: int):
    """Smooth each 2-D slice in a 2-D or 3-D dRSA result array."""
    import numpy as np

    array = np.asarray(matrix, dtype=float)
    if array.ndim == 2:
        return smooth_matrix_for_plot(array, window_samples)
    slices = [
        smooth_matrix_for_plot(array[:, :, i_model], window_samples)
        for i_model in range(array.shape[2])
    ]
    return np.stack(slices, axis=2)


def iter_matrix_slices(
    matrix: Any,
    slice_labels: Sequence[str] | None = None,
) -> Iterable[tuple[int, str, Any]]:
    """Yield ``(slice_index, slice_label, matrix_2d)`` from a 2-D/3-D array."""
    import numpy as np

    array = np.asarray(matrix, dtype=float)
    if array.ndim == 2:
        yield 0, "matrix", array
        return
    for i_model in range(array.shape[2]):
        label = (
            str(slice_labels[i_model])
            if slice_labels is not None and i_model < len(slice_labels)
            else f"model{i_model + 1:02d}"
        )
        yield i_model, label, array[:, :, i_model]


def time_axis_for_matrix(
    matrix_2d: Any,
    *,
    time_axis: Any | None,
    sampling_rate_hz: float,
):
    """Return a model/neural time axis matching one square dRSA matrix."""
    import numpy as np

    n_time = int(matrix_2d.shape[0])
    if time_axis is not None:
        axis = np.asarray(time_axis, dtype=float).squeeze()
        if axis.ndim == 1 and axis.size == n_time:
            return axis
    return np.arange(n_time, dtype=float) / float(sampling_rate_hz)


def sampling_step_seconds(time_axis_sec: Any, fallback_hz: float) -> float:
    """Estimate sample spacing in seconds from a time axis or fallback rate."""
    import numpy as np

    axis = np.asarray(time_axis_sec, dtype=float).squeeze()
    if axis.ndim == 1 and axis.size > 1:
        diffs = np.diff(axis)
        diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
        if diffs.size:
            return float(np.median(diffs))
    return 1.0 / float(fallback_hz)


def condition_from_path(path: Path) -> str:
    """Infer the dRSA condition from a Stage-07 result path."""
    parts = [part.lower() for part in path.parts]
    name = path.name.lower()
    for condition in CONDITIONS:
        if condition.lower() in parts or condition.lower() in name:
            return condition
    return "unknown"


def analysis_family_from_path(path: Path) -> str:
    """Infer a readable analysis-family label from a Stage-07 result path."""
    text = str(path).lower()
    if "a_drsa_alt" in text:
        return "A_dRSA_alt"
    if "07_drsa" in text:
        return "A_dRSA"
    if "model_vs_model" in text:
        return "simulations_model_vs_model"
    if "neural_vs_neural" in text:
        return "simulations_neural_vs_neural"
    return "dRSA"


def drsa_type_from_path(path: Path) -> str:
    """Infer corr/PCR from path components when available."""
    for part in path.parts:
        lowered = part.lower()
        if lowered in {"corr", "pcr"}:
            return lowered
    name = path.name.lower()
    if "_pcr" in name:
        return "pcr"
    return "corr" if "_corr" in name else "unknown"


def model_source_from_path(path: Path) -> str:
    """Infer observed/predicted/prediction-error source metadata from a path."""
    parts = [part.lower() for part in path.parts]
    if "prediction_error" in parts:
        return "prediction_error"
    if "predicted" in parts:
        return "predicted"
    if "observed" in parts:
        return "observed"
    if "neural_vs_neural" in parts:
        return "neural"
    return "observed"


def _slice_number(slice_label: Any) -> int | None:
    """Return the 1-based model slice number from model01/_model01 labels."""
    match = re.search(r"model(\d+)", str(slice_label).lower())
    return int(match.group(1)) if match else None


def model_source_for_matrix_slice(
    condition: str,
    model_source: str,
    variable: str,
    slice_label: Any = "",
    analysis_family: str = "",
) -> str:
    """
    Resolve the effective model source for one matrix slice.

    Main A_dRSA outputs store multiple model slices in one MAT file. For those
    files, the slice index is the only place where observed, predicted, and
    prediction-error models differ. Single-model A_dRSA_alt outputs instead
    encode observed/predicted in the folder path.
    """
    family = str(analysis_family).lower()
    variable_name = str(variable)
    lowered_label = str(slice_label).lower()
    if variable_name == "dRSA_neural_correlation":
        if "prediction error" in lowered_label or "pe " in lowered_label:
            return "prediction_error"
        if "expected" in lowered_label or "predicted" in lowered_label:
            return "predicted"
        if "actual" in lowered_label or "observed" in lowered_label:
            return "observed"
    slice_num = _slice_number(slice_label)
    if "a_drsa" in family and variable_name == "dRSA_neural_correlation":
        if str(condition) == "occluded_deviant":
            if slice_num in {3, 4}:
                return "predicted"
            if slice_num in {5, 6}:
                return "prediction_error"
        return "observed"
    return str(model_source or "observed").lower()


def model_feature_label(
    variable: str,
    slice_label: Any = "",
    analysis_family: str = "",
) -> str:
    """
    Resolve the represented feature for one dRSA matrix.

    The returned value is intentionally short because it is used in axis labels,
    report captions, and file navigation contexts.
    """
    variable_name = str(variable)
    lowered_label = str(slice_label).lower()
    if variable_name == "dRSA_neural_correlation":
        if "direction" in lowered_label:
            return "direction"
        if "position" in lowered_label:
            return "position"
    slice_num = _slice_number(slice_label)
    if variable_name == "dRSA_neural_correlation":
        if slice_num in {2, 4, 6}:
            return "direction"
        return "position"
    if variable_name in {"dRSA_neural_position", "dRSA_model_pos_pos"}:
        return "position"
    if variable_name == "dRSA_model_pos_dir":
        return "direction"
    if variable_name == "dRSA_neural_self":
        return "neural signal"
    return sanitize_label(variable_name).replace("-", " ")


def stimulus_model_label(
    condition: str,
    model_source: str,
    variable: str,
    slice_label: Any = "",
    analysis_family: str = "",
) -> str:
    """
    Return the human-readable stimulus/comparison model label for a plot.

    Labels are condition-aware so every plot states which stimulus trajectory is
    represented, not only whether a path was technically ``observed`` or
    ``predicted``.
    """
    source = model_source_for_matrix_slice(
        condition,
        model_source,
        variable,
        slice_label=slice_label,
        analysis_family=analysis_family,
    )
    feature = model_feature_label(
        variable,
        slice_label=slice_label,
        analysis_family=analysis_family,
    )
    condition_name = str(condition)
    if source == "predicted":
        prefix = "expected/predicted"
    elif source == "prediction_error":
        prefix = "prediction error"
    elif source == "neural":
        return "neural signal"
    elif source == "model":
        prefix = "stimulus model"
    elif condition_name == "always_visible":
        prefix = "actual visible"
    elif condition_name == "occluded_nondeviant":
        prefix = "actual nondeviant"
    elif condition_name == "occluded_deviant":
        prefix = "actual deviant"
    else:
        prefix = "actual/observed"
    return f"{prefix} {feature}".strip()


def model_time_axis_label(
    condition: str,
    model_source: str,
    variable: str,
    slice_label: Any = "",
    analysis_family: str = "",
) -> str:
    """Return the y/x-axis label for model-time-like plot dimensions."""
    label = stimulus_model_label(
        condition,
        model_source,
        variable,
        slice_label=slice_label,
        analysis_family=analysis_family,
    )
    prefix = "comparison time" if "neural signal" in label else "model time"
    return f"{prefix} (s) - {label}"


def sensor_group_from_path(path: Path) -> str:
    """
    Infer the sensor group represented by a Stage-07 result path.

    Main A_dRSA real-neural outputs are whole-brain by construction. Model-only
    simulations use the synthetic ``model`` group. Sensor-group analyses encode
    the group in the result filename suffix.
    """
    text = str(path).lower()
    if "model_vs_model" in text:
        return "model"
    if "07_drsa" in text:
        return "wholebrain"
    for group in KNOWN_SENSOR_GROUPS:
        if f"_{group}.mat" in path.name.lower() or f"_{group}_" in path.name.lower():
            return group
    return "unspecified"


def normalize_subject_label(value: Any) -> str:
    """Normalize 01/sub-01 style input to a BIDS subject label."""
    text = str(value).strip()
    if not text:
        raise ValueError("Subject labels must not be empty.")
    if text.startswith("sub-"):
        return text
    return f"sub-{int(text):02d}" if text.isdigit() else f"sub-{text}"


def normalize_sensor_group(value: Any) -> str:
    """Normalize sensor group selectors to the Stage-07 lowercase tags."""
    text = str(value).strip().lower().replace("_", "")
    if text == "wholebrain":
        return "wholebrain"
    return text


def html_escape_join(values: Iterable[Any], separator: str = ", ") -> str:
    """Join values after HTML-escaping each item."""
    return html.escape(separator.join(str(value) for value in values))
