#!/usr/bin/env python3
"""
06_export_neural_for_drsa.py

Purpose
Prepare epoched DAD MEG data for dRSA by exporting condition-specific neural
tensors averaged by `trial_sequence_identity`.

High-level flow
1) load the step-04 epoched FIF for one subject/session,
2) read each surviving epoch's `trial_sequence_identity` from epochs metadata,
3) crop epochs to the post-onset analysis window [0, 2.667] seconds,
4) downsample to 120 Hz without applying an explicit low-pass filter,
5) build/load a canonical MEG channel-label list shared across exports,
6) average epochs within each condition by identity 1..20,
7) write averages into the full canonical channel axis and leave bad/missing
   channels as NaN,
8) save one MAT file per condition containing both the neural tensor and
   ordered channel labels, and
9) write figures plus an MNE HTML report documenting counts, mapping checks,
   downsampling traces, and final tensor shapes.

Usage example
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/code/preprocessing_basic/06_export_neural_for_drsa.py

Sandbox/testing example
DAD_MEG_PROJECT_ROOT=/tmp/dad_drsa_sandbox python /tmp/dad_drsa_sandbox/code/preprocessing_basic/06_export_neural_for_drsa.py

Inputs
- Epoched FIF from step 04:
  derivatives/analysis/sub-01/ses-01/meg/sub-01_ses-01_task-Attention_run-all_epo.fif
  This FIF must contain `epochs.metadata['trial_sequence_identity']`.

Outputs
- MAT tensors under:
  derivatives/analysis/sub-01/ses-01/drsa/
  Each MAT file contains:
  - `neural`: shape 20 x n_canonical_meg_channels x n_timepoints
  - `channel_labels`: MEG channel labels in the exact order used by
    `neural` channel dimension.
- Static figures under:
  derivatives/preprocessing/sub-01/ses-01/figures/
- HTML report:
  derivatives/reports/sub-01/ses-01/06_export_neural_for_drsa_report.html

Key assumptions
- The exported signal uses MEG channels only (magnetometers and gradiometers)
  and excludes non-MEG channels such as STI/EOG.
- Channels marked bad in `epochs.info['bads']` (or missing from the current
  FIF but present in the canonical list) are kept in `channel_labels` and
  exported as NaNs in `neural`.
- Canonical channel labels are loaded from:
  `code/preprocessing_basic/drsa_canonical_meg_channels.txt`.
  If this file is absent, it is initialized from the first run's MEG channel
  order and then reused for all later exports.
- Identity mapping is owned by step 04 and is read from the saved epochs
  metadata. Missing metadata is treated as a hard error so stale epoch files do
  not silently fall back to older reconstruction logic.
- Any missing, duplicate, or ambiguous epoch-to-identity mapping is treated as
  a hard error before MAT export.
"""

# %%
# -------------------------------------------------------------------------
# Imports
# -------------------------------------------------------------------------
import html
import os
import os.path as op
from collections import Counter

import matplotlib.pyplot as plt
import mne
import numpy as np
from scipy.io import savemat

# %%
# -------------------------------------------------------------------------
# Subject/session defaults and analysis constants
# -------------------------------------------------------------------------
subj_str = "01"
ses_str = "01"

# Environment override is used only to keep sandbox tests from writing into the
# real derivatives tree. Without this variable, the local project path is used.
project_root = os.environ.get(
    "DAD_MEG_PROJECT_ROOT",
    "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG",
)

# Condition mapping inherited from step 04.
event_id = {
    "always_visible": 31,
    "occluded_nondeviant": 41,
    "occluded_deviant": 51,
}
code_to_condition = {code: name for name, code in event_id.items()}

# dRSA export settings.
identity_values = list(range(1, 21))
export_tmin = 0.0
export_tmax = 2.667
export_sfreq = 120.0

# %%
# -------------------------------------------------------------------------
# Paths and output layout
# -------------------------------------------------------------------------
preproc_dir = op.join(
    project_root,
    "derivatives",
    "preprocessing",
    f"sub-{subj_str}",
    f"ses-{ses_str}",
)

analysis_meg_dir = op.join(
    project_root,
    "derivatives",
    "analysis",
    f"sub-{subj_str}",
    f"ses-{ses_str}",
    "meg",
)
epochs_fif = op.join(
    analysis_meg_dir,
    f"sub-{subj_str}_ses-{ses_str}_task-Attention_run-all_epo.fif",
)
drsa_output_dir = op.join(
    project_root,
    "derivatives",
    "analysis",
    f"sub-{subj_str}",
    f"ses-{ses_str}",
    "drsa",
)
figure_output_dir = op.join(preproc_dir, "figures")
out_report = op.join(
    project_root,
    "derivatives",
    "reports",
    f"sub-{subj_str}",
    f"ses-{ses_str}",
    "06_export_neural_for_drsa_report.html",
)
canonical_channels_txt = op.join(
    project_root,
    "code",
    "preprocessing_basic",
    "drsa_canonical_meg_channels.txt",
)

os.makedirs(drsa_output_dir, exist_ok=True)
os.makedirs(figure_output_dir, exist_ok=True)
os.makedirs(op.dirname(out_report), exist_ok=True)

# %%
# -------------------------------------------------------------------------
# Shared report/figure helpers
# -------------------------------------------------------------------------
report = mne.Report(title=f"sub-{subj_str} | 06 Export neural for dRSA")


def savefig_logged(fig, filename, **kwargs):
    """
    Save a matplotlib figure into the shared preprocessing figure folder.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save.
    filename : str
        Figure filename relative to `figure_output_dir`.
    kwargs : dict
        Additional keyword arguments forwarded to `Figure.savefig`.

    Returns
    -------
    out_path : str
        Absolute output figure path.
    """
    out_path = op.join(figure_output_dir, filename)
    fig.savefig(out_path, **kwargs)
    print(f"Saved figure to: {out_path}")
    return out_path


def log_to_console_and_report(text, *, title, section, tags=("log", "qc")):
    """
    Print text to console and store the same content in the MNE report.

    Parameters
    ----------
    text : str
        Text block to log.
    title : str
        Report item title.
    section : str
        Report section label.
    tags : tuple of str
        MNE report tags.
    """
    print(text)
    report.add_html(
        f"<pre>{html.escape(text)}</pre>",
        title=title,
        section=section,
        tags=tags,
    )


def save_report_before_error(message):
    """
    Save the report before raising a fatal mapping/export error.

    Parameters
    ----------
    message : str
        Fatal error message.
    """
    log_to_console_and_report(
        message,
        title="Fatal validation error",
        section="Identity mapping",
        tags=("error", "identity"),
    )
    report.save(out_report, overwrite=True, open_browser=False)
    print(f"Saved HTML report to: {out_report}")


def load_or_initialize_canonical_meg_channels(*, canonical_path, current_meg_channels):
    """
    Load canonical MEG channel labels, or initialize them from current data.

    Parameters
    ----------
    canonical_path : str
        Absolute path to the canonical channel text file.
    current_meg_channels : list of str
        MEG channel labels observed in the current epochs file.

    Returns
    -------
    canonical_channels : list of str
        Canonical MEG channel labels in fixed export order.
    file_created : bool
        True if the canonical file was created in this run.
    """
    if op.exists(canonical_path):
        with open(canonical_path, "r", encoding="utf-8") as handle:
            canonical_channels = [
                line.strip()
                for line in handle
                if line.strip() and not line.lstrip().startswith("#")
            ]
        if not canonical_channels:
            raise RuntimeError(
                "Canonical MEG channel file exists but has no channel labels.\n"
                f"canonical_channels_txt={canonical_path}"
            )
        duplicate_channels = [
            ch
            for ch, count in Counter(canonical_channels).items()
            if int(count) > 1
        ]
        if duplicate_channels:
            raise RuntimeError(
                "Canonical MEG channel file contains duplicates.\n"
                f"canonical_channels_txt={canonical_path}\n"
                f"duplicate_channels={sorted(duplicate_channels)}"
            )
        return canonical_channels, False

    os.makedirs(op.dirname(canonical_path), exist_ok=True)
    with open(canonical_path, "w", encoding="utf-8") as handle:
        handle.write("# Canonical MEG channel labels for dRSA MAT exports\n")
        handle.write("# One channel label per line. Order is export order.\n")
        handle.write("# Auto-generated from the first successful run.\n")
        for ch_name in current_meg_channels:
            handle.write(f"{ch_name}\n")
    return list(current_meg_channels), True


# %%
# -------------------------------------------------------------------------
# Identity-metadata helper
# -------------------------------------------------------------------------
def identities_from_metadata(epochs):
    """
    Recover identities from `epochs.metadata` when available and complete.

    Parameters
    ----------
    epochs : mne.Epochs
        Step-04 epochs object.

    Returns
    -------
    mapping_rows : list of dict or None
        Mapping rows if metadata contains valid `trial_sequence_identity`;
        otherwise None.
    mapping_summary : dict
        Mapping provenance summary.
    mapping_problems : list of str
        Fatal validation problems if metadata is present but invalid.
    """
    if epochs.metadata is None or "trial_sequence_identity" not in epochs.metadata:
        return None, {"mapping_source": "metadata_absent"}, []

    values = np.asarray(epochs.metadata["trial_sequence_identity"])
    mapping_problems = []
    if len(values) != len(epochs.events):
        mapping_problems.append(
            "metadata trial_sequence_identity length does not match epochs: "
            f"metadata={len(values)}, epochs={len(epochs.events)}"
        )

    mapping_rows = []
    for epoch_index, (value, event_row) in enumerate(zip(values, epochs.events)):
        try:
            identity = int(value)
        except Exception:
            mapping_problems.append(
                f"epoch_index={epoch_index}: metadata identity is not an int: {value!r}"
            )
            continue
        if identity not in identity_values:
            mapping_problems.append(
                f"epoch_index={epoch_index}: metadata identity outside 1..20: {identity}"
            )
            continue
        code = int(event_row[2])
        if code not in code_to_condition:
            mapping_problems.append(
                f"epoch_index={epoch_index}: unexpected event code {code}"
            )
            continue
        mapping_rows.append(
            {
                "epoch_index": epoch_index,
                "concat_sample": int(event_row[0]),
                "condition_code": code,
                "condition_name": code_to_condition[code],
                "source_fif": "",
                "raw_sample": None,
                "identity": identity,
                "identity_sample_raw": None,
                "identity_source": "epochs_metadata",
            }
        )

    mapping_summary = {
        "mapping_source": "epochs_metadata",
        "n_mapping_rows": len(mapping_rows),
        "n_mapping_problems": len(mapping_problems),
    }
    return mapping_rows, mapping_summary, mapping_problems


# %%
# -------------------------------------------------------------------------
# Load epochs and derive identity mapping
# -------------------------------------------------------------------------
if not op.exists(epochs_fif):
    raise FileNotFoundError(f"Required epochs file not found: {epochs_fif}")

epochs = mne.read_epochs(epochs_fif, preload=False, verbose=True)
if epochs.event_id != event_id:
    raise RuntimeError(
        "Epoch event_id does not match expected dRSA conditions.\n"
        f"Expected: {event_id}\nObserved: {epochs.event_id}"
    )

meg_picks_all = mne.pick_types(
    epochs.info,
    meg=True,
    eeg=False,
    eog=False,
    stim=False,
    misc=False,
    exclude=[],
)
if len(meg_picks_all) == 0:
    raise RuntimeError("No MEG channels found in epochs; cannot export neural tensors.")
marked_bad_channels = list(dict.fromkeys(str(ch) for ch in epochs.info.get("bads", [])))
meg_channel_names_all = [epochs.ch_names[idx] for idx in meg_picks_all]

# Canonical channel labels are fixed across subjects/runs to guarantee stable
# tensor shape. Channels unavailable or marked bad in a given file are exported
# as NaN later during per-condition tensor assembly.
canonical_meg_channel_names, canonical_file_created = (
    load_or_initialize_canonical_meg_channels(
        canonical_path=canonical_channels_txt,
        current_meg_channels=meg_channel_names_all,
    )
)
present_canonical_meg_channels = [
    ch for ch in canonical_meg_channel_names if ch in meg_channel_names_all
]
missing_canonical_meg_channels = [
    ch for ch in canonical_meg_channel_names if ch not in meg_channel_names_all
]
unexpected_meg_channels = [
    ch for ch in meg_channel_names_all if ch not in canonical_meg_channel_names
]
bad_meg_channels_marked = [
    ch for ch in canonical_meg_channel_names if ch in marked_bad_channels
]
good_meg_channels_present = [
    ch for ch in present_canonical_meg_channels if ch not in marked_bad_channels
]
if len(present_canonical_meg_channels) == 0:
    raise RuntimeError(
        "No canonical MEG channels were found in this epochs file.\n"
        f"canonical_channels_txt={canonical_channels_txt}\n"
        f"n_canonical_meg_channels={len(canonical_meg_channel_names)}\n"
        f"n_meg_channels_available={len(meg_channel_names_all)}"
    )
if len(good_meg_channels_present) == 0:
    raise RuntimeError(
        "All canonical MEG channels present in this file are marked bad.\n"
        f"marked_bad_channels={marked_bad_channels}\n"
        "Cannot compute identity averages without at least one valid MEG channel."
    )

meg_selected_indices = [epochs.ch_names.index(ch) for ch in present_canonical_meg_channels]
channel_type_counts = Counter(mne.channel_type(epochs.info, idx) for idx in meg_selected_indices)

log_to_console_and_report(
    "dRSA export setup complete.\n"
    f"project_root={project_root}\n"
    f"subj_str={subj_str}\n"
    f"ses_str={ses_str}\n"
    f"epochs_fif={epochs_fif}\n"
    f"event_id={event_id}\n"
    f"input_epochs_sfreq={epochs.info['sfreq']}\n"
    f"input_epoch_window=[{epochs.tmin}, {epochs.tmax}] sec\n"
    f"export_window=[{export_tmin}, {export_tmax}] sec\n"
    f"export_sfreq={export_sfreq}\n"
    "explicit_low_pass_applied=False\n"
    f"n_epochs={len(epochs)}\n"
    f"n_channels_total={epochs.info['nchan']}\n"
    f"n_meg_channels_available={len(meg_channel_names_all)}\n"
    f"canonical_channels_txt={canonical_channels_txt}\n"
    f"canonical_channel_file_created={canonical_file_created}\n"
    f"n_canonical_meg_channels={len(canonical_meg_channel_names)}\n"
    f"n_canonical_meg_channels_present_in_epochs={len(present_canonical_meg_channels)}\n"
    f"n_canonical_meg_channels_missing_from_epochs={len(missing_canonical_meg_channels)}\n"
    f"canonical_meg_channels_missing_from_epochs={missing_canonical_meg_channels}\n"
    f"n_unexpected_meg_channels_not_in_canonical={len(unexpected_meg_channels)}\n"
    f"unexpected_meg_channels_not_in_canonical={unexpected_meg_channels}\n"
    f"n_good_meg_channels_present={len(good_meg_channels_present)}\n"
    f"n_bad_meg_channels_marked={len(bad_meg_channels_marked)}\n"
    f"bad_meg_channels_marked={bad_meg_channels_marked}\n"
    f"meg_channel_type_counts={dict(sorted(channel_type_counts.items()))}\n"
    f"n_marked_bad_channels={len(marked_bad_channels)}\n"
    f"marked_bad_channels={marked_bad_channels}\n"
    f"drsa_output_dir={drsa_output_dir}\n"
    f"figure_output_dir={figure_output_dir}\n"
    f"out_report={out_report}",
    title="Setup summary",
    section="Preparation",
)

mapping_rows, mapping_summary, mapping_problems = identities_from_metadata(epochs)
if mapping_rows is None:
    problem_text = (
        "Epoch metadata is missing trial_sequence_identity; MAT export will not run.\n"
        "Re-run 04_epoching.py to regenerate the epoched FIF with identity metadata."
    )
    save_report_before_error(problem_text)
    raise RuntimeError(problem_text)

if mapping_problems:
    problem_text = (
        "Identity mapping failed; MAT export will not run.\n"
        f"n_mapping_problems={len(mapping_problems)}\n"
        + "\n".join(f"  - {line}" for line in mapping_problems[:50])
    )
    if len(mapping_problems) > 50:
        problem_text += f"\n  ... truncated {len(mapping_problems) - 50} additional problems"
    save_report_before_error(problem_text)
    raise RuntimeError(problem_text)

if len(mapping_rows) != len(epochs.events):
    message = (
        "Identity mapping row count mismatch after validation.\n"
        f"n_mapping_rows={len(mapping_rows)}\n"
        f"n_epochs={len(epochs.events)}"
    )
    save_report_before_error(message)
    raise RuntimeError(message)

mapping_counts = {
    condition_name: {
        identity: int(
            sum(
                row["condition_name"] == condition_name
                and int(row["identity"]) == identity
                for row in mapping_rows
            )
        )
        for identity in identity_values
    }
    for condition_name in event_id
}

log_to_console_and_report(
    "Identity mapping completed successfully.\n"
    f"mapping_source={mapping_summary['mapping_source']}\n"
    f"n_mapping_rows={len(mapping_rows)}\n"
    f"n_mapping_problems=0\n"
    f"mapping_summary={mapping_summary}\n"
    "missing_or_ambiguous_identity_mappings=None",
    title="Identity mapping summary",
    section="Identity mapping",
    tags=("identity", "events", "qc"),
)

# %%
# -------------------------------------------------------------------------
# Identity count plots and report table
# -------------------------------------------------------------------------
html_rows = [
    "<table>",
    "<thead><tr><th>condition</th>"
    + "".join(f"<th>{identity}</th>" for identity in identity_values)
    + "<th>total</th></tr></thead>",
    "<tbody>",
]
for condition_name in event_id:
    counts = mapping_counts[condition_name]
    total = sum(counts.values())
    html_rows.append(
        "<tr><td>{cond}</td>{cells}<td>{total}</td></tr>".format(
            cond=html.escape(condition_name),
            cells="".join(f"<td>{counts[identity]}</td>" for identity in identity_values),
            total=total,
        )
    )
html_rows.extend(["</tbody>", "</table>"])
report.add_html(
    "\n".join(html_rows),
    title="Trial counts per identity",
    section="Identity mapping",
    tags=("identity", "counts"),
)

for condition_name in event_id:
    counts = mapping_counts[condition_name]
    fig_counts, ax = plt.subplots(figsize=(10, 4))
    ax.bar(identity_values, [counts[identity] for identity in identity_values])
    ax.set_xlabel("trial_sequence_identity")
    ax.set_ylabel("surviving epochs")
    ax.set_title(f"{condition_name}: trial counts per identity")
    ax.set_xticks(identity_values)
    ax.set_xlim(0.25, 20.75)
    ax.grid(axis="y", alpha=0.25)
    report.add_figure(
        fig=fig_counts,
        title=f"Counts by identity ({condition_name})",
        section="Identity mapping",
        tags=("identity", "counts", "qc"),
    )
    savefig_logged(
        fig_counts,
        f"06_export_neural_for_drsa_counts_{condition_name}.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig_counts)

# %%
# -------------------------------------------------------------------------
# Before/after downsample example trace figure
# -------------------------------------------------------------------------
example_condition = next(name for name in event_id if len(epochs[name]) > 0)
example_channels = good_meg_channels_present[: min(3, len(good_meg_channels_present))]
example_epochs = epochs[example_condition][:1].copy()
example_epochs.load_data()
example_epochs.pick(example_channels)
example_epochs.crop(tmin=export_tmin, tmax=export_tmax)
example_before = example_epochs.get_data(copy=False)[0]
example_times_before = example_epochs.times.copy()
example_epochs_ds = example_epochs.copy().resample(export_sfreq, npad="auto")
example_after = example_epochs_ds.get_data(copy=False)[0]
example_times_after = example_epochs_ds.times.copy()

fig_ds, axes = plt.subplots(
    len(example_channels),
    1,
    figsize=(10, 2.2 * len(example_channels)),
    sharex=True,
)
if len(example_channels) == 1:
    axes = [axes]
for ch_idx, (axis, ch_name) in enumerate(zip(axes, example_channels)):
    axis.plot(
        example_times_before,
        example_before[ch_idx],
        color="0.65",
        linewidth=0.8,
        label=f"before ({epochs.info['sfreq']:.0f} Hz)",
    )
    axis.plot(
        example_times_after,
        example_after[ch_idx],
        color="tab:blue",
        linewidth=1.0,
        label=f"after ({export_sfreq:.0f} Hz)",
    )
    axis.set_ylabel(ch_name)
    axis.grid(axis="x", alpha=0.2)
axes[-1].set_xlabel("time from condition onset (s)")
axes[0].legend(loc="upper right")
fig_ds.suptitle(
    f"Before/after downsample example: {example_condition}, first surviving epoch"
)
fig_ds.tight_layout()

report.add_figure(
    fig=fig_ds,
    title="Before/after downsample example traces",
    caption=(
        f"Condition={example_condition}; channels={example_channels}; "
        f"explicit_low_pass_applied=False."
    ),
    section="Signal checks",
    tags=("resample", "trace", "qc"),
)
savefig_logged(
    fig_ds,
    "06_export_neural_for_drsa_downsample_example.png",
    dpi=150,
    bbox_inches="tight",
)
plt.close(fig_ds)

log_to_console_and_report(
    "Generated downsampling trace example.\n"
    f"example_condition={example_condition}\n"
    f"example_channels={example_channels}\n"
    f"before_sfreq={epochs.info['sfreq']}\n"
    f"before_n_times={len(example_times_before)}\n"
    f"after_sfreq={example_epochs_ds.info['sfreq']}\n"
    f"after_n_times={len(example_times_after)}",
    title="Downsample example summary",
    section="Signal checks",
    tags=("resample", "trace", "qc"),
)

# %%
# -------------------------------------------------------------------------
# Per-condition export: crop, downsample, average by identity, save MAT
# -------------------------------------------------------------------------
export_summaries = []
for condition_name, condition_code in event_id.items():
    condition_rows = [
        row for row in mapping_rows if row["condition_name"] == condition_name
    ]
    condition_epochs = epochs[condition_name].copy()

    observed_samples = [int(sample) for sample in condition_epochs.events[:, 0]]
    expected_samples = [int(row["concat_sample"]) for row in condition_rows]
    if observed_samples != expected_samples:
        message = (
            "Condition epoch order mismatch before export.\n"
            f"condition={condition_name}\n"
            f"observed_first={observed_samples[:10]}\n"
            f"expected_first={expected_samples[:10]}"
        )
        save_report_before_error(message)
        raise RuntimeError(message)

    identity_by_epoch = np.asarray(
        [int(row["identity"]) for row in condition_rows],
        dtype=int,
    )
    missing_identities = [
        identity
        for identity in identity_values
        if int(np.sum(identity_by_epoch == identity)) == 0
    ]
    if missing_identities:
        message = (
            "Cannot export exactly 20 identity averages because some identities "
            "have no surviving epochs.\n"
            f"condition={condition_name}\n"
            f"missing_identities={missing_identities}"
        )
        save_report_before_error(message)
        raise RuntimeError(message)

    # Load only this condition, then restrict to currently available canonical
    # channels and crop before downsampling. No explicit `filter` call is made
    # in this export step.
    condition_epochs.load_data()
    condition_epochs.pick(present_canonical_meg_channels)
    # Guardrail: after picking, in-file channel order must match canonical
    # order restricted to channels physically present in this FIF.
    if list(condition_epochs.ch_names) != list(present_canonical_meg_channels):
        raise RuntimeError(
            "Condition-specific channel order mismatch before export.\n"
            f"condition={condition_name}\n"
            f"expected_first10={present_canonical_meg_channels[:10]}\n"
            f"observed_first10={condition_epochs.ch_names[:10]}"
        )
    condition_epochs.crop(tmin=export_tmin, tmax=export_tmax)
    condition_epochs.resample(export_sfreq, npad="auto")
    data = condition_epochs.get_data(copy=False)

    # Build per-identity averages for channels with available samples.
    tensor_present = np.empty((len(identity_values), data.shape[1], data.shape[2]), dtype=float)
    counts_by_identity = {}
    for identity in identity_values:
        epoch_indices = np.flatnonzero(identity_by_epoch == identity)
        counts_by_identity[identity] = int(len(epoch_indices))
        tensor_present[identity - 1, :, :] = data[epoch_indices, :, :].mean(axis=0)

    # Expand to canonical channel axis and keep bad/missing channels as NaN.
    tensor = np.full(
        (len(identity_values), len(canonical_meg_channel_names), data.shape[2]),
        np.nan,
        dtype=float,
    )
    present_idx_by_name = {name: idx for idx, name in enumerate(condition_epochs.ch_names)}
    n_channels_with_signal = 0
    nan_channels_this_condition = []
    filled_channels_this_condition = []
    for canonical_idx, ch_name in enumerate(canonical_meg_channel_names):
        in_data = ch_name in present_idx_by_name
        marked_bad = ch_name in marked_bad_channels
        if in_data and not marked_bad:
            tensor[:, canonical_idx, :] = tensor_present[:, present_idx_by_name[ch_name], :]
            n_channels_with_signal += 1
            filled_channels_this_condition.append(ch_name)
        else:
            nan_channels_this_condition.append(ch_name)

    mat_path = op.join(
        drsa_output_dir,
        (
            f"sub-{subj_str}_ses-{ses_str}_task-Attention_"
            f"condition-{condition_name}_desc-identityAveraged_neural.mat"
        ),
    )
    channel_labels = np.asarray(canonical_meg_channel_names, dtype=object).reshape(-1, 1)
    savemat(
        mat_path,
        {"neural": tensor, "channel_labels": channel_labels},
        do_compression=True,
    )
    print(f"Saved dRSA MAT tensor to: {mat_path}")

    export_summaries.append(
        {
            "condition": condition_name,
            "condition_code": condition_code,
            "mat_path": mat_path,
            "tensor_shape": tuple(int(value) for value in tensor.shape),
            "n_input_epochs": int(len(condition_epochs)),
            "counts_by_identity": counts_by_identity,
            "sfreq": float(condition_epochs.info["sfreq"]),
            "n_times": int(len(condition_epochs.times)),
            "time_start": float(condition_epochs.times[0]),
            "time_end": float(condition_epochs.times[-1]),
            "n_meg_channels": int(len(canonical_meg_channel_names)),
            "n_meg_channels_present_in_epochs": int(len(condition_epochs.ch_names)),
            "n_meg_channels_with_signal": int(n_channels_with_signal),
            "n_meg_channels_nan": int(len(nan_channels_this_condition)),
            "first_channels_with_signal": filled_channels_this_condition[:10],
            "first_channels_nan": nan_channels_this_condition[:10],
            "mat_variables": ["neural", "channel_labels"],
        }
    )

    del condition_epochs
    del data
    del tensor_present
    del tensor

# %%
# -------------------------------------------------------------------------
# Final output summaries and report save
# -------------------------------------------------------------------------
summary_lines = ["dRSA neural export completed successfully."]
summary_lines.append(
    "Export channel set: MEG only (canonical fixed list); "
    f"n_meg_channels={len(canonical_meg_channel_names)}; "
    f"channel_type_counts={dict(sorted(channel_type_counts.items()))}"
)
summary_lines.append(f"Marked bad channels (epochs): {marked_bad_channels}")
summary_lines.append(f"Canonical channel file: {canonical_channels_txt}")
summary_lines.append(f"Canonical channels missing from this epochs file: {missing_canonical_meg_channels}")
summary_lines.append(f"Unexpected MEG channels not present in canonical list: {unexpected_meg_channels}")
summary_lines.append(
    "Bad/missing channels are preserved in channel_labels and exported as NaN in neural."
)
for item in export_summaries:
    summary_lines.extend(
        [
            "",
            f"Condition: {item['condition']} ({item['condition_code']})",
            f"  mat_path={item['mat_path']}",
            f"  mat_variables={item['mat_variables']}",
            f"  tensor_shape={item['tensor_shape']}",
            f"  n_input_epochs={item['n_input_epochs']}",
            f"  counts_by_identity={item['counts_by_identity']}",
            f"  sfreq={item['sfreq']}",
            f"  n_times={item['n_times']}",
            f"  time_range=[{item['time_start']}, {item['time_end']}] sec",
            f"  n_meg_channels={item['n_meg_channels']}",
            f"  n_meg_channels_present_in_epochs={item['n_meg_channels_present_in_epochs']}",
            f"  n_meg_channels_with_signal={item['n_meg_channels_with_signal']}",
            f"  n_meg_channels_nan={item['n_meg_channels_nan']}",
            f"  first_channels_with_signal={item['first_channels_with_signal']}",
            f"  first_channels_nan={item['first_channels_nan']}",
        ]
    )

log_to_console_and_report(
    "\n".join(summary_lines),
    title="Final export summary",
    section="Outputs",
    tags=("outputs", "drsa"),
)

report.save(out_report, overwrite=True, open_browser=False)
print(f"Saved HTML report to: {out_report}")
