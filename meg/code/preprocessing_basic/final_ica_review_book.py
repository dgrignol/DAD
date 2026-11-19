#!/usr/bin/env python3
"""Build a compact, offline component navigator for one FINAL ICA model.

Usage inside 03_C_ICA.py, after diagnostics are saved
---------------------------------------------------
    coverage = write_component_book(payload, diagnostic_images, output, raw_root)

Inputs: the current final model report and ONLY the image paths produced for
that model in this invocation. Output: one self-contained HTML, with component
selection, keyboard navigation, all maps/properties and per-run native context.
No source data, preliminary selections, component decisions or tuning are read
or changed. Navigation and display do not persist review or accept exclusions.
Component identity is the participant-level final model SHA-256; preliminary
per-run annotation-aid indices have no correspondence to this decomposition.

Data flow: saved current-model PNGs + fit/sample/convergence evidence -> embedded
offline navigator. Coverage is explicit; missing properties block final review
acceptance instead of being filled by stale PNGs left by an earlier model.
"""

import base64
import html
import json
from pathlib import Path

from preprocessing_stage_utils import write_text


# -----------------------------------------------------------------------------
# Embed current diagnostics with explicit model identity and coverage
# -----------------------------------------------------------------------------

def review_guide_html(payload):
    """Explain the human component-review procedure inside the offline book."""
    scope = payload.get("diagnostic_scope", "legacy_first_fit_stretch")
    warning = ('<p class="alert"><strong>Legacy property panels cover only the first retained fit stretch. '
        'Use their full-fit commands, or a report-only regeneration, for the corrected display. The fitted model used all recorded fit samples.</strong></p>'
        if scope != "full_fit" else '')
    return warning + '''<details id="guide" open><summary><strong>How to review this participant</strong></summary>
<ol><li>Check model identity, fitted/mapped sample counts, dimension and convergence. An iteration-limit warning needs reconciliation before final decisions.</li>
<li>Inspect EVERY final component: mag/grad map, spectrum, arbitrary-segment image/mean, variance and temporal source traces. Keep a table of IC index, keep/exclude/defer, hypothesis, confidence and evidence.</li>
<li>For a suspected artifact, open the native source browser and frontal/MISC context. Check several examples and ordinary controls across acquisitions. The small native heatmap is the first labelled ocular example, not a systematic BAD-enrichment analysis. Labels and exact window bounds are below it. Its within-window IC z-scores cannot rank component strengths.</li>
<li>Inspect BOTH the default native continuous traces and the task-locked trial images/condition means for EVERY component. Continuous examples retain BAD periods; task views use accepted Stage 2/3 eligibility and corrected motion-onset samples with the configured -0.5/+2.667 s window and BAD rejection. These transient diagnostics precede cleaning; they are not production epochs. Task locking alone distinguishes neither brain nor artifact: ocular activity can also be task-locked. Inspect spatial, spectral, temporal and MISC evidence together.</li>
<li>Ordinary property panels are NOT task-locked: MNE divides the protected 1-40 Hz/200 Hz fit copies into arbitrary 2 s segments. Their mean and confidence band are not a task ERF. Fit data omit all BAD categories, including routine ocular periods. Missing artifacts here do not establish that a component is neural.</li>
<li>Discuss mixed/uncertain components explicitly; defer when evidence is insufficient. Never transfer preliminary per-run frontal ICA indices. A browser click affects only the temporary in-memory object; these commands save no exclusions or review receipts.</li>
<li>After explicit participant-level approval, record the exact final-model exclusion list with review_cleaning.py --stage 4C. An explicitly reviewed empty list is valid. Application, native before/after QC, P10's slow-contamination decision and cleaned-output acceptance follow separately; Stage 5 remains later.</li></ol>
<p><strong>Time axes:</strong> property panels show 0-1.995 s within arbitrary fit segments; fit-source browsing uses compressed fit time; native examples use acquisition seconds; TASK plots use seconds from accepted motion onset (codes 71/81/91), not uncorrected digital trial start.</p>
<p><strong>Commands:</strong> run local commands in the existing MNE environment. Cluster commands require the dad-meg environment and a GUI for interactive browsing; add --save /tmp/name.png for a static diagnostic. Every command pins this model SHA-256. Reproducing full-fit properties reads/rebuilds display copies but never refits ICA.</p>
<p><strong>Precise references and interpretation:</strong>
<a href="https://mne.tools/1.10/generated/mne.preprocessing.ICA.html#mne.preprocessing.ICA.plot_properties">MNE plot_properties API</a> defines the map/image/mean/spectrum/variance panels and confidence bands.
<a href="https://mne.tools/1.10/auto_tutorials/preprocessing/40_artifact_correction_ica.html#selecting-ica-components-manually">MNE manual component selection</a> demonstrates combining time courses, maps and spectra; variance extremes are investigation leads, not automatic exclusions.
<a href="https://labeling.ucsd.edu/tutorial/labels">ICLabel: Telling Components Apart</a>, especially Brain, Muscle, Eye and Heart, gives morphology examples. Its EEG scalp/dipole residual-variance thresholds must not be imported into this MEG/SSS workflow. An alpha peak, a dipolar-looking map or absence of an ERP alone is insufficient to classify an IC.</p></details>'''


def write_component_book(payload, images, output, raw_root):
    """Return diagnostic completeness while writing one participant's review book."""
    images = [Path(path) for path in images]
    n_components = int(payload["n_components_fitted"])
    groups = {i: [p for p in images if f"_ICA{i:03d}_properties_" in p.name] for i in range(n_components)}
    missing = [i for i, paths in groups.items() if not paths or not all(p.is_file() for p in paths)]
    continuous = {i: [p for p in images if f"_ICA{i:03d}_continuous_" in p.name] for i in range(n_components)}
    task = {i: [p for p in images if f"_ICA{i:03d}_task_" in p.name] for i in range(n_components)}
    missing_continuous = [i for i, paths in continuous.items() if not paths or not all(p.is_file() for p in paths)]
    missing_task = [i for i, paths in task.items() if not paths or not all(p.is_file() for p in paths)]
    native_required = payload.get('review_interface_version', 1) >= 3
    summary = {"complete": not missing and bool(payload.get("component_properties_requested")),
        "missing_properties": missing, "n_components": n_components,
        "model_sha256": payload["prefit_model_sha256"]}
    if native_required:
        summary.update(missing_continuous=missing_continuous, missing_task=missing_task,
            complete=summary['complete'] and not missing_continuous and not missing_task)

    def figure(path):
        """Use native PNGs so offline review needs no server or remote assets."""
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        recipe = payload.get("plot_recipes", {}).get(path.name, {})
        caption = '<p>' + html.escape(recipe.get("description", "Saved current-model diagnostic.")) + '</p>'
        commands = recipe.get("commands", {})
        if commands:
            caption += '<strong>Reproduce locally / open interactively</strong><pre class="command"><code>' + html.escape(commands["local"]) + '</code></pre>'
            caption += '<details><summary>Cluster command</summary><pre class="command"><code>' + html.escape(commands["cluster"]) + '</code></pre></details>'
        for key, label in (("full_commands", "Corrected full-fit properties (same model; no refit)"),
                           ("task_commands", "Inspect this component locked to accepted motion onset")):
            if recipe.get(key):
                caption += '<details><summary>' + label + '</summary><pre class="command"><code>' + html.escape(recipe[key]["local"]) + '</code></pre></details>'
        if recipe.get('native_browser_commands'):
            caption += '<details><summary>Scroll the full native acquisitions for this component</summary>'
            for command in recipe['native_browser_commands']:
                caption += '<p>' + html.escape(command['source']) + '</p><pre class="command"><code>' + html.escape(command['local']) + '</code></pre>'
                caption += '<details><summary>Cluster</summary><pre class="command"><code>' + html.escape(command['cluster']) + '</code></pre></details>'
            caption += '</details>'
        if recipe.get("metadata"):
            caption += '<details><summary>Exact scope / BAD labels / native window selection</summary><pre>' + html.escape(json.dumps(recipe["metadata"], indent=2, default=str, sort_keys=True)) + '</pre></details>'
        return f'<figure><img loading="lazy" src="data:image/png;base64,{encoded}" alt="{html.escape(path.stem)}"><figcaption>{html.escape(path.stem)}{caption}</figcaption></figure>'

    common = [p for p in images if not any(tag in p.name for tag in ('_properties_', '_continuous_', '_task_'))]
    sections = []
    for i, paths in groups.items():
        sections.append(f'<section class="component" id="component-{i}" hidden><h2>ICA{i:03d}</h2>' +
            ("".join(figure(p) for p in paths) if paths else '<p class="alert">Missing current-model property diagnostics: review incomplete.</p>') +
            ''.join(figure(p) for p in continuous[i] + task[i]) +
            ('<p class="alert">Missing continuous or task-locked diagnostic: review incomplete.</p>'
                if native_required and (i in missing_continuous or i in missing_task) else '') + '</section>')
    model = html.escape(payload["prefit_model_sha256"])
    options = ''.join(f'<option value="{i}">ICA{i:03d}</option>' for i in range(n_components))
    evidence = {key: payload.get(key) for key in ("fit_execution", "fit_input_summaries", "fingerprint", "optimizer_restart")}
    fit = payload.get("fit_execution") or {}
    fit_summary = html.escape(f"Fit samples: {fit.get('n_samples_fitted', 'not reported')}; "
        f"mapped samples: {fit.get('mapped_fit_samples', 'not reported')}; "
        f"iterations: {fit.get('n_iterations', 'not reported')} / {fit.get('max_iter', 'not reported')}.")
    convergence = ('<p class="alert"><strong>Iteration limit reached: reconcile convergence before component decisions.</strong></p>'
        if fit.get("iteration_limit_reached") else '')
    restart = payload.get("optimizer_restart")
    if restart:
        convergence += ('<p><strong>Numerically reconciled deterministic restart.</strong> '
            'The original seed-97 solution cycled. Its full PCA basis was reinitialized from its update-curvature eigenvectors; '
            'standard parallel/logcosh FastICA and its original tolerance were retained. '
            f'Original phase: {int(restart["original_phase_iterations"])} iterations; '
            f'restart phase: {int(restart["restart_phase_iterations"])} iterations. '
            'This was not a fresh random initialization. All component decisions require this current model; '
            'numerical convergence grants no scientific acceptance. Full history is in Fit evidence.</p>')
    text = f'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(payload['subject_bids_id'])} final ICA review</title>
<style>body{{font:16px system-ui;margin:24px auto;max-width:1250px;padding:0 18px;color:#162333;background:#f8fafc}}
nav{{position:sticky;top:0;background:#e6eef7;padding:12px;z-index:2}}button,select{{font:inherit;margin:4px;padding:8px}}
img{{max-width:100%;height:auto}}figure{{margin:15px 0;background:white}}figcaption{{font-size:12px;overflow-wrap:anywhere}}
code{{overflow-wrap:anywhere}}pre{{white-space:pre-wrap;overflow-wrap:anywhere}}.command{{background:#edf2f7;padding:12px;font-size:12px}}.alert{{color:#922}}[hidden]{{display:none!important}}</style>
<h1>{html.escape(payload['subject_bids_id'])}: final ICA component review</h1>
<p>Model SHA-256: <code>{model}</code>. {n_components} components. Exclusion decisions remain explicit and model-specific.</p>
<p>{fit_summary}</p>{convergence}
<p>Inspect every map, property plot and source trace, then use native run context. Preliminary frontal ICA indices must never be transferred here.
Fit data omit every operative BAD category. These 1-40 Hz plots cannot resolve the open P10 approximately 0.33 Hz concern; use native before/after application QC.</p>
<p class="alert">Diagnostic coverage: {'complete' if summary['complete'] else 'INCOMPLETE'}. Execution and navigation do not establish scientific acceptance.</p>
{review_guide_html(payload)}
<nav><button id="previous" type="button">Previous</button><label>Component <select id="picker">{options}</select></label>
<button id="next" type="button">Next</button><a href="#guide">Review guide</a> | <a href="#context">Maps, traces and run context</a> | <a href="#evidence">Fit evidence</a>
<span> Keyboard: j/k or left/right arrows.</span></nav>
{''.join(sections)}<details id="context" open><summary>All component maps, fit source traces and native run context</summary>{''.join(figure(p) for p in common)}</details>
<details id="evidence"><summary>Exact fit preparation, sample accounting and convergence evidence</summary><pre>{html.escape(json.dumps(evidence, indent=2, default=str, sort_keys=True))}</pre></details>
<script>const picker=document.getElementById('picker');function show(i){{i=Math.max(0,Math.min({n_components-1},Number(i)));picker.value=String(i);
document.querySelectorAll('.component').forEach(e=>e.hidden=e.id!=='component-'+i);}}
picker.addEventListener('change',()=>show(picker.value));document.getElementById('previous').onclick=()=>show(Number(picker.value)-1);
document.getElementById('next').onclick=()=>show(Number(picker.value)+1);document.addEventListener('keydown',e=>{{
if(['SELECT','INPUT','TEXTAREA'].includes(e.target.tagName))return;
if(['j','ArrowRight','k','ArrowLeft'].includes(e.key)){{e.preventDefault();show(Number(picker.value)+(['j','ArrowRight'].includes(e.key)?1:-1));}}}});show(0);</script></html>'''
    write_text(output, text, raw_root=raw_root, overwrite=True)
    return summary
